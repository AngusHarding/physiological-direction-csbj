# =============================================================================
# READ ME FIRST: 2_prepare_nhanes_iii/nhanes3_io.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Readers for the NHANES III fixed-width files (column positions
#               from the published SAS layouts), the plausibility rule that sets
#               impossible values to missing, and distribution reports.
#
# METHODS       Data handling for Supp. Section Zero; estimates nothing.
# =============================================================================
"""Layout-driven readers for NHANES III fixed-width files.

SPECIFICATION: see SPECIFICATION.md section 6a, Decision 6.

NHANES III ships three large fixed-width text files and a SAS layout for each.
There is nothing in a .dat file that identifies where one variable ends and the
next begins -- a record is simply ~2,000 characters. Column positions therefore
come from the published layout at runtime and are never hard-coded here, so this
module cannot silently disagree with the layout it claims to implement.

Two failure modes drive the design.

  1. Off-by-one slicing is SILENT. A field read one character out often keeps a
     plausible shape while sitting in the wrong place. Nothing in the data
     objects. Detection is external: range assertions (layer 2) and the
     side-by-side comparison against continuous NHANES (layer 4).

  2. Missing values are FILL CHARACTERS, not blanks. An unavailable numeric is
     filled with 8s or 9s across the whole field width -- 88888 in a five-wide
     field, 8888 in four. Read naively, haemoglobin averages several thousand.
     `read_fixed_width` clears these by construction.

Verified 9 August 2026: the INPUT blocks contain ZERO decimal specifiers and
values carry their own decimal point (haemoglobin '014.7', calcium '09.5'), so
fields parse straight to float. Do NOT apply implied-decimal scaling.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from nhanes3_locations import CONFIG, RAW_NHANES3  # public_code: folders and write guard

# The mortality endpoint is defined once, for the development cohort and for
# NHANES III alike, in endpoints.py (this folder).  # public_code
from endpoints import attach_endpoints  # noqa: E402


@dataclass(frozen=True)
class Field:
    name: str
    start: int          # 1-indexed, inclusive, as printed in the SAS layout
    end: int            # 1-indexed, inclusive
    label: str = ""

    @property
    def width(self) -> int:
        return self.end - self.start + 1

    @property
    def slice(self) -> slice:
        return slice(self.start - 1, self.end)


def load_mortality() -> pd.DataFrame:
    """Mortality-eligible NHANES III participants with the shared endpoints.

    Parses the NCHS 2019 public-use linked mortality file for NHANES III. The
    field positions are those of the continuous-NHANES files except that SEQN
    occupies five characters, not six.

    The NHANES III file releases all ten cause categories, so accidents are
    identifiable throughout (`accident_identifiable=True`) and both the
    all-cause and the non-accidental endpoints are defined for every
    participant. Not every death has a cause code: 13 of the 8,623 eligible
    deaths have none. Those 13 are events under the all-cause endpoint and are
    censored, as deaths of unknown cause, under the non-accidental endpoint.

    Returns one row per eligible participant with `mortstat`, `ucod` (the
    three-character recode as released, used by the cause-specific analysis),
    `ucod_leading` (the same recode as a number), `permth_exm`,
    `followup_years`, and the endpoint columns added by `attach_endpoints`.
    """
    rows = []
    with open(RAW_NHANES3 / "NHANES_III_MORT_2019_PUBLIC.dat", encoding="latin-1") as fh:  # public_code
        for line in fh:
            rows.append({"SEQN": int(line[0:5]), "eligstat": int(line[14:15]),
                         "mortstat": line[15:16].strip(), "ucod": line[16:19].strip(),
                         "permth_exm": line[45:48].strip()})
    m = pd.DataFrame(rows)
    for column in ("mortstat", "permth_exm"):
        m[column] = pd.to_numeric(m[column], errors="coerce")
    m = m[m.eligstat.eq(1)].copy()
    m["ucod_leading"] = pd.to_numeric(m.ucod, errors="coerce")
    attach_endpoints(m, accident_identifiable=True)
    m["followup_years"] = m.permth_exm / 12.0
    return m


def read_layout(sas_path: Path) -> dict[str, Field]:
    """Parse a NHANES III .sas layout into {NAME: Field}.

    The layout has three blocks; we need INPUT (positions) and LABEL
    (descriptions). LENGTH is storage size and is irrelevant to reading text.
    """
    text = Path(sas_path).read_text(encoding="latin-1")
    upper = text.upper()
    i, j = upper.find("INPUT"), upper.find("LABEL")
    if i < 0 or j < 0 or j < i:
        raise AssertionError(f"{sas_path.name}: could not locate INPUT/LABEL blocks")

    positions: dict[str, tuple[int, int]] = {}
    for m in re.finditer(
        r"^\s*([A-Z0-9_]+)\s+(\d+)(?:-(\d+))?\s*$", text[i + 5 : j], re.M | re.I
    ):
        name = m.group(1).upper()
        start = int(m.group(2))
        end = int(m.group(3) or m.group(2))
        positions[name] = (start, end)

    labels: dict[str, str] = {}
    for m in re.finditer(
        r'^\s*([A-Z0-9_]+)\s*=\s*"([^"]*)"', text[j:], re.M | re.I
    ):
        labels[m.group(1).upper()] = m.group(2)

    if not positions:
        raise AssertionError(f"{sas_path.name}: no column positions parsed")
    return {
        n: Field(n, a, b, labels.get(n, "")) for n, (a, b) in positions.items()
    }


def read_fixed_width(
    dat_path: Path, layout: dict[str, Field], wanted: list[str]
) -> pd.DataFrame:
    """Read selected variables from a NHANES III .dat file.

    Fill codes are cleared per field using that field's own width, which is why
    the layout must be consulted rather than assuming a global sentinel.
    """
    missing = [v for v in wanted if v.upper() not in layout]
    if missing:
        raise AssertionError(f"{dat_path.name}: not in layout: {missing}")
    fields = [layout[v.upper()] for v in wanted]

    columns: dict[str, list[str]] = {f.name: [] for f in fields}
    with open(dat_path, "r", encoding="latin-1") as handle:
        for line in handle:
            for f in fields:
                columns[f.name].append(line[f.slice].strip())

    frame = pd.DataFrame(columns)
    for f in fields:
        values = pd.to_numeric(frame[f.name], errors="coerce")
        # NHANES III fills unavailable numerics with the field width in 8s or 9s.
        for fill_char in ("8", "9"):
            sentinel = float(fill_char * f.width)
            values = values.mask(values.eq(sentinel))
        frame[f.name] = values
    return frame


def load_marker_map() -> pd.DataFrame:
    """The marker registry: which NHANES III variable supplies which marker,
    in which file, with the plausible ranges used by the layer-2 assertions."""
    path = CONFIG / "nhanes3_marker_map.csv"  # public_code
    spec = pd.read_csv(path)
    required = {
        "marker", "nhanes3_var", "source_file", "units",
        "median_lo", "median_hi", "hard_lo", "hard_hi", "in_primary_set",
    }
    if not required.issubset(spec.columns):
        raise AssertionError(f"marker map missing columns: {required - set(spec.columns)}")
    if spec.marker.duplicated().any():
        raise AssertionError("duplicate marker in marker map")
    return spec


def apply_plausibility(
    frame: pd.DataFrame, spec: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Set physiologically impossible values to missing, and report what was set.

    This follows the main pipeline: "confirmed implausible values set to missing"
    (SUPPLEMENTARY section S0.6). It is NOT trimming of extremes --
    TRIMMING_ASSESSMENT.md rejected that, because excluded participants were
    mortality-enriched at a death-rate ratio of 2.42, so deletion would remove
    informative pathological observations.

    The distinction is deliberate and the bounds encode it: `hard_lo`/`hard_hi`
    are set at physiological IMPOSSIBILITY, not clinical rarity. A cholesterol of
    702 mg/dL is familial hypercholesterolaemia and is kept; an MCHC of 52 g/dL
    exceeds haemoglobin saturation and cannot be a red cell.
    """
    out = frame.copy()
    rows = []
    for r in spec.itertuples():
        v = r.nhanes3_var
        if v not in out.columns:
            continue
        impossible = (out[v] < r.hard_lo) | (out[v] > r.hard_hi)
        n = int(impossible.sum())
        if n:
            rows.append({
                "marker": r.marker, "nhanes3_var": v, "n_set_missing": n,
                "hard_lo": r.hard_lo, "hard_hi": r.hard_hi,
                "values": ", ".join(str(x) for x in sorted(out.loc[impossible, v].unique())[:6]),
            })
            out.loc[impossible, v] = np.nan
    return out, pd.DataFrame(rows)


def distribution_report(
    frame: pd.DataFrame, spec: pd.DataFrame, sex_col: str = "HSSEX"
) -> pd.DataFrame:
    """Layer 3. Shape and moments per marker, overall and by sex.

    Skewness is what transformation selection minimises, so this table is both a
    quality check and the input to the next stage. D'Agostino K2 tests the
    combined skew/kurtosis departure from normality.
    """
    from scipy import stats

    rows = []
    groups = [("both", frame)]
    if sex_col in frame.columns:
        groups += [("men", frame[frame[sex_col].eq(1)]),
                   ("women", frame[frame[sex_col].eq(2)])]
    for r in spec.itertuples():
        v = r.nhanes3_var
        if v not in frame.columns:
            continue
        for label, g in groups:
            x = g[v].dropna()
            if len(x) < 20:
                continue
            k2, p = stats.normaltest(x)
            rows.append({
                "marker": r.marker, "sex": label, "n": int(len(x)),
                "missing_pct": round(100 * g[v].isna().mean(), 1),
                "mean": round(float(x.mean()), 3),
                "median": round(float(x.median()), 3),
                "sd": round(float(x.std(ddof=1)), 3),
                "min": round(float(x.min()), 3), "max": round(float(x.max()), 3),
                "iqr": round(float(x.quantile(.75) - x.quantile(.25)), 3),
                "skew": round(float(stats.skew(x)), 3),
                "excess_kurtosis": round(float(stats.kurtosis(x)), 3),
                "dagostino_k2": round(float(k2), 1),
                "dagostino_p": float(p),
                "abs_skew_gt_1": bool(abs(stats.skew(x)) > 1),
            })
    return pd.DataFrame(rows)


def range_report(frame: pd.DataFrame, spec: pd.DataFrame) -> pd.DataFrame:
    """Layer 2. Median must sit inside [median_lo, median_hi]; values outside
    [hard_lo, hard_hi] are physiologically impossible and counted.

    Returns one row per marker with a boolean `passes`. The caller decides what
    to do about failures -- per Decision 6 the build halts rather than warns.
    """
    rows = []
    for r in spec.itertuples():
        if r.nhanes3_var not in frame.columns:
            rows.append({"marker": r.marker, "nhanes3_var": r.nhanes3_var,
                         "status": "absent", "passes": False})
            continue
        x = frame[r.nhanes3_var].dropna()
        median = float(x.median()) if len(x) else np.nan
        impossible = int(((x < r.hard_lo) | (x > r.hard_hi)).sum())
        median_ok = bool(r.median_lo <= median <= r.median_hi)
        rows.append({
            "marker": r.marker, "nhanes3_var": r.nhanes3_var, "units": r.units,
            "n": int(len(x)), "coverage_pct": round(100 * len(x) / len(frame), 1),
            "median": round(median, 3),
            "median_lo": r.median_lo, "median_hi": r.median_hi,
            "median_in_range": median_ok,
            "impossible_values": impossible,
            "status": "ok" if (median_ok and impossible == 0) else "FAIL",
            "passes": bool(median_ok and impossible == 0),
        })
    return pd.DataFrame(rows)
