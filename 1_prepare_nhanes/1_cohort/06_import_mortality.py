# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/1_cohort/06_import_mortality.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Reads the NCHS public-use linked mortality files (deaths
#               through 31 December 2019) for the ten cycles and keeps
#               linkage-eligible records: vital status, underlying cause group,
#               and months of follow-up from examination.
#
# METHODS       "Experimental design": "each linked to National Death Index
#               records through 31 December 2019". The endpoint itself is
#               defined in 3_results/data.py (principal_endpoint): death with a
#               recorded cause other than accident (ucod_leading 4).
#
# PSEUDOCODE
#   read each fixed-width file at the published column positions
#   keep eligstat == 1 (eligible for linkage)
#   keep mortstat, ucod_leading, permth_exm (months from examination), ...
#   Participants with no row here have no linkage; 3_results drops them
#   (follow-up must be > 0): 12 people aged 52-79 with complete panels.
# =============================================================================
"""Reconstruct the lossless NHANES 2019 linked-mortality dataframe.

This is a transparent Python translation of the corrected historical
``NHANES_mort_to_csv_lossless.R`` parser. It reads the ten public-use fixed-width
files, retains mortality-eligible records, and preserves causes, MCOD flags and
both follow-up clocks. Endpoint definitions belong downstream.
"""

from hashlib import sha256
import os
from pathlib import Path

import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
OUT = ROOT / "outputs/01_imported"
PROVENANCE = ROOT / "outputs/provenance"  # public_code: records go with the outputs
COLSPECS = [(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21),
            (42, 45), (45, 48)]
NAMES = ["SEQN", "eligstat", "mortstat", "ucod_leading", "diabetes_mcod",
         "hyperten_mcod", "permth_int", "permth_exm"]
UCOD = {1: "Diseases of heart", 2: "Malignant neoplasms",
        3: "Chronic lower respiratory diseases",
        4: "Accidents (unintentional injuries)", 5: "Cerebrovascular diseases",
        6: "Alzheimer's disease", 7: "Diabetes mellitus",
        8: "Influenza and pneumonia",
        9: "Nephritis, nephrotic syndrome and nephrosis",
        10: "All other causes (residual)"}


def mortality_root() -> Path:
    value = os.environ.get("NHANES_MORTALITY_ROOT")
    if not value:
        raise RuntimeError("NHANES_MORTALITY_ROOT is not set; see RUN_ORDER.md")
    root = Path(value).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)
    return root


def fingerprint(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest().upper()


def cause_label(row) -> str:
    if pd.notna(row.ucod_leading):
        code = int(row.ucod_leading)
        if code not in UCOD:
            raise AssertionError(f"unexpected UCOD_LEADING value: {code}")
        return UCOD[code]
    return "Deceased, cause not recorded" if row.mortstat == 1 else "Alive"


def main() -> None:
    sources = pd.read_csv(ROOT / "config/mortality_sources.csv")
    source_root = mortality_root()
    frames, manifest = [], []
    for row in sources.itertuples(index=False):
        path = source_root / row.filename
        frame = pd.read_fwf(path, colspecs=COLSPECS, names=NAMES,
                            na_values=["", "."])
        if frame.SEQN.isna().any() or frame.SEQN.duplicated().any():
            raise AssertionError(f"invalid SEQN in mortality cycle {row.cycle}")
        frame["cycle"] = row.cycle
        frames.append(frame)
        manifest.append({"cycle": row.cycle, "filename": row.filename,
                         "bytes": path.stat().st_size, "sha256": fingerprint(path),
                         "raw_rows": len(frame)})

    raw = pd.concat(frames, ignore_index=True)
    mortality = raw.loc[raw.eligstat.eq(1)].copy()
    mortality["ucod_label"] = mortality.apply(cause_label, axis=1)
    mortality = mortality[["SEQN", "mortstat", "ucod_leading", "ucod_label",
                           "diabetes_mcod", "hyperten_mcod", "permth_int",
                           "permth_exm", "cycle"]].reset_index(drop=True)
    if mortality.SEQN.duplicated().any():
        raise AssertionError("duplicate SEQN across mortality cycles")
    mortality["SEQN"] = mortality.SEQN.astype(int)
    for column in ["mortstat", "ucod_leading", "diabetes_mcod", "hyperten_mcod",
                   "permth_int", "permth_exm"]:
        mortality[column] = mortality[column].astype("Int64")

    deaths = (mortality.loc[mortality.mortstat.eq(1)]
              .groupby(["ucod_leading", "ucod_label"]).size()
              .reset_index(name="deaths").sort_values("ucod_leading"))
    OUT.mkdir(parents=True, exist_ok=True)
    mortality.to_csv(OUT / "mortality.csv", index=False)
    deaths.to_csv(OUT / "deaths_by_cause.csv", index=False)
    pd.DataFrame(manifest).to_csv(PROVENANCE / "mortality_raw_manifest.csv", index=False)
    pd.DataFrame([{"raw_rows": len(raw), "mortality_eligible_rows": len(mortality),
                   "alive": int(mortality.mortstat.eq(0).sum()),
                   "deaths": int(mortality.mortstat.eq(1).sum()),
                   "duplicate_seqn": int(mortality.SEQN.duplicated().sum())}]).to_csv(
        OUT / "mortality_import_summary.csv", index=False)
    print(f"mortality: raw={len(raw):,}; eligible={len(mortality):,}; "
          f"deaths={mortality.mortstat.eq(1).sum():,}")


if __name__ == "__main__":
    main()
