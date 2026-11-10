# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/8_adverse_distance/05_build_deficit_severities.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Severity = the adverse excursion itself (DEAD_ZONE = 0, no
#               cap), so DA is continuous and unbounded. The older bounded rule
#               is kept only for a sensitivity analysis.
#
# METHODS       "Distances": "Adverse distance (DA): the mean excursion beyond
#               the reference median in each marker's prespecified adverse
#               direction, on the transformed scale and in units of the
#               reference interquartile range divided by 1.349."
#               Supp. Section Zero, "Adverse-excursion coordinates":
#               "rij = oj [Tj(xij) - mjs] / (IQRjs / 1.349) ... Each marker's
#               prespecified adverse direction then retains the excursion above
#               the median, below it or on both sides, and DA averages the
#               retained absolute excursions across markers."
#
# NAMES         "FI branch" and "deficits" are this code's historical names for
#               adverse distance (DA) and its per-marker excursions; the paper
#               calls them adverse distance and adverse excursions.
#
# THE SIX STEPS (this folder; DA itself is the mean in 3_results/measures.py)
#   01 lock the reference members (the 1,238 from 7_marker_panel/04)
#   02 reference median and IQR/1.349 per marker and sex
#   03 r = orientation x (transformed value - median) / (IQR/1.349)
#   04 keep only the adverse side: high_only max(r,0); low_only max(-r,0);
#      two_sided |r|
#   05 severity = the excursion itself (no dead zone, no cap)
#   06 put the sign back (+ high, - low) and save the "x" coordinates
#               Methods: DA "is constructed like a laboratory frailty index
#               [27] but is continuous and unbounded."
# =============================================================================
"""Convert adverse excursions to FI deficit severities.

METHODS: see FI_METHODS.md Step 6 (version 2, amended 5 August 2026).

Version 2: severity IS the adverse excursion, d = e. No dead zone, no cap, no
free parameters. Severity is exactly proportional to displacement into the
adverse tail across the whole range.

The version-1 rule d = clip((e - 1) / 2, 0, 1) was flat at both ends: it assigned
zero to any excursion within one robust-scale unit and one to every excursion
beyond three. An index flat at both ends is piecewise-thresholded, not linear,
which contradicted this instrument's definition. Both constants were removed as a
specification defect (FI_METHODS.md Amendment 1). The cap in particular flattened
68.7% of waist-to-height values, making the largest single contributor to the
index effectively binary for two thirds of the cohort.

The bounded form is retained here as a parameterised option so the effect of
bounding can be measured as a sensitivity rather than argued.

This step does not construct signed deficits or aggregate marker deficits into
FI scores.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
EXCURSIONS = ROOT / "outputs/11_fi/fi_adverse_excursions.csv"
SCORING_SPEC = ROOT / "config/fi_scoring_spec.csv"
OUT = ROOT / "outputs/11_fi"

EXPECTED_SHA256 = {
    EXCURSIONS: "7D8FC71445DB1AEA89F426CF38ADDB873FC71F52B33D75A3229C5B21978548B1",
    SCORING_SPEC: "2E096673A2AC23932A566FF93DC48C8755AD77CC86F635DF30DDDC13C10304FD",
}
META = ["SEQN", "RIAGENDR", "RIDAGEYR", "is_fi_reference"]
# Version 2 (FI_METHODS.md Amendment 1): both grading constants removed.
# DEAD_ZONE = 1.0 and CAP = 3.0 were the superseded version-1 values; passing
# them to grade_excursion reproduces the version-1 rule for the sensitivity.
DEAD_ZONE = 0.0
CAP = None  # None means no cap: severity is the excursion itself


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def grade_excursion(
    values: np.ndarray, dead_zone: float = DEAD_ZONE, cap: float | None = CAP
) -> np.ndarray:
    """Grade adverse excursions into deficit severities.

    Version 2 (dead_zone=0, cap=None) returns the excursion unchanged, so
    severity is exactly proportional to displacement at every magnitude.
    Supplying the version-1 constants (1.0, 3.0) reproduces the superseded
    bounded rule for the prespecified sensitivity.
    """
    values = np.asarray(values, dtype=float)
    assert np.isfinite(dead_zone) and dead_zone >= 0
    if cap is None:
        return np.maximum(values - dead_zone, 0.0)
    assert np.isfinite(cap) and dead_zone < cap
    return np.clip((values - dead_zone) / (cap - dead_zone), 0.0, 1.0)


def build_severities(excursions: pd.DataFrame, scoring_spec: pd.DataFrame) -> pd.DataFrame:
    """Apply the same locked grading rule to all marker excursions."""
    markers = scoring_spec.marker.tolist()
    assert list(excursions.columns) == META + markers, "excursion schema or order changed"
    assert len(excursions) == 29053 and excursions.SEQN.is_unique
    values = excursions[markers].to_numpy(float)
    assert np.isfinite(values).all() and (values >= 0).all()

    severities = excursions[META].copy()
    severities[markers] = grade_excursion(values)
    graded = severities[markers].to_numpy(float)
    assert np.isfinite(graded).all()
    assert (graded >= 0).all()  # version 2: no upper bound
    return severities


def build_checks(
    excursions: pd.DataFrame, severities: pd.DataFrame, markers: list[str]
) -> pd.DataFrame:
    """Audit each piece of the grading function over all values."""
    e = excursions[markers].to_numpy(float)
    d = severities[markers].to_numpy(float)
    zero = e <= 0
    positive = e > 0
    expected = grade_excursion(e)
    # Values that WOULD have been flattened by the superseded version-1 rule,
    # reported so the size of the correction is visible in the audit trail.
    former_dead_zone = (e > 0) & (e <= 1.0)
    former_cap = e >= 3.0
    return pd.DataFrame([{
        "participant_marker_values": int(e.size),
        "zero_excursion_values": int(zero.sum()),
        "positive_excursion_values": int(positive.sum()),
        "minimum_severity": float(d.min()),
        "maximum_severity": float(d.max()),
        "maximum_excursion": float(e.max()),
        "maximum_absolute_formula_difference": float(np.max(np.abs(d - expected))),
        "maximum_absolute_zero_excursion_severity": float(np.max(np.abs(d[zero]))),
        # No dead zone: the smallest severity among non-zero excursions must be
        # strictly positive. Under version 1 this was exactly 0 up to one unit.
        "minimum_positive_excursion_severity": float(np.min(d[positive])),
        # No cap: the maximum severity must equal the maximum excursion rather
        # than a ceiling. Under version 1 this ratio was far below 1.
        "max_severity_over_max_excursion": float(d.max() / e.max()),
        "distinct_severity_values": int(np.unique(d).size),
        "values_freed_from_former_dead_zone": int(former_dead_zone.sum()),
        "values_freed_from_former_cap": int(former_cap.sum()),
    }])


def main() -> None:
    for path, expected_hash in EXPECTED_SHA256.items():
        assert fingerprint(path) == expected_hash, f"locked input hash changed: {path.name}"

    excursions = pd.read_csv(EXCURSIONS)
    scoring_spec = pd.read_csv(SCORING_SPEC)
    markers = scoring_spec.marker.tolist()
    severities = build_severities(excursions, scoring_spec)
    checks = build_checks(excursions, severities, markers)
    assert checks.maximum_absolute_formula_difference.iloc[0] == 0
    assert checks.maximum_absolute_zero_excursion_severity.iloc[0] == 0
    # Version-2 gates. No dead zone: every non-zero excursion receives a strictly
    # positive severity. No cap: the maximum severity is the maximum excursion,
    # not a ceiling.
    assert checks.minimum_positive_excursion_severity.iloc[0] > 0, "dead zone present"
    assert checks.max_severity_over_max_excursion.iloc[0] == 1.0, "cap present"

    OUT.mkdir(parents=True, exist_ok=True)
    severities.to_csv(OUT / "fi_absolute_deficit_severities.csv", index=False)
    checks.to_csv(OUT / "fi_deficit_severity_checks.csv", index=False)
    print("FI absolute deficit severities built: 29,053 participants x 38 biomarkers")
    print(checks.to_string(index=False))


if __name__ == "__main__":
    main()
