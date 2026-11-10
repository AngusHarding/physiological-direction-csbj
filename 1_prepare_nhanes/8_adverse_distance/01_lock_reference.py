# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/8_adverse_distance/01_lock_reference.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Confirms the young-healthy reference members used for adverse
#               distance are the same 1,238 (547 men, 691 women) as everywhere
#               else, all under 30 with complete markers. Estimates nothing.
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
# =============================================================================
"""Lock and validate healthy-reference membership for the FI branch.

This step does not estimate medians, robust scales, deficits, or FI scores.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
TRANSFORMED = ROOT / "outputs/06_transform/final_38_panel_transformed_unstandardised_v2.csv"
REFERENCE = ROOT / "outputs/08_reference/healthy_reference_members.csv"
SCORING_SPEC = ROOT / "config/fi_scoring_spec.csv"
OUT = ROOT / "outputs/11_fi"

EXPECTED_SHA256 = {
    TRANSFORMED: "E2F404E78414F9C33B9DCB8DA5BC905C3AD25F4C47472E406F13564CB9243EDC",
    REFERENCE: "9CFDA42A07F3DB275547A7014EFF7C11EBA8BD2F59DF93332B30CB218FEC15B2",
    SCORING_SPEC: "2E096673A2AC23932A566FF93DC48C8755AD77CC86F635DF30DDDC13C10304FD",
}
EXPECTED_ROWS = 29053
EXPECTED_REFERENCE_BY_SEX = {1.0: 547, 2.0: 691}


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def build_membership(
    transformed: pd.DataFrame, reference: pd.DataFrame, markers: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the FI membership artifact and its sex-specific audit summary."""
    assert len(transformed) == EXPECTED_ROWS, "transformed cohort row count changed"
    assert transformed.SEQN.is_unique, "transformed cohort contains duplicate SEQN"
    assert reference.SEQN.is_unique, "reference artifact contains duplicate SEQN"
    assert len(reference) == len(transformed), "reference and transformed row counts differ"
    assert list(transformed.columns) == ["SEQN", "RIAGENDR", "RIDAGEYR"] + markers, (
        "transformed panel schema or marker order changed"
    )
    assert set(reference.columns) == {"SEQN", "RIAGENDR", "RIDAGEYR", "is_reference"}, (
        "unexpected reference membership schema"
    )

    joined = transformed[["SEQN", "RIAGENDR", "RIDAGEYR"]].merge(
        reference, on="SEQN", how="left", validate="one_to_one", sort=False,
        suffixes=("", "_reference")
    )
    assert joined.SEQN.equals(transformed.SEQN), "reference join changed participant order"
    assert joined.is_reference.notna().all(), "reference membership is missing for a participant"
    assert joined.RIAGENDR.equals(joined.RIAGENDR_reference), "sex differs between artifacts"
    assert joined.RIDAGEYR.equals(joined.RIDAGEYR_reference), "age differs between artifacts"
    assert set(joined.RIAGENDR.unique()) == set(EXPECTED_REFERENCE_BY_SEX), "unexpected sex code"

    is_reference = joined.is_reference.astype(bool)
    assert transformed[markers].notna().all().all(), "FI transformed panel contains missing values"
    assert np.isfinite(transformed[markers].to_numpy(float)).all(), "FI transformed panel is non-finite"
    assert transformed.loc[is_reference, markers].notna().all().all(), "reference marker is missing"
    assert transformed.loc[is_reference, "RIDAGEYR"].lt(30).all(), "reference contains age 30 or older"

    counts = joined.loc[is_reference].groupby("RIAGENDR").size().to_dict()
    assert counts == EXPECTED_REFERENCE_BY_SEX, "sex-specific reference counts changed"

    membership = transformed[["SEQN", "RIAGENDR", "RIDAGEYR"]].copy()
    membership["is_fi_reference"] = is_reference.to_numpy()
    summary = pd.DataFrame([
        {
            "sex": "male" if sex_code == 1.0 else "female",
            "sex_code": int(sex_code),
            "cohort_n": int(transformed.RIAGENDR.eq(sex_code).sum()),
            "reference_n": int(counts[sex_code]),
            "biomarkers_complete": len(markers),
            "minimum_reference_age": float(transformed.loc[is_reference & transformed.RIAGENDR.eq(sex_code), "RIDAGEYR"].min()),
            "maximum_reference_age": float(transformed.loc[is_reference & transformed.RIAGENDR.eq(sex_code), "RIDAGEYR"].max()),
        }
        for sex_code in EXPECTED_REFERENCE_BY_SEX
    ])
    return membership, summary


def main() -> None:
    for path, expected_hash in EXPECTED_SHA256.items():
        assert fingerprint(path) == expected_hash, f"locked input hash changed: {path.name}"

    transformed = pd.read_csv(TRANSFORMED)
    reference = pd.read_csv(REFERENCE)
    markers = pd.read_csv(SCORING_SPEC).marker.tolist()
    membership, summary = build_membership(transformed, reference, markers)

    OUT.mkdir(parents=True, exist_ok=True)
    membership.to_csv(OUT / "fi_reference_membership.csv", index=False)
    summary.to_csv(OUT / "fi_reference_validation.csv", index=False)
    print("FI reference locked: 1,238 participants (547 male; 691 female)")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
