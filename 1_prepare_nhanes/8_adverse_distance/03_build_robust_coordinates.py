# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/8_adverse_distance/03_build_robust_coordinates.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  For every participant and marker, r = orientation x
#               (transformed value - same-sex reference median) / (IQR/1.349).
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
"""Build raw-oriented robust coordinates for the FI branch.

This step does not apply adverse-tail rules, dead zones, caps, or FI scoring.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
TRANSFORMED = ROOT / "outputs/06_transform/final_38_panel_transformed_unstandardised_v2.csv"
MEMBERSHIP = ROOT / "outputs/11_fi/fi_reference_membership.csv"
PARAMETERS = ROOT / "outputs/11_fi/fi_reference_robust_parameters.csv"
SCORING_SPEC = ROOT / "config/fi_scoring_spec.csv"
OUT = ROOT / "outputs/11_fi"

EXPECTED_SHA256 = {
    TRANSFORMED: "E2F404E78414F9C33B9DCB8DA5BC905C3AD25F4C47472E406F13564CB9243EDC",
    MEMBERSHIP: "6D2D62F7C390A1CFB1A88F6C7FC6746F0670D1D7926DECE4AFCFD977C98ADC6A",
    PARAMETERS: "77852C65DD36434E2806BE8F387BE0E71AD51A5512282C7D6FA0E62688839219",
    SCORING_SPEC: "2E096673A2AC23932A566FF93DC48C8755AD77CC86F635DF30DDDC13C10304FD",
}
SEX_LABELS = {1: "male", 2: "female"}
EXPECTED_ROWS = 29053
EXPECTED_REFERENCE_BY_SEX = {1: 547, 2: 691}


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


# -- CALCULATION: the robust coordinate r ------------------------------------
# Supp. Section Zero: "The orientation oj is +1 when Tj increases with x and -1
#   when it decreases, so r always increases with the raw value."
def build_coordinates(
    transformed: pd.DataFrame,
    membership: pd.DataFrame,
    parameters: pd.DataFrame,
    scoring_spec: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate orientation * (transformed value - median) / robust scale."""
    markers = scoring_spec.marker.tolist()
    assert len(transformed) == EXPECTED_ROWS
    assert list(transformed.columns) == ["SEQN", "RIAGENDR", "RIDAGEYR"] + markers
    assert membership.SEQN.equals(transformed.SEQN), "FI membership order changed"
    assert membership.RIAGENDR.equals(transformed.RIAGENDR), "FI membership sex changed"
    assert membership.RIDAGEYR.equals(transformed.RIDAGEYR), "FI membership age changed"
    assert parameters.shape[0] == len(markers) * len(SEX_LABELS)
    assert not parameters.duplicated(["marker", "sex_code"]).any()

    coordinates = membership.copy()
    for marker in markers:
        coordinates[marker] = np.nan
        spec_row = scoring_spec.loc[scoring_spec.marker.eq(marker)].iloc[0]
        for sex_code in SEX_LABELS:
            sex = transformed.RIAGENDR.eq(sex_code)
            parameter = parameters.loc[
                parameters.marker.eq(marker) & parameters.sex_code.eq(sex_code)
            ]
            assert len(parameter) == 1, f"missing or duplicate parameter: {marker}/{sex_code}"
            parameter = parameter.iloc[0]
            assert parameter.selected_transform == spec_row.selected_transform
            assert int(parameter.orientation) == int(spec_row.orientation)
            median = float(parameter.reference_median_transformed)
            scale = float(parameter.reference_robust_scale)
            assert np.isfinite(median) and np.isfinite(scale) and scale > 0
            coordinates.loc[sex, marker] = (
                int(spec_row.orientation)
                * (transformed.loc[sex, marker] - median)
                / scale
            )

    assert coordinates.shape == (EXPECTED_ROWS, 4 + len(markers))
    assert np.isfinite(coordinates[markers].to_numpy(float)).all(), "coordinates are non-finite"
    return coordinates


def build_checks(coordinates: pd.DataFrame, markers: list[str]) -> pd.DataFrame:
    """Check median-zero and IQR-1.349 behaviour in each reference sex."""
    rows = []
    for sex_code, sex_label in SEX_LABELS.items():
        reference = coordinates.is_fi_reference.astype(bool) & coordinates.RIAGENDR.eq(sex_code)
        values = coordinates.loc[reference, markers]
        medians = values.median()
        iqrs = values.quantile(0.75, interpolation="linear") - values.quantile(
            0.25, interpolation="linear"
        )
        rows.append({
            "sex": sex_label,
            "reference_n": int(reference.sum()),
            "biomarkers": len(markers),
            "maximum_absolute_reference_median": float(medians.abs().max()),
            "maximum_absolute_reference_iqr_minus_1_349": float((iqrs - 1.349).abs().max()),
            "nonfinite_coordinates": int((~np.isfinite(values.to_numpy(float))).sum()),
        })
        assert int(reference.sum()) == EXPECTED_REFERENCE_BY_SEX[sex_code]
    return pd.DataFrame(rows)


def main() -> None:
    for path, expected_hash in EXPECTED_SHA256.items():
        assert fingerprint(path) == expected_hash, f"locked input hash changed: {path.name}"

    scoring_spec = pd.read_csv(SCORING_SPEC)
    markers = scoring_spec.marker.tolist()
    coordinates = build_coordinates(
        pd.read_csv(TRANSFORMED),
        pd.read_csv(MEMBERSHIP),
        pd.read_csv(PARAMETERS),
        scoring_spec,
    )
    checks = build_checks(coordinates, markers)
    assert checks.maximum_absolute_reference_median.max() < 1e-12
    assert checks.maximum_absolute_reference_iqr_minus_1_349.max() < 1e-12

    OUT.mkdir(parents=True, exist_ok=True)
    coordinates.to_csv(OUT / "fi_raw_oriented_robust_coordinates.csv", index=False)
    checks.to_csv(OUT / "fi_robust_coordinate_checks.csv", index=False)
    print("FI robust coordinates built: 29,053 participants x 38 biomarkers")
    print(checks.to_string(index=False))


if __name__ == "__main__":
    main()
