"""Estimate sex-specific healthy-reference parameters for the FI branch.

This step estimates transformed-scale quartiles, medians, and robust scales. It
does not calculate participant coordinates, deficits, or FI scores.
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
SCORING_SPEC = ROOT / "config/fi_scoring_spec.csv"
OUT = ROOT / "outputs/11_fi"

EXPECTED_SHA256 = {
    TRANSFORMED: "E2F404E78414F9C33B9DCB8DA5BC905C3AD25F4C47472E406F13564CB9243EDC",
    MEMBERSHIP: "6D2D62F7C390A1CFB1A88F6C7FC6746F0670D1D7926DECE4AFCFD977C98ADC6A",
    SCORING_SPEC: "2E096673A2AC23932A566FF93DC48C8755AD77CC86F635DF30DDDC13C10304FD",
}
SEX_LABELS = {1.0: "male", 2.0: "female"}
EXPECTED_REFERENCE_BY_SEX = {1.0: 547, 2.0: 691}
NORMAL_IQR = 1.349


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def estimate_parameters(
    transformed: pd.DataFrame, membership: pd.DataFrame, scoring_spec: pd.DataFrame
) -> pd.DataFrame:
    """Estimate Q1, median, Q3, and IQR/1.349 by marker and sex."""
    markers = scoring_spec.marker.tolist()
    assert list(transformed.columns) == ["SEQN", "RIAGENDR", "RIDAGEYR"] + markers
    assert list(membership.columns) == ["SEQN", "RIAGENDR", "RIDAGEYR", "is_fi_reference"]
    assert membership.SEQN.equals(transformed.SEQN), "FI membership order changed"
    assert membership.RIAGENDR.equals(transformed.RIAGENDR), "FI membership sex changed"
    assert membership.RIDAGEYR.equals(transformed.RIDAGEYR), "FI membership age changed"

    rows = []
    for sex_code, sex_label in SEX_LABELS.items():
        reference = membership.is_fi_reference.astype(bool) & membership.RIAGENDR.eq(sex_code)
        assert int(reference.sum()) == EXPECTED_REFERENCE_BY_SEX[sex_code]
        for spec_row in scoring_spec.itertuples(index=False):
            values = transformed.loc[reference, spec_row.marker].to_numpy(float)
            assert len(values) == EXPECTED_REFERENCE_BY_SEX[sex_code]
            assert np.isfinite(values).all(), f"non-finite reference value: {spec_row.marker}/{sex_label}"
            q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75], method="linear")
            iqr = q3 - q1
            robust_scale = iqr / NORMAL_IQR
            assert np.isfinite(robust_scale) and robust_scale > 0, (
                f"invalid robust scale: {spec_row.marker}/{sex_label}"
            )
            rows.append({
                "marker": spec_row.marker,
                "sex": sex_label,
                "sex_code": int(sex_code),
                "selected_transform": spec_row.selected_transform,
                "orientation": int(spec_row.orientation),
                "adverse_direction": spec_row.adverse_direction,
                "physiological_system": spec_row.physiological_system,
                "reference_q1_transformed": float(q1),
                "reference_median_transformed": float(median),
                "reference_q3_transformed": float(q3),
                "reference_iqr_transformed": float(iqr),
                "reference_robust_scale": float(robust_scale),
                "reference_n": len(values),
                "quantile_method": "linear",
                "normal_iqr_divisor": NORMAL_IQR,
            })

    parameters = pd.DataFrame(rows)
    assert parameters.shape == (76, 15)
    assert not parameters.duplicated(["marker", "sex_code"]).any()
    assert parameters.groupby("sex_code", sort=False).marker.apply(list).tolist() == [markers, markers]
    return parameters


def build_checks(parameters: pd.DataFrame) -> pd.DataFrame:
    """Summarise parameter completeness and scale validity by sex."""
    return pd.DataFrame([
        {
            "sex": sex_label,
            "reference_n": int(group.reference_n.iloc[0]),
            "biomarkers": len(group),
            "nonpositive_robust_scales": int(group.reference_robust_scale.le(0).sum()),
            "nonfinite_parameter_values": int((~np.isfinite(group[[
                "reference_q1_transformed", "reference_median_transformed",
                "reference_q3_transformed", "reference_iqr_transformed",
                "reference_robust_scale",
            ]].to_numpy(float))).sum()),
            "minimum_robust_scale": float(group.reference_robust_scale.min()),
            "maximum_robust_scale": float(group.reference_robust_scale.max()),
        }
        for sex_label, group in parameters.groupby("sex", sort=False)
    ])


def main() -> None:
    for path, expected_hash in EXPECTED_SHA256.items():
        assert fingerprint(path) == expected_hash, f"locked input hash changed: {path.name}"

    parameters = estimate_parameters(
        pd.read_csv(TRANSFORMED), pd.read_csv(MEMBERSHIP), pd.read_csv(SCORING_SPEC)
    )
    checks = build_checks(parameters)
    OUT.mkdir(parents=True, exist_ok=True)
    parameters.to_csv(OUT / "fi_reference_robust_parameters.csv", index=False)
    checks.to_csv(OUT / "fi_reference_parameter_checks.csv", index=False)
    print("FI reference parameters estimated: 38 biomarkers x 2 sexes")
    print(checks.to_string(index=False))


if __name__ == "__main__":
    main()
