"""Standardise locked transformed biomarkers to the healthy same-sex reference.

This stage creates coordinates for the HD branch only. It does not calculate
Mahalanobis distance, covariance, FI deficits, trimming, or winsorisation.
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
RECIPE = ROOT / "config/final_38_transformations_v2.csv"
OUT = ROOT / "outputs/09_standardised"

EXPECTED_TRANSFORMED_SHA256 = "E2F404E78414F9C33B9DCB8DA5BC905C3AD25F4C47472E406F13564CB9243EDC"
EXPECTED_REFERENCE_SHA256 = "9CFDA42A07F3DB275547A7014EFF7C11EBA8BD2F59DF93332B30CB218FEC15B2"
SEX_LABELS = {1.0: "male", 2.0: "female"}


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    if fingerprint(TRANSFORMED) != EXPECTED_TRANSFORMED_SHA256:
        raise AssertionError("locked transformed-panel hash changed")
    if fingerprint(REFERENCE) != EXPECTED_REFERENCE_SHA256:
        raise AssertionError("locked reference-membership hash changed")

    transformed = pd.read_csv(TRANSFORMED)
    reference = pd.read_csv(REFERENCE)
    recipe = pd.read_csv(RECIPE)
    markers = recipe.marker.tolist()
    expected = ["SEQN", "RIAGENDR", "RIDAGEYR"] + markers
    if list(transformed.columns) != expected or transformed.shape != (29053, 41):
        raise AssertionError("transformed panel structure changed")
    if reference.SEQN.duplicated().any() or len(reference) != len(transformed):
        raise AssertionError("reference membership structure changed")

    joined = transformed.merge(reference[["SEQN", "is_reference"]], on="SEQN",
                               how="left", validate="one_to_one", sort=False)
    if not joined.SEQN.equals(transformed.SEQN) or joined.is_reference.isna().any():
        raise AssertionError("reference join changed participant identity")

    standardised = joined[["SEQN", "RIAGENDR", "RIDAGEYR", "is_reference"]].copy()
    parameter_rows = []
    for marker in markers:
        standardised[marker] = np.nan
        for sex_code, sex_label in SEX_LABELS.items():
            sex = joined.RIAGENDR.eq(sex_code)
            ref = sex & joined.is_reference.astype(bool)
            values = joined.loc[ref, marker]
            mean = float(values.mean())
            sd = float(values.std(ddof=1))
            if not np.isfinite(mean) or not np.isfinite(sd) or sd <= 0:
                raise AssertionError(f"invalid reference parameters: {marker}/{sex_label}")
            standardised.loc[sex, marker] = (joined.loc[sex, marker] - mean) / sd
            recipe_row = recipe.loc[recipe.marker.eq(marker)].iloc[0]
            parameter_rows.append({"marker": marker, "sex": sex_label,
                                   "sex_code": int(sex_code),
                                   "transform": recipe_row.selected_transform,
                                   "orientation": int(recipe_row.orientation),
                                   "reference_mean_transformed": mean,
                                   "reference_sd_transformed": sd,
                                   "reference_n": int(ref.sum())})

    if standardised[markers].isna().any().any() or not np.isfinite(
            standardised[markers].to_numpy()).all():
        raise AssertionError("standardisation created missing or infinite values")

    parameters = pd.DataFrame(parameter_rows)
    checks = []
    for sex_code, sex_label in SEX_LABELS.items():
        ref = standardised.RIAGENDR.eq(sex_code) & standardised.is_reference.astype(bool)
        checks.append({"sex": sex_label, "reference_n": int(ref.sum()),
                       "maximum_absolute_marker_mean": float(standardised.loc[ref, markers].mean().abs().max()),
                       "maximum_absolute_sd_minus_one": float((standardised.loc[ref, markers].std(ddof=1) - 1).abs().max())})

    OUT.mkdir(parents=True, exist_ok=True)
    standardised.to_csv(OUT / "final_38_panel_transformed_standardised_hd.csv", index=False)
    parameters.to_csv(OUT / "hd_reference_standardisation_parameters.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT / "hd_reference_standardisation_checks.csv", index=False)
    print(f"HD-standardised panel: {len(standardised):,} participants x {len(markers)} biomarkers")
    print(parameters.groupby("sex").reference_n.first().to_string())
    print(pd.DataFrame(checks).to_string(index=False))


if __name__ == "__main__":
    main()
