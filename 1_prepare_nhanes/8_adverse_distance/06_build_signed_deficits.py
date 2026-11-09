"""Construct the signed graded-deficit matrix for the FI branch.

This step restores adverse-tail sign to the deficit severities. It does not
aggregate markers into participant-level or subsystem FI scores.

Version 2 (FI_METHODS.md Amendment 1): severities are no longer bounded to
[0, 1], so signed deficits are no longer bounded to [-1, 1]. The invariant that
matters is unchanged and still enforced: |signed| == severity exactly.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
COORDINATES = ROOT / "outputs/11_fi/fi_raw_oriented_robust_coordinates.csv"
SEVERITIES = ROOT / "outputs/11_fi/fi_absolute_deficit_severities.csv"
SCORING_SPEC = ROOT / "config/fi_scoring_spec.csv"
OUT = ROOT / "outputs/11_fi"

EXPECTED_SHA256 = {
    COORDINATES: "4F5E71E2C50A4307858AEA66C0AE5CBDAE064CAFB5AB66BED090E5194B7A028C",
    SEVERITIES: "F3BA6949DE3E517F5733D85ED3DBC4194CA0EA93622E4A1821710B388834874A",
    SCORING_SPEC: "2E096673A2AC23932A566FF93DC48C8755AD77CC86F635DF30DDDC13C10304FD",
}
META = ["SEQN", "RIAGENDR", "RIDAGEYR", "is_fi_reference"]


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def signed_deficit(
    coordinates: np.ndarray, severities: np.ndarray, direction: str
) -> np.ndarray:
    """Restore the direction-specific sign while preserving severity magnitude."""
    coordinates = np.asarray(coordinates, dtype=float)
    severities = np.asarray(severities, dtype=float)
    if direction == "high_only":
        return severities
    if direction == "low_only":
        return -severities
    if direction == "two_sided":
        return np.sign(coordinates) * severities
    raise ValueError(f"invalid adverse direction: {direction}")


def build_signed_deficits(
    coordinates: pd.DataFrame, severities: pd.DataFrame, scoring_spec: pd.DataFrame
) -> pd.DataFrame:
    """Build signed deficits for every participant-marker pair."""
    markers = scoring_spec.marker.tolist()
    expected_columns = META + markers
    assert list(coordinates.columns) == expected_columns, "coordinate schema or order changed"
    assert list(severities.columns) == expected_columns, "severity schema or order changed"
    assert coordinates[META].equals(severities[META]), "coordinate and severity identities differ"
    assert len(coordinates) == 29053 and coordinates.SEQN.is_unique

    signed = coordinates[META].copy()
    for row in scoring_spec.itertuples(index=False):
        signed[row.marker] = signed_deficit(
            coordinates[row.marker].to_numpy(float),
            severities[row.marker].to_numpy(float),
            row.adverse_direction,
        )

    values = signed[markers].to_numpy(float)
    absolute = severities[markers].to_numpy(float)
    assert np.isfinite(values).all(), "signed deficits are non-finite"
    # Version 2: no [-1, 1] bound; severities are unbounded above by design.
    assert np.array_equal(np.abs(values), absolute), "signed magnitude differs from severity"
    return signed


def build_checks(
    coordinates: pd.DataFrame,
    severities: pd.DataFrame,
    signed: pd.DataFrame,
    scoring_spec: pd.DataFrame,
) -> pd.DataFrame:
    """Audit sign identities separately for each adverse-direction class."""
    rows = []
    for direction in ("high_only", "low_only", "two_sided"):
        markers = scoring_spec.loc[
            scoring_spec.adverse_direction.eq(direction), "marker"
        ].tolist()
        z = coordinates[markers].to_numpy(float)
        d = severities[markers].to_numpy(float)
        g = signed[markers].to_numpy(float)
        expected = signed_deficit(z, d, direction)
        nonzero = d > 0
        rows.append({
            "adverse_direction": direction,
            "biomarkers": len(markers),
            "participant_marker_values": int(g.size),
            "nonzero_signed_deficits": int(nonzero.sum()),
            "negative_signed_deficits": int((g < 0).sum()),
            "positive_signed_deficits": int((g > 0).sum()),
            "maximum_absolute_formula_difference": float(np.max(np.abs(g - expected))),
            "maximum_absolute_magnitude_difference": float(np.max(np.abs(np.abs(g) - d))),
        })
    return pd.DataFrame(rows)


def main() -> None:
    for path, expected_hash in EXPECTED_SHA256.items():
        assert fingerprint(path) == expected_hash, f"locked input hash changed: {path.name}"

    coordinates = pd.read_csv(COORDINATES)
    severities = pd.read_csv(SEVERITIES)
    scoring_spec = pd.read_csv(SCORING_SPEC)
    signed = build_signed_deficits(coordinates, severities, scoring_spec)
    checks = build_checks(coordinates, severities, signed, scoring_spec)
    assert checks.maximum_absolute_formula_difference.max() == 0
    assert checks.maximum_absolute_magnitude_difference.max() == 0

    OUT.mkdir(parents=True, exist_ok=True)
    signed.to_csv(OUT / "fi_signed_graded_deficits.csv", index=False)
    checks.to_csv(OUT / "fi_signed_deficit_checks.csv", index=False)
    print("FI signed graded deficits built: 29,053 participants x 38 biomarkers")
    print(checks.to_string(index=False))


if __name__ == "__main__":
    main()
