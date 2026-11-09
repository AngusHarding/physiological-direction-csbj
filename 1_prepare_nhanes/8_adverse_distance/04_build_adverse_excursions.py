"""Map raw-oriented robust coordinates to non-negative adverse excursions.

This step applies marker directionality only. It does not apply the dead zone,
cap excursions, construct signed deficits, or calculate FI scores.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
COORDINATES = ROOT / "outputs/11_fi/fi_raw_oriented_robust_coordinates.csv"
SCORING_SPEC = ROOT / "config/fi_scoring_spec.csv"
OUT = ROOT / "outputs/11_fi"

EXPECTED_SHA256 = {
    COORDINATES: "4F5E71E2C50A4307858AEA66C0AE5CBDAE064CAFB5AB66BED090E5194B7A028C",
    SCORING_SPEC: "2E096673A2AC23932A566FF93DC48C8755AD77CC86F635DF30DDDC13C10304FD",
}
META = ["SEQN", "RIAGENDR", "RIDAGEYR", "is_fi_reference"]
VALID_DIRECTIONS = {"high_only", "low_only", "two_sided"}


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def adverse_excursion(values: np.ndarray, direction: str) -> np.ndarray:
    """Return adverse distance from zero under one locked direction rule."""
    values = np.asarray(values, dtype=float)
    if direction == "high_only":
        return np.maximum(0.0, values)
    if direction == "low_only":
        return np.maximum(0.0, -values)
    if direction == "two_sided":
        return np.abs(values)
    raise ValueError(f"invalid adverse direction: {direction}")


def build_excursions(coordinates: pd.DataFrame, scoring_spec: pd.DataFrame) -> pd.DataFrame:
    """Apply each marker's direction rule without grading or capping."""
    markers = scoring_spec.marker.tolist()
    assert list(coordinates.columns) == META + markers, "coordinate schema or order changed"
    assert len(coordinates) == 29053 and coordinates.SEQN.is_unique
    assert scoring_spec.marker.is_unique and len(scoring_spec) == 38
    assert set(scoring_spec.adverse_direction) <= VALID_DIRECTIONS
    assert np.isfinite(coordinates[markers].to_numpy(float)).all()

    excursions = coordinates[META].copy()
    for row in scoring_spec.itertuples(index=False):
        excursions[row.marker] = adverse_excursion(
            coordinates[row.marker].to_numpy(float), row.adverse_direction
        )

    values = excursions[markers].to_numpy(float)
    assert np.isfinite(values).all(), "adverse excursions are non-finite"
    assert (values >= 0).all(), "adverse excursions must be non-negative"
    return excursions


def build_checks(
    coordinates: pd.DataFrame, excursions: pd.DataFrame, scoring_spec: pd.DataFrame
) -> pd.DataFrame:
    """Summarise direction-specific identities and protective-side zeros."""
    rows = []
    for direction in ("high_only", "low_only", "two_sided"):
        markers = scoring_spec.loc[
            scoring_spec.adverse_direction.eq(direction), "marker"
        ].tolist()
        z = coordinates[markers].to_numpy(float)
        e = excursions[markers].to_numpy(float)
        if direction == "high_only":
            expected = np.maximum(0.0, z)
            protective_zero = e[z < 0]
        elif direction == "low_only":
            expected = np.maximum(0.0, -z)
            protective_zero = e[z > 0]
        else:
            expected = np.abs(z)
            protective_zero = np.array([], dtype=float)
        rows.append({
            "adverse_direction": direction,
            "biomarkers": len(markers),
            "coordinate_values": int(z.size),
            "maximum_absolute_formula_difference": float(np.max(np.abs(e - expected))),
            "protective_side_values": int(protective_zero.size),
            "maximum_absolute_protective_side_excursion": (
                float(np.max(np.abs(protective_zero))) if protective_zero.size else 0.0
            ),
            "maximum_excursion": float(e.max()),
        })
    return pd.DataFrame(rows)


def main() -> None:
    for path, expected_hash in EXPECTED_SHA256.items():
        assert fingerprint(path) == expected_hash, f"locked input hash changed: {path.name}"

    coordinates = pd.read_csv(COORDINATES)
    scoring_spec = pd.read_csv(SCORING_SPEC)
    excursions = build_excursions(coordinates, scoring_spec)
    checks = build_checks(coordinates, excursions, scoring_spec)
    assert checks.maximum_absolute_formula_difference.max() == 0
    assert checks.maximum_absolute_protective_side_excursion.max() == 0

    OUT.mkdir(parents=True, exist_ok=True)
    excursions.to_csv(OUT / "fi_adverse_excursions.csv", index=False)
    checks.to_csv(OUT / "fi_adverse_excursion_checks.csv", index=False)
    print("FI adverse excursions built: 29,053 participants x 38 biomarkers")
    print(checks.to_string(index=False))


if __name__ == "__main__":
    main()
