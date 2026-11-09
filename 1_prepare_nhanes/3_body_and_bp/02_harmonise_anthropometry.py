"""Retain any participant with at least one measured body dimension."""

from pathlib import Path
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
MEASURED = ["BMXWT", "BMXHT", "BMXWAIST"]


def main() -> None:
    raw = pd.read_csv(ROOT / "outputs/01_imported/anthropometry_raw.csv")
    # Do not require all three here: BMI needs weight+height, whereas WHtR needs
    # waist+height. Keeping partial rows preserves every potentially usable case.
    clean = raw.dropna(subset=MEASURED, how="all").reset_index(drop=True)
    output = ROOT / "outputs/02_harmonised"
    output.mkdir(parents=True, exist_ok=True)
    clean.to_csv(output / "anthropometry.csv", index=False)
    clean.groupby("cycle").size().rename("retained_rows").reset_index().to_csv(output / "anthropometry_by_cycle.csv", index=False)
    pd.DataFrame([{"imported_rows": len(raw), "fully_empty_removed": len(raw)-len(clean),
                   "retained_rows": len(clean)}]).to_csv(output / "anthropometry_retention_log.csv", index=False)
    print(f"harmonised anthropometry: {len(clean):,} rows; removed {len(raw)-len(clean):,} fully empty exams")


if __name__ == "__main__":
    main()
