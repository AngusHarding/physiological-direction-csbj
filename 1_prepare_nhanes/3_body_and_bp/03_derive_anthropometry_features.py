"""Derive BMI for reference screening and WHtR for biomarker scoring."""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)


def main() -> None:
    raw = pd.read_csv(ROOT / "outputs/02_harmonised/anthropometry.csv")
    derived = raw[["SEQN", "cycle"]].copy()
    derived["BMI"] = raw.BMXWT / (raw.BMXHT / 100) ** 2
    derived["WHtR"] = raw.BMXWAIST / raw.BMXHT
    derived = derived.replace([np.inf, -np.inf], np.nan)
    output = ROOT / "outputs/03_derived"
    output.mkdir(parents=True, exist_ok=True)
    derived.to_csv(output / "anthropometry_derived_features.csv", index=False)
    pd.DataFrame([{"feature": "BMI", "complete_values": int(derived.BMI.notna().sum()), "role": "reference_screen"},
                  {"feature": "WHtR", "complete_values": int(derived.WHtR.notna().sum()), "role": "scored_biomarker"}]).to_csv(output / "anthropometry_derivation_summary.csv", index=False)
    print(f"anthropometry derivation: BMI={derived.BMI.notna().sum():,}; WHtR={derived.WHtR.notna().sum():,}")


if __name__ == "__main__":
    main()
