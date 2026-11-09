"""Canonicalise HDL names and complete-case the two lipid source strata."""

from pathlib import Path
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)


def main() -> None:
    imported = ROOT / "outputs/01_imported"
    full = pd.read_csv(imported / "cholesterol_hdl_raw.csv")
    fasting = pd.read_csv(imported / "fasting_triglycerides_raw.csv")
    full["LBDHDD"] = full.pop("HDL_source_value")
    clean = full.dropna(subset=["LBXTC", "LBDHDD"])[["SEQN", "cycle", "LBXTC", "LBDHDD"]].reset_index(drop=True)
    fasting_clean = fasting.dropna(subset=["LBXTR"])[["SEQN", "cycle", "LBXTR"]].reset_index(drop=True)
    output = ROOT / "outputs/02_harmonised"
    output.mkdir(parents=True, exist_ok=True)
    clean.to_csv(output / "cholesterol_hdl.csv", index=False)
    fasting_clean.to_csv(output / "fasting_triglycerides.csv", index=False)
    clean.groupby("cycle").size().rename("complete_rows").reset_index().to_csv(output / "cholesterol_hdl_by_cycle.csv", index=False)
    fasting_clean.groupby("cycle").size().rename("complete_rows").reset_index().to_csv(output / "fasting_triglycerides_by_cycle.csv", index=False)
    print(f"lipid harmonisation: {len(clean):,} complete TC/HDL; {len(fasting_clean):,} fasting triglycerides")


if __name__ == "__main__":
    main()
