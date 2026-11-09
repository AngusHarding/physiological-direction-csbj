"""Apply the frozen error screens and complete-case rule to measured CBC data."""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
RAW = ROOT / "outputs" / "01_imported" / "cbc_raw.csv"
OUTPUT = ROOT / "outputs" / "02_harmonised"
MEASURED = ["LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT",
            "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXMCVSI", "LBXRDW",
            "LBXPLTSI", "LBXMPSI"]
BOUNDS = {"LBXPLTSI": (1, 999), "LBXMCVSI": (50, 130), "LBXRDW": (10, 45),
          "LBXLYPCT": (0, 100), "LBXMOPCT": (0, 100), "LBXNEPCT": (0, 100),
          "LBXEOPCT": (0, 100), "LBXBAPCT": (0, 100)}


def main() -> None:
    frame = pd.read_csv(RAW)
    log = []
    # Deliberately leave WBC, RBC, haemoglobin and MPV unbounded: their extreme
    # values can represent genuine disease and should not be silently clipped.
    for variable, (low, high) in BOUNDS.items():
        invalid = frame[variable].notna() & ~frame[variable].between(low, high)
        log.append({"variable": variable, "lower_inclusive": low,
                    "upper_inclusive": high, "values_to_missing": int(invalid.sum())})
        frame.loc[invalid, variable] = np.nan
    complete = frame[MEASURED].notna().all(axis=1)
    clean = frame.loc[complete, ["SEQN", "cycle"] + MEASURED].reset_index(drop=True)
    if clean.SEQN.duplicated().any():
        raise ValueError("harmonised CBC contains duplicate SEQN")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    clean.to_csv(OUTPUT / "whole_blood_count.csv", index=False)
    pd.DataFrame(log).to_csv(OUTPUT / "cbc_error_screen_log.csv", index=False)
    clean.groupby("cycle").size().rename("complete_rows").reset_index().to_csv(OUTPUT / "cbc_by_cycle.csv", index=False)
    print(f"harmonised measured CBC: {len(clean):,} complete rows from {len(frame):,} imported participants")


if __name__ == "__main__":
    main()
