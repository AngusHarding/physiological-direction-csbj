# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/3_body_and_bp/05_harmonise_blood_pressure.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Discards reading 1, sets zero readings to missing, averages
#               readings 2-4 (at least one usable), and applies the pressure
#               and pulse bounds (pulse as released, BPXPLS, in beats/min).
#
# METHODS       Supp. Section Zero, "Blood pressure and pulse" and "Cleaning".
# =============================================================================
"""Reconstruct the frozen cardiovascular biomarkers with explicit QC logs."""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
RAW = ROOT / "outputs" / "01_imported" / "cardiovascular_raw.csv"
OUTPUT = ROOT / "outputs" / "02_harmonised"
SYS = ["BPXSY2", "BPXSY3", "BPXSY4"]
DIA = ["BPXDI2", "BPXDI3", "BPXDI4"]
BOUNDS = {"systolic": (60, 290), "diastolic": (30, 150), "BPXPLS": (30, 220)}


def main() -> None:
    frame = pd.read_csv(RAW)
    changes = []

    # A literal zero is meaningful in some laboratory assays, but it cannot be
    # a valid resting blood pressure. Replace it before calculating a mean.
    for variable in SYS + DIA:
        zero = frame[variable].eq(0)
        changes.append({"rule_id": "CARDIO-02", "variable": variable,
                        "values_changed": int(zero.sum()), "detail": "zero -> missing"})
        frame.loc[zero, variable] = np.nan

    # count() records how many of readings 2, 3 and 4 contribute to each mean.
    frame["systolic_readings_used"] = frame[SYS].count(axis=1)
    frame["diastolic_readings_used"] = frame[DIA].count(axis=1)
    frame["systolic"] = frame[SYS].mean(axis=1)
    frame["diastolic"] = frame[DIA].mean(axis=1)

    # map(), unlike an if/else lambda, leaves missing values missing.
    before_codes = frame["BPXPULS"].value_counts(dropna=False)
    frame["BPXPULS"] = frame["BPXPULS"].map({1.0: 0.0, 2.0: 1.0})
    for label, count in before_codes.items():
        changes.append({"rule_id": "CARDIO-04", "variable": "BPXPULS",
                        "values_changed": int(count), "detail": f"raw code {label}"})

    # These are deliberately broad error bounds, not healthy-range exclusions.
    for variable, (low, high) in BOUNDS.items():
        invalid = frame[variable].notna() & ~frame[variable].between(low, high)
        changes.append({"rule_id": f"CARDIO-bound-{variable}", "variable": variable,
                        "values_changed": int(invalid.sum()),
                        "detail": f"outside inclusive [{low}, {high}] -> missing"})
        frame.loc[invalid, variable] = np.nan

    biomarkers = ["systolic", "diastolic", "BPXPLS"]
    complete = frame[biomarkers].notna().all(axis=1)
    clean = frame.loc[complete, ["SEQN", "cycle"] + biomarkers + ["BPXPULS"]].reset_index(drop=True)
    if clean.SEQN.duplicated().any():
        raise ValueError("harmonised cardiovascular output contains duplicate SEQN")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    clean.to_csv(OUTPUT / "cardiovascular.csv", index=False)
    pd.DataFrame(changes).to_csv(OUTPUT / "cardiovascular_change_log.csv", index=False)
    frame.groupby(["cycle", "systolic_readings_used", "diastolic_readings_used"], dropna=False).size().rename("participants").reset_index().to_csv(OUTPUT / "cardiovascular_reading_counts.csv", index=False)
    clean.groupby("cycle").size().rename("complete_rows").reset_index().to_csv(OUTPUT / "cardiovascular_by_cycle.csv", index=False)
    print(f"harmonised cardiovascular: {len(clean):,} complete rows from {len(frame):,} imported participants")
    print(f"BPXPULS missing but deliberately retained: {clean.BPXPULS.isna().sum():,}")


if __name__ == "__main__":
    main()
