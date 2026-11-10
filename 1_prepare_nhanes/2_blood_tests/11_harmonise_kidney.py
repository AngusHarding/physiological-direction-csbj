# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/2_blood_tests/11_harmonise_kidney.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Keeps complete urine pairs and codes the kidney-failure
#               answer as yes / no / unknown.
#
# METHODS       Supp. Section Zero, "Renal, body-size and lipid measures".
# =============================================================================
"""Complete-case urine inputs and materialise valid binary kidney flags."""

from pathlib import Path
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
IMPORTED = ROOT / "outputs" / "01_imported"
OUTPUT = ROOT / "outputs" / "02_harmonised"


def main() -> None:
    urine = pd.read_csv(IMPORTED / "kidney_urine_raw.csv")
    flags = pd.read_csv(IMPORTED / "kidney_failure_flag_raw.csv")
    urine_clean = urine.dropna(subset=["URXUMA", "URXUCR"]).reset_index(drop=True)
    flags_clean = flags.dropna(subset=["kidney_fail"]).reset_index(drop=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    urine_clean.to_csv(OUTPUT / "urine_albumin_creatinine.csv", index=False)
    flags_clean.to_csv(OUTPUT / "kidney_failure_flag.csv", index=False)
    urine_clean.groupby("cycle").size().rename("complete_rows").reset_index().to_csv(OUTPUT / "kidney_urine_by_cycle.csv", index=False)
    flags_clean.groupby(["cycle", "kidney_fail"]).size().rename("participants").reset_index().to_csv(OUTPUT / "kidney_flag_counts.csv", index=False)
    print(f"kidney harmonisation: {len(urine_clean):,} complete urine pairs; {len(flags_clean):,} valid flags")


if __name__ == "__main__":
    main()
