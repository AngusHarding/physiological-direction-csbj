# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/2_blood_tests/12_derive_kidney_features.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Derives urine ACR (100 x albumin / creatinine) and eGFR (CKD-EPI
#               2021, race-free, from calibrated creatinine, age and sex).
#
# METHODS       Supp. Section Zero, "Renal, body-size and lipid measures".
# =============================================================================
"""Derive urine ACR and race-free CKD-EPI 2021 eGFR.

ACR is available for every complete urine pair. eGFR uses the already verified
eligible demographics and assay-corrected serum creatinine tables, ensuring
that calibration, age, sex and eligibility each have one authoritative source.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
OUTPUT = ROOT / "outputs" / "03_derived"


def egfr_ckdepi_2021(creatinine, age, sex):
    """Return CKD-EPI 2021 creatinine eGFR; NHANES sex codes 1=male, 2=female."""
    creatinine = np.asarray(creatinine, dtype=float)
    female = np.asarray(sex) == 2
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    ratio = creatinine / kappa
    estimate = (142.0 * np.minimum(ratio, 1.0) ** alpha
                * np.maximum(ratio, 1.0) ** -1.200
                * 0.9938 ** np.asarray(age, dtype=float))
    return np.where(female, estimate * 1.012, estimate)


def main() -> None:
    urine = pd.read_csv(ROOT / "outputs" / "02_harmonised" / "urine_albumin_creatinine.csv")
    urine_derived = urine[["SEQN", "cycle"]].copy()
    urine_derived["Alb_Cr"] = urine.URXUMA / urine.URXUCR * 100
    urine_derived = urine_derived.replace([np.inf, -np.inf], np.nan)

    eligible = pd.read_csv(ROOT / "outputs" / "04_cohorts" / "eligible_adult_cohort.csv")
    biochemistry = pd.read_csv(ROOT / "outputs" / "02_harmonised" / "standard_biochemistry.csv",
                              usecols=["SEQN", "LBXSCR"])
    kidney = eligible[["SEQN", "RIAGENDR", "RIDAGEYR"]].merge(
        biochemistry, on="SEQN", how="inner").merge(
        urine_derived[["SEQN", "Alb_Cr"]], on="SEQN", how="inner")
    kidney["ACR"] = kidney.pop("Alb_Cr")
    kidney["eGFR"] = egfr_ckdepi_2021(kidney.LBXSCR, kidney.RIDAGEYR, kidney.RIAGENDR)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    urine_derived.to_csv(OUTPUT / "urine_derived_features.csv", index=False)
    kidney.to_csv(OUTPUT / "kidney_features_eligible.csv", index=False)
    print(f"kidney derivation: {len(urine_derived):,} ACR values; {len(kidney):,} eligible participants with ACR and eGFR")


if __name__ == "__main__":
    main()
