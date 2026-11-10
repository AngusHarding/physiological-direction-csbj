# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/4_lifestyle/06_harmonise_physical_activity.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Classifies physical activity within each questionnaire era.
#
# METHODS       Table 1 ("physically inactive", era-specific definitions) and
#               Table 4 (recreational activity, 2007-2014).
# =============================================================================
"""Harmonize physical activity while retaining the questionnaire-era distinction.

The all-cycle flag means broadly nonsedentary, not guideline-compliant exercise.
For 1999–2006 it comes from usual daily activity; for 2007–2018 it means any
reported GPAQ work, transport, or recreational MET dose. A separate ≥600
MET-min/week flag is available only in the GPAQ era.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
GATES = ["PAQ605", "PAQ620", "PAQ635", "PAQ650", "PAQ665"]
DOMAINS = [("PAQ610", "PAD615", 8, "vigorous_work"),
           ("PAQ625", "PAD630", 4, "moderate_work"),
           ("PAQ640", "PAD645", 4, "transport"),
           ("PAQ655", "PAD660", 8, "vigorous_recreation"),
           ("PAQ670", "PAD675", 4, "moderate_recreation")]
HEALTHY_MET = 600
MAX_MET = 10_000


def derive_early(raw: pd.DataFrame) -> pd.DataFrame:
    early = raw.loc[raw.instrument_era.eq("PAQ180")].copy()
    early["physically_active"] = early.PAQ180.map({1: 0, 2: 1, 3: 1, 4: 1})
    early["met_minutes_week"] = np.nan
    early["healthy_exercise"] = np.nan
    early["gates_valid"] = np.nan
    early["met_above_audit_limit"] = False
    return early


def derive_gpaq(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    original = raw.loc[raw.instrument_era.eq("GPAQ")].copy()
    original["gates_valid"] = np.logical_and.reduce([original[g].le(2) for g in GATES])
    # Historical output retains only respondents with valid yes/no answers to all
    # five domain gates. The prefilter audit remains available separately.
    gpaq = original.loc[original.gates_valid].copy()
    total = pd.Series(0.0, index=gpaq.index)
    for days_col, minutes_col, mets, name in DOMAINS:
        days = gpaq[days_col].where(gpaq[days_col].between(1, 7))
        minutes = gpaq[minutes_col].where(gpaq[minutes_col].between(1, 1440))
        component = days.fillna(0) * minutes.fillna(0) * mets
        gpaq[f"met_{name}"] = component
        total = total + component
    gpaq["met_minutes_week"] = total
    gpaq["met_above_audit_limit"] = gpaq.met_minutes_week.gt(MAX_MET)
    above_limit_n = int(gpaq.met_above_audit_limit.sum())
    # Inclusion-preserving alternative: retain a respondent's clear activity
    # classification while limiting only the continuous dose used downstream.
    corrected = gpaq.copy()
    corrected["met_minutes_week_raw"] = corrected.met_minutes_week
    corrected["met_minutes_week"] = corrected.met_minutes_week.clip(upper=MAX_MET)
    corrected["healthy_exercise"] = corrected.met_minutes_week_raw.ge(HEALTHY_MET).astype(int)
    corrected["physically_active"] = corrected.met_minutes_week_raw.gt(0).astype(int)
    corrected["met_extreme_above_10000_flag"] = corrected.met_minutes_week_raw.gt(MAX_MET).astype("Int64")
    corrected["met_quality_flag"] = np.where(
        corrected.met_extreme_above_10000_flag.eq(1),
        "extreme_self_report_review_recommended",
        "within_prespecified_range",
    )
    gpaq = gpaq.loc[~gpaq.met_above_audit_limit].copy()
    gpaq["healthy_exercise"] = gpaq.met_minutes_week.ge(HEALTHY_MET).astype(int)
    gpaq["physically_active"] = gpaq.met_minutes_week.gt(0).astype(int)
    return gpaq, corrected, original, above_limit_n


def main() -> None:
    raw = pd.read_csv(ROOT / "outputs/01_imported/physical_activity_raw.csv")
    early = derive_early(raw)
    gpaq, gpaq_corrected, gpaq_prefilter, above_limit_n = derive_gpaq(raw)
    columns = ["SEQN", "cycle", "physically_active", "met_minutes_week", "healthy_exercise"]
    clean = pd.concat([early[columns], gpaq[columns]], ignore_index=True)
    # Early-cycle participants have no comparable MET dose, so their warning fields
    # are explicitly not applicable rather than False.
    early_corrected = early[columns].copy()
    early_corrected["met_minutes_week_raw"] = np.nan
    early_corrected["met_extreme_above_10000_flag"] = pd.Series(pd.NA, index=early_corrected.index, dtype="Int64")
    early_corrected["met_quality_flag"] = "not_applicable_PAQ180_era"
    corrected_columns = columns + ["met_minutes_week_raw", "met_extreme_above_10000_flag", "met_quality_flag"]
    corrected = pd.concat([early_corrected[corrected_columns], gpaq_corrected[corrected_columns]], ignore_index=True)

    out = ROOT / "outputs/02_harmonised"
    out.mkdir(parents=True, exist_ok=True)
    clean.to_csv(out / "physical_activity.csv", index=False)
    corrected.to_csv(out / "physical_activity_inclusion_preserving.csv", index=False)
    early.to_csv(out / "physical_activity_early_audit.csv", index=False)
    gpaq.to_csv(out / "physical_activity_gpaq_audit.csv", index=False)
    gpaq_prefilter[["SEQN", "cycle", "gates_valid"]].to_csv(out / "physical_activity_gpaq_gate_audit.csv", index=False)

    extreme = gpaq_corrected.loc[gpaq_corrected.met_extreme_above_10000_flag.eq(1), [
        "SEQN", "cycle", "met_minutes_week_raw", "met_minutes_week",
        "met_extreme_above_10000_flag", "met_quality_flag",
    ]].copy()
    # The GPAQ calculation uses 4 MET for moderate and 8 MET for vigorous activity.
    # These translations are equivalent-duration illustrations, not claims about
    # the respondent's actual mix of activity domains.
    extreme["equivalent_moderate_hours_week_at_4_MET"] = extreme.met_minutes_week_raw / 4 / 60
    extreme["equivalent_vigorous_hours_week_at_8_MET"] = extreme.met_minutes_week_raw / 8 / 60
    extreme.to_csv(out / "physical_activity_extreme_met_audit.csv", index=False)
    extreme[["met_minutes_week_raw", "equivalent_moderate_hours_week_at_4_MET",
             "equivalent_vigorous_hours_week_at_8_MET"]].describe(
                 percentiles=[.01, .25, .5, .75, .9, .95, .99]
             ).T.reset_index(names="measure").to_csv(
                 out / "physical_activity_extreme_met_summary.csv", index=False
             )

    summary = clean.groupby("cycle").agg(
        output_n=("SEQN", "size"), observed_primary_n=("physically_active", "count"),
        active_n=("physically_active", lambda x: x.eq(1).sum()),
        sedentary_n=("physically_active", lambda x: x.eq(0).sum()),
        met_observed_n=("met_minutes_week", "count"),
        median_met=("met_minutes_week", "median"),
        guideline_active_n=("healthy_exercise", lambda x: x.eq(1).sum()),
    ).reset_index()
    summary.to_csv(out / "physical_activity_by_cycle.csv", index=False)
    pd.DataFrame([{
        "raw_participants": len(raw), "early_raw_n": int(raw.instrument_era.eq("PAQ180").sum()),
        "gpaq_raw_n": len(gpaq_prefilter), "gpaq_invalid_gate_n": int((~gpaq_prefilter.gates_valid).sum()),
        "gpaq_above_10000_met_n": above_limit_n,
        "harmonised_rows": len(clean), "primary_observed_n": int(clean.physically_active.notna().sum()),
    }]).to_csv(out / "physical_activity_change_log.csv", index=False)
    print(f"Physical activity: {len(clean):,} harmonized rows; "
          f"{clean.physically_active.notna().sum():,} with primary flag")


if __name__ == "__main__":
    main()
