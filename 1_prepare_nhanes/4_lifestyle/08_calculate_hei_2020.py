# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/4_lifestyle/08_calculate_hei_2020.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Computes the Healthy Eating Index 2020 for each participant,
#               following heiscore 0.1.4.
#
# METHODS       Table 4 and Results 6 ("the lowest third of diet quality").
# =============================================================================
"""Reproduce individual-level simple HEI-2020 scores component by component.

The implementation follows heiscore 0.1.4 exactly, including its two-day-weight
eligibility, day-one plus day-two energy-density calculation, and linear component
scoring. It retains substantially more provenance than the historical total-only
file. This stage replicates the frozen package; official FPED validation is separate.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
COMPONENTS = ["F_TOTAL", "FWHOLEFRT", "VTOTALLEG", "VDRKGRLEG", "G_WHOLE",
              "D_TOTAL", "PFALLPROTLEG", "PFSEAPLANTLEG", "TFACIDS",
              "G_REFINED", "TSODI", "ADD_SUGARS", "TSFAT"]


def source_columns(frame: pd.DataFrame, component: str) -> list[str]:
    if component == "TFACIDS":
        return ["TOT_TFACIDS"]
    return [column for column in frame.columns if component in column]


def score_ratio(ratio: pd.Series, standard: pd.Series) -> pd.Series:
    maximum, zero, points = standard.max_amount, standard.zero_score, standard.max_points
    if standard.component_type == "adequacy":
        score = (ratio - zero).div(maximum - zero).mul(points)
        return score.clip(lower=0, upper=points)
    score = (zero - ratio).div(zero - maximum).mul(points)
    return score.clip(lower=0, upper=points)


def main() -> None:
    raw = pd.read_csv(ROOT / "outputs/01_imported/hei_package_fped_inputs.csv")
    standards = pd.read_csv(ROOT / "config/hei_2020_scoring_standards.csv").set_index("component")

    # heiscore::cleanDataset requires a nonzero two-day dietary weight and age >=1;
    # the historical call's default demographic filter further restricts age 2-100.
    eligible = (raw.WTDR2D.notna() & raw.WTDR2D.ne(0) & raw.AGE.between(2, 100) &
                raw.SEX.isin(["Female", "Male"]) &
                raw.RACE_ETH.isin(["Asian", "White", "Black", "Other", "Mexican American", "Other Hispanic"]))
    work = raw.loc[eligible].copy()
    # The archived source uses ordinary addition here: if either recall day's fat
    # fields are absent, the fatty-acid ratio is unresolved. The historical total
    # then treated that missing component as a zero contribution (eight children,
    # no adults in this study); both facts are retained explicitly below.
    work["TOT_TFACIDS"] = (work.DR1_MONOPOLY + work.DR2_MONOPOLY).div(work.DR1TSFAT + work.DR2TSFAT)
    energy = work[["DR1TKCAL", "DR2TKCAL"]].sum(axis=1, skipna=True)

    result = work[["SEQN", "cycle", "cycle_code", "WTDRD1", "WTDR2D", "SEX", "AGE", "RACE_ETH", "FAMINC"]].copy()
    result["one_recall_day_energy_missing"] = work[["DR1TKCAL", "DR2TKCAL"]].isna().sum(axis=1).eq(1)
    for component in COMPONENTS:
        columns = source_columns(work, component)
        amount = work[columns].sum(axis=1, skipna=True) if component != "TFACIDS" else work.TOT_TFACIDS
        if component == "TSODI":
            ratio = amount.div(energy)
        elif component == "TSFAT":
            ratio = amount.mul(9).div(energy).mul(100)
        elif component == "ADD_SUGARS":
            ratio = amount.mul(16).div(energy).mul(100)
        elif component == "TFACIDS":
            ratio = amount
        else:
            ratio = amount.div(energy).mul(1000)
        result[f"amount_{component}"] = amount
        result[f"ratio_{component}"] = ratio.replace([np.inf, -np.inf], np.nan)
        result[f"score_{component}"] = score_ratio(result[f"ratio_{component}"], standards.loc[component])

    score_columns = [f"score_{component}" for component in COMPONENTS]
    result["missing_component_count"] = result[score_columns].isna().sum(axis=1)
    result["HEI_2020_total"] = result[score_columns].sum(axis=1, skipna=True)
    result["all_components_scored"] = result[score_columns].notna().all(axis=1)
    # Historical package-compatible total. The strict total is retained as the
    # preferable general-purpose definition, requiring every component.
    result["HEI_2020_total_strict"] = result[score_columns].sum(axis=1, min_count=len(score_columns))
    # Package output excludes rows where no component can be scored at all, while
    # retaining the eight partial rows described above.
    final = result.loc[result.missing_component_count.lt(len(score_columns))].copy()

    out = ROOT / "outputs/02_harmonised"
    out.mkdir(parents=True, exist_ok=True)
    final.to_csv(out / "hei_2020_component_scores.csv", index=False)
    final[["SEQN", "cycle", "HEI_2020_total", "HEI_2020_total_strict",
           "missing_component_count"]].to_csv(out / "hei_2020_total.csv", index=False)

    audit = raw[["SEQN", "cycle", "WTDR2D", "AGE", "SEX", "RACE_ETH"]].copy()
    audit["eligible_two_day_package_score"] = eligible
    complete_ids = result.loc[result.all_components_scored, "SEQN"]
    audit["scored_all_components"] = audit.SEQN.isin(complete_ids)
    audit["historical_total_available"] = audit.SEQN.isin(final.SEQN)
    audit["exclusion_reason"] = np.select(
        [raw.WTDR2D.isna() | raw.WTDR2D.eq(0), ~raw.AGE.between(2, 100),
         ~raw.SEX.isin(["Female", "Male"]),
         ~raw.RACE_ETH.isin(["Asian", "White", "Black", "Other", "Mexican American", "Other Hispanic"]),
         eligible & ~audit.scored_all_components],
        ["two_day_weight_missing_or_zero", "age_outside_2_100", "sex_outside_default",
         "race_ethnicity_outside_default", "one_or_more_components_unscored"],
        default="scored")
    audit.to_csv(out / "hei_2020_scoring_audit.csv", index=False)
    final.groupby("cycle").HEI_2020_total.agg(["count", "mean", "std", "min", "median", "max"]).reset_index().to_csv(
        out / "hei_2020_by_cycle.csv", index=False)
    print(f"HEI-2020: {len(final):,}/{len(raw):,} participants scored")


if __name__ == "__main__":
    main()
