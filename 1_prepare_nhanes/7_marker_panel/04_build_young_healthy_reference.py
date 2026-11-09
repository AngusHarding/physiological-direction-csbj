"""Independently construct strict and historical-compatible healthy references."""

from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
MASTER = ROOT / "outputs/07_master/analysis_master_raw.csv"
OUT = ROOT / "outputs/08_reference"
DISEASES = ["arthritis", "congestive_heart_failure", "coronary_heart_disease",
            "angina", "heart_attack", "stroke", "emphysema", "thyroid_disease",
            "chronic_bronchitis", "liver_condition", "cancer", "diabetes", "kidney_fail"]
PFQ = ["PFQ051", "PFQ054", "PFQ057"]


def three_state(condition_pass: pd.Series, condition_fail: pd.Series) -> pd.Series:
    """Encode an auditable pass/fail/unknown state without filling missing data."""
    return pd.Series(np.select([condition_pass, condition_fail], ["pass", "fail"],
                               default="unknown"), index=condition_pass.index)


def build_states(master: pd.DataFrame) -> pd.DataFrame:
    states = pd.DataFrame({"SEQN": master.SEQN, "RIAGENDR": master.RIAGENDR,
                           "RIDAGEYR": master.RIDAGEYR})
    states["age_state"] = three_state(master.RIDAGEYR.lt(30), master.RIDAGEYR.ge(30))
    states["bmi_state"] = three_state(master.BMI.between(18.5, 24.9),
                                       master.BMI.notna() & ~master.BMI.between(18.5, 24.9))

    disease = master[DISEASES]
    any_disease = disease.eq(1).any(axis=1)
    all_known_absent = disease.eq(0).all(axis=1)
    states["disease_state"] = three_state(all_known_absent, any_disease)
    states["disease_unknown_items"] = disease.isna().sum(axis=1)
    # Several MCQ disease items were not administered at ages 18-19 in later
    # cycles. Non-administration is not evidence of health, but excluding it
    # would erase those ages by design. Preserve it as an explicit N/A pass.
    disease_structural_na = master.RIDAGEYR.lt(20) & ~any_disease
    states["disease_available_state"] = three_state(
        all_known_absent | disease_structural_na, any_disease)
    states["disease_structural_not_applicable"] = disease_structural_na & ~all_known_absent

    states["smoking_state"] = three_state(master.Smoking_status.eq("Never"),
        master.Smoking_status.notna() & ~master.Smoking_status.eq("Never"))
    states["alcohol_state"] = three_state(master.drinks_per_week.le(7),
        master.drinks_per_week.notna() & master.drinks_per_week.gt(7))

    pfq_valid_no = master[PFQ].eq(2).all(axis=1)
    pfq_any_yes = master[PFQ].eq(1).any(axis=1)
    under_20 = master.RIDAGEYR.lt(20)
    states["functional_state"] = three_state(under_20 | (~under_20 & pfq_valid_no),
                                               ~under_20 & pfq_any_yes)
    pfq_unavailable_cycle = master.cycle.isin(["1999-2000", "2001-2002"])
    states["functional_available_state"] = three_state(
        under_20 | pfq_unavailable_cycle | pfq_valid_no,
        ~under_20 & ~pfq_unavailable_cycle & pfq_any_yes)
    states["functional_not_applicable"] = under_20
    states["functional_structural_not_collected"] = pfq_unavailable_cycle
    states["functional_unknown_items"] = master[PFQ].isna().sum(axis=1)
    return states


def cascade(states: pd.DataFrame, reference_column: str, mode: str) -> pd.DataFrame:
    rows = []
    mask = pd.Series(True, index=states.index)
    disease_column = "disease_available_state" if mode == "availability_aware" else "disease_state"
    functional_column = ("functional_available_state" if mode == "availability_aware"
                         else "functional_state")
    for label, column in (("master", None), ("age_under_30", "age_state"),
                          ("healthy_BMI", "bmi_state"), ("disease_free", disease_column),
                          ("never_smoker", "smoking_state"), ("alcohol_le_7", "alcohol_state"),
                          ("no_functional_impairment", functional_column)):
        if column:
            if mode == "strict":
                mask &= states[column].eq("pass")
            else:
                # Historical code was permissive only for missing disease and
                # functional items. BMI, smoking and alcohol missingness failed.
                if column in {"disease_state", "functional_state"}:
                    mask &= ~states[column].eq("fail")
                else:
                    mask &= states[column].eq("pass")
        rows.append({"definition": mode, "step": label, "total": int(mask.sum()),
                     "men": int((mask & states.RIAGENDR.eq(1)).sum()),
                     "women": int((mask & states.RIAGENDR.eq(2)).sum())})
    states[reference_column] = mask
    return pd.DataFrame(rows)


def first_exclusion(states: pd.DataFrame, mode: str) -> pd.Series:
    disease_column = "disease_available_state" if mode == "availability_aware" else "disease_state"
    functional_column = ("functional_available_state" if mode == "availability_aware"
                         else "functional_state")
    order = [("age", "age_state"), ("bmi", "bmi_state"), ("disease", disease_column),
             ("smoking", "smoking_state"), ("alcohol", "alcohol_state"),
             ("functional", functional_column)]
    result = pd.Series("reference", index=states.index, dtype=object)
    active = pd.Series(True, index=states.index)
    for label, column in order:
        allowed = states[column].eq("pass")
        if mode == "historical_compatible" and column in {"disease_state", "functional_state"}:
            allowed = ~states[column].eq("fail")
        excluded = active & ~allowed
        result.loc[excluded] = label + "_" + states.loc[excluded, column]
        active &= allowed
    return result


def main() -> None:
    master = pd.read_csv(MASTER)
    if len(master) != 29053 or master.SEQN.duplicated().any():
        raise AssertionError("master cohort lock failed")
    states = build_states(master)
    strict_cascade = cascade(states, "is_reference_strict", "strict")
    primary_cascade = cascade(states, "is_reference_primary", "availability_aware")
    historical_cascade = cascade(states, "is_reference_historical_compatible",
                                  "historical_compatible")
    states["strict_first_exclusion"] = first_exclusion(states, "strict")
    states["primary_first_exclusion"] = first_exclusion(states, "availability_aware")
    states["historical_first_exclusion"] = first_exclusion(states, "historical_compatible")

    criterion_rows = []
    for column in [c for c in states.columns if c.endswith("_state")]:
        for sex, label in ((None, "all"), (1, "men"), (2, "women")):
            subset = states if sex is None else states.loc[states.RIAGENDR.eq(sex)]
            counts = subset[column].value_counts()
            criterion_rows.append({"criterion": column.removesuffix("_state"), "sex": label,
                                   "pass": int(counts.get("pass", 0)),
                                   "fail": int(counts.get("fail", 0)),
                                   "unknown": int(counts.get("unknown", 0))})

    OUT.mkdir(parents=True, exist_ok=True)
    states.to_csv(OUT / "healthy_reference_flags.csv", index=False)
    states[["SEQN", "RIAGENDR", "RIDAGEYR", "is_reference_primary"]].rename(
        columns={"is_reference_primary": "is_reference"}).to_csv(
        OUT / "healthy_reference_members.csv", index=False)
    pd.concat([primary_cascade, strict_cascade, historical_cascade], ignore_index=True).to_csv(
        OUT / "healthy_reference_cascade.csv", index=False)
    pd.DataFrame(criterion_rows).to_csv(OUT / "healthy_reference_criterion_counts.csv", index=False)

    print(pd.concat([primary_cascade.tail(1), strict_cascade.tail(1),
                     historical_cascade.tail(1)]).to_string(index=False))
    # public_code: the original code ended by comparing the historical-compatible
    # members with the stored list from the first analysis. That file is not part of
    # the public release; on 6 Nov 2026 the comparison found exact agreement.


if __name__ == "__main__":
    main()
