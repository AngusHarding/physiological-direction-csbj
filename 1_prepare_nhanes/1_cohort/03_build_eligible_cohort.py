"""Build explicit eligibility flags and reproduce the certified adult cohort.

This is the heart of Chunk 1.  It does not merely save the participants who
pass.  It also saves one row per demographic participant with individual True/
False flags explaining each rule.  That larger table lets a reviewer ask why a
particular SEQN was excluded.

Important pandas ideas used below
---------------------------------
``merge`` joins two tables using a shared key, like a database join or VLOOKUP.
``True``/``False`` columns are masks: they select rows satisfying a condition.
``axis=1`` means 'work across columns within each row'.
"""

from pathlib import Path

import numpy as np
import pandas as pd


import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
IMPORTED = ROOT / "outputs" / "01_imported"
OUTPUT = ROOT / "outputs" / "04_cohorts"


def primary_exclusion_reason(frame: pd.DataFrame) -> pd.Series:
    """Assign one deterministic first-failing reason for human audit."""
    # Conditions are deliberately ordered.  If somebody fails several rules,
    # np.select reports the first failure in this list as the primary reason.
    conditions = [
        ~frame["eligible_demographic_complete"],
        ~frame["eligible_recent_illness_matched"],
        ~frame["eligible_recent_illness_complete"],
        ~frame["eligible_age_18_79"],
        ~frame["eligible_not_acute"],
        ~frame["eligible_not_pregnant"],
    ]
    labels = [
        "missing_required_demographic",
        "no_recent_illness_record",
        "missing_recent_illness_response",
        "age_outside_18_79",
        "recent_acute_illness",
        "pregnant",
    ]
    return pd.Series(np.select(conditions, labels, default="eligible"), index=frame.index)


def main() -> None:
    # Load the two minimally interpreted tables produced by the import scripts.
    demo = pd.read_csv(IMPORTED / "demographics.csv")
    illness = pd.read_csv(IMPORTED / "recent_illness.csv")
    # Both tables contain 'cycle'. Rename the HSQ copy before merging so pandas
    # does not create confusing cycle_x and cycle_y column names.
    illness = illness.rename(columns={"cycle": "illness_cycle"})
    # how='left' keeps every demographic participant. validate='one_to_one'
    # makes pandas fail if either side contains duplicate SEQN values.
    frame = demo.merge(illness, on="SEQN", how="left", validate="one_to_one")

    # Each expression below creates a Boolean (True/False) column.  Keeping
    # rules separate is more auditable than one long filter expression.
    frame["eligible_demographic_complete"] = frame[["SEQN", "RIAGENDR", "RIDAGEYR"]].notna().all(axis=1)
    # A non-missing illness_cycle proves that SEQN appeared in the HSQ table.
    frame["eligible_recent_illness_matched"] = frame["illness_cycle"].notna()
    frame["eligible_recent_illness_complete"] = frame[["HSQ500", "HSQ510", "HSQ520"]].notna().all(axis=1)
    # ge = greater-than-or-equal; lt = strictly less than.
    frame["eligible_age_18_79"] = frame["RIDAGEYR"].ge(18) & frame["RIDAGEYR"].lt(80)
    # NHANES code 2 means 'No'. All three cells must equal 2.
    frame["eligible_not_acute"] = frame[["HSQ500", "HSQ510", "HSQ520"]].eq(2).all(axis=1)
    # This deliberately reproduces the certified rule. Structural RIDEXPRG missingness
    # is not pregnancy; only the affirmative pregnancy code is excluded.
    # ne(1) means 'not equal to the affirmative pregnancy code 1'.  Structural
    # missingness (for participants not asked) is therefore not called pregnant.
    frame["eligible_not_pregnant"] = frame["RIDEXPRG"].ne(1)
    required_flags = [
        "eligible_demographic_complete",
        "eligible_recent_illness_matched",
        "eligible_recent_illness_complete",
        "eligible_age_18_79",
        "eligible_not_acute",
        "eligible_not_pregnant",
    ]
    # Final eligibility is True only when every rule column is True.
    frame["eligible_adult"] = frame[required_flags].all(axis=1)
    frame["primary_exclusion_reason"] = primary_exclusion_reason(frame)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    # Save the full audit table first, including excluded participants.
    frame.to_csv(OUTPUT / "eligibility_flags.csv", index=False)
    # Then materialise the compact scientific cohort used downstream.
    eligible = frame.loc[frame.eligible_adult, ["SEQN", "RIAGENDR", "RIDAGEYR"]].reset_index(drop=True)
    eligible.to_csv(OUTPUT / "eligible_adult_cohort.csv", index=False)

    steps = [
        ("raw_demographics", pd.Series(True, index=frame.index)),
        ("demographic_complete", frame.eligible_demographic_complete),
        ("matched_recent_illness", frame.eligible_recent_illness_matched),
        ("complete_recent_illness", frame.eligible_recent_illness_complete),
        ("age_18_79", frame.eligible_age_18_79),
        ("not_acute", frame.eligible_not_acute),
        ("not_pregnant", frame.eligible_not_pregnant),
    ]
    # mask starts with everybody. At each loop it is narrowed using &=, meaning
    # 'passed every previous rule AND this rule'. This creates a true cascade.
    mask = pd.Series(True, index=frame.index)
    cascade = []
    for step, condition in steps:
        mask &= condition
        cascade.append(
            {
                "step": step,
                "total": int(mask.sum()),
                "men": int((mask & frame.RIAGENDR.eq(1)).sum()),
                "women": int((mask & frame.RIAGENDR.eq(2)).sum()),
                "removed_at_step": int((~condition & mask.shift(fill_value=True)).sum()) if False else "",
            }
        )
    cascade_frame = pd.DataFrame(cascade)
    cascade_frame["removed_from_previous"] = cascade_frame.total.shift(fill_value=cascade_frame.total.iloc[0]) - cascade_frame.total
    cascade_frame = cascade_frame.drop(columns="removed_at_step")
    cascade_frame.to_csv(OUTPUT / "eligibility_cascade.csv", index=False)

    # A separate summary makes cycle and sex imbalances visible without opening
    # the large participant-level flags file.
    grouped = (
        frame.groupby(["cycle", "RIAGENDR"], dropna=False)
        .agg(total=("SEQN", "size"), eligible=("eligible_adult", "sum"))
        .reset_index()
    )
    grouped["sex"] = grouped.RIAGENDR.map({1.0: "men", 2.0: "women"}).fillna("unknown")
    grouped[["cycle", "sex", "total", "eligible"]].to_csv(
        OUTPUT / "eligibility_by_cycle_and_sex.csv", index=False
    )
    print(
        f"eligible adult cohort: {len(eligible):,} "
        f"(men={(eligible.RIAGENDR == 1).sum():,}, women={(eligible.RIAGENDR == 2).sum():,})"
    )


if __name__ == "__main__":
    main()
