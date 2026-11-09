"""Harmonise doctor-diagnosed condition responses to 1=yes, 0=no, NA=unknown.

No missing response is converted to absence.  Diabetes code 3 (borderline) is
left unresolved, matching the manuscript-era rule.  Asked/observed indicators
make intermittent questionnaire coverage auditable downstream.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)

VARIABLES = {
    "asthma": ["MCQ010"],
    "blood_transfusion": ["MCQ092"],
    "arthritis": ["MCQ160A"],
    "congestive_heart_failure": ["MCQ160B"],
    "coronary_heart_disease": ["MCQ160C"],
    "angina": ["MCQ160D"],
    "heart_attack": ["MCQ160E"],
    "stroke": ["MCQ160F"],
    "emphysema": ["MCQ160G"],
    "thyroid_disease": ["MCQ160I", "MCD160M", "MCQ160M"],
    "chronic_bronchitis": ["MCQ160K"],
    "liver_condition": ["MCQ160L"],
    "gout": ["MCQ160N"],
    "copd": ["MCQ160O"],
    "cancer": ["MCQ220"],
    "diabetes": ["DIQ010"],
    "kidney_fail": ["KIQ020", "KIQ022"],
    "hypertension": ["BPQ020"],
    "osteoporosis": ["OSQ060"],
}

PRIMARY_13 = [
    "arthritis", "congestive_heart_failure", "coronary_heart_disease",
    "angina", "heart_attack", "stroke", "emphysema", "thyroid_disease",
    "chronic_bronchitis", "liver_condition", "cancer", "diabetes",
    "kidney_fail",
]


def coalesce(frame: pd.DataFrame, candidates: list[str]) -> pd.Series:
    """Combine renamed variables and reject overlapping contradictory answers."""
    present = [name for name in candidates if name in frame]
    if not present:
        return pd.Series(np.nan, index=frame.index)
    values = frame[present]
    contradictions = values.nunique(axis=1, dropna=True).gt(1)
    if contradictions.any():
        raise AssertionError(f"contradictory renamed fields: {candidates}")
    return values.bfill(axis=1).iloc[:, 0]


def binary(raw: pd.Series) -> pd.Series:
    return raw.map({1.0: 1.0, 2.0: 0.0})


def main() -> None:
    raw = pd.read_csv(ROOT / "outputs/01_imported/self_reported_disease_raw.csv")
    flags = raw[["SEQN", "cycle"]].copy()

    for condition, fields in VARIABLES.items():
        response = coalesce(raw, fields)
        flags[condition] = binary(response)
        flags[f"{condition}_observed"] = response.isin([1, 2])

    # A complete historical burden requires an observed yes/no for every item.
    observed_columns = [f"{name}_observed" for name in PRIMARY_13]
    flags["historical_13_complete"] = flags[observed_columns].all(axis=1)
    flags["historical_13_count"] = flags[PRIMARY_13].sum(axis=1, min_count=len(PRIMARY_13))
    flags["historical_13_any"] = flags["historical_13_count"].gt(0).where(
        flags.historical_13_complete
    )

    out = ROOT / "outputs/02_harmonised"
    out.mkdir(parents=True, exist_ok=True)
    flags.to_csv(out / "self_reported_disease_flags.csv", index=False)
    flags[["SEQN", "cycle"] + PRIMARY_13 + ["historical_13_complete",
          "historical_13_count", "historical_13_any"]].to_csv(
        out / "self_reported_disease_primary_13.csv", index=False
    )

    audit = []
    for condition in VARIABLES:
        for cycle, group in flags.groupby("cycle", sort=False):
            answered = group[f"{condition}_observed"]
            audit.append({
                "condition": condition, "cycle": cycle, "cycle_rows": len(group),
                "answered_n": int(answered.sum()),
                "yes_n": int(group.loc[answered, condition].eq(1).sum()),
                "prevalence_answered_pct": group.loc[answered, condition].mean() * 100
                    if answered.any() else np.nan,
            })
    pd.DataFrame(audit).to_csv(out / "self_reported_disease_by_cycle.csv", index=False)
    print(f"Disease harmonisation: {len(flags):,} participants")
    print(f"Complete historical 13: {flags.historical_13_complete.sum():,}")


if __name__ == "__main__":
    main()
