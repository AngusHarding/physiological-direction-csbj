"""Derive mutually exclusive Never, Former, and Current smoking status.

NHANES asks SMQ040 only of participants who report at least 100 lifetime
cigarettes. A missing SMQ040 is therefore structurally correct for never-smokers,
not missing evidence. Refused, unknown, inconsistent, and genuinely missing
responses remain unclassified rather than being guessed.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)


def classify(frame: pd.DataFrame) -> pd.Series:
    status = pd.Series(pd.NA, index=frame.index, dtype="string")
    status.loc[frame.SMQ020.eq(2)] = "Never"
    status.loc[frame.SMQ020.eq(1) & frame.SMQ040.eq(3)] = "Former"
    status.loc[frame.SMQ020.eq(1) & frame.SMQ040.isin([1, 2])] = "Current"
    return status


def exclusion_reason(frame: pd.DataFrame, status: pd.Series) -> pd.Series:
    reason = pd.Series("classified", index=frame.index, dtype="string")
    unresolved = status.isna()
    reason.loc[unresolved & frame.SMQ020.isin([7, 9])] = "SMQ020 refused_or_unknown"
    reason.loc[unresolved & frame.SMQ020.isna()] = "SMQ020 missing"
    reason.loc[unresolved & frame.SMQ020.eq(1) & frame.SMQ040.isin([7, 9])] = "SMQ040 refused_or_unknown"
    reason.loc[unresolved & frame.SMQ020.eq(1) & frame.SMQ040.isna()] = "SMQ040 missing_after_ever_smoker"
    reason.loc[unresolved & reason.eq("classified")] = "unrecognised_or_inconsistent_code"
    return reason


def main() -> None:
    raw = pd.read_csv(ROOT / "outputs/01_imported/smoking_raw.csv")
    status = classify(raw)
    reason = exclusion_reason(raw, status)

    audit = raw.copy()
    audit["Smoking_status"] = status
    audit["classification_reason"] = reason
    audit["smoking_status_observed"] = status.notna()

    clean = audit.loc[audit.smoking_status_observed, ["SEQN", "cycle", "Smoking_status"]].copy()
    clean["Smoking_status"] = pd.Categorical(
        clean.Smoking_status, categories=["Never", "Former", "Current"], ordered=True
    )

    out = ROOT / "outputs/02_harmonised"
    out.mkdir(parents=True, exist_ok=True)
    clean.to_csv(out / "smoking_status.csv", index=False)
    audit.to_csv(out / "smoking_classification_audit.csv", index=False)

    by_cycle = (
        audit.groupby(["cycle", "Smoking_status"], dropna=False, observed=False)
        .size().rename("n").reset_index()
    )
    by_cycle.to_csv(out / "smoking_by_cycle.csv", index=False)
    reason.value_counts(dropna=False).rename_axis("classification_reason").reset_index(name="n").to_csv(
        out / "smoking_change_log.csv", index=False
    )
    print(f"Smoking harmonisation: {len(clean):,}/{len(raw):,} classified")
    print(clean.Smoking_status.value_counts().reindex(["Never", "Former", "Current"]).to_string())


if __name__ == "__main__":
    main()

