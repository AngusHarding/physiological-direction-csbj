# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/4_lifestyle/04_harmonise_alcohol.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Converts the alcohol answers to drinks per week (screened
#               non-drinkers and zero frequency = 0; refused or unknown stay
#               missing).
#
# METHODS       "Young-healthy reference" ("drank lightly or not at all";
#               Supp. Section Zero: "no more than seven alcoholic drinks per week"); Table 1.
# =============================================================================
"""Harmonize past-year average alcohol intake as drinks per week.

Older cycles record a frequency count and unit; 2017-2018 records a frequency
category. Screened lifetime/past-year nondrinkers and zero-frequency respondents
receive zero. Refused, unknown, or insufficient responses remain missing.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)

UNIT_TO_YEAR = {1: 52.0, 2: 12.0, 3: 1.0}
ALQ121_TO_DAYS = {0: 0.0, 1: 365.0, 2: 350.0, 3: 182.0, 4: 104.0,
                  5: 52.0, 6: 30.0, 7: 12.0, 8: 9.0, 9: 4.5, 10: 1.5}


def clean_quantity(series: pd.Series, cap: bool = True) -> pd.Series:
    """Retain plausible numeric answers; 777/999-style codes become missing."""
    result = series.where(series < 100)
    if cap:
        result = result.clip(upper=15)
    return result


def derive(frame: pd.DataFrame, cap_quantity: bool = True) -> pd.DataFrame:
    result = frame[["SEQN", "cycle", "instrument_era"]].copy()
    days = pd.Series(np.nan, index=frame.index, dtype=float)
    quantity = clean_quantity(frame.ALQ130, cap=cap_quantity)
    nondrinker = pd.Series(False, index=frame.index)

    old = frame.instrument_era.eq("count_unit")
    unit = frame.ALQ120U.map(UNIT_TO_YEAR)
    valid_count = frame.ALQ120Q.between(0, 365)
    days.loc[old & frame.ALQ120Q.eq(0)] = 0.0
    convert = old & valid_count & unit.notna() & ~frame.ALQ120Q.eq(0)
    converted_days = pd.to_numeric(frame.loc[convert, "ALQ120Q"], errors="coerce").astype(float) * unit.loc[convert]
    days.loc[convert] = converted_days.astype(float)
    nondrinker.loc[old] = (frame.loc[old, "ALQ101"].eq(2) |
                            frame.loc[old, "ALQ110"].eq(2)).fillna(False)

    new = frame.instrument_era.eq("categorical_frequency")
    days.loc[new] = frame.loc[new, "ALQ121"].map(ALQ121_TO_DAYS)
    nondrinker.loc[new] = frame.loc[new, "ALQ111"].eq(2).fillna(False)

    # No respondent can drink on more than 365 distinct days in one year.
    drinks = days.clip(upper=365).div(52).mul(quantity)
    drinks = drinks.mask(days.eq(0), 0.0)
    drinks = drinks.mask(nondrinker, 0.0)
    result["days_per_year"] = days
    result["drinks_per_drinking_day"] = quantity
    result["screened_nondrinker"] = nondrinker
    result["drinks_per_week"] = drinks
    return result


def main() -> None:
    raw = pd.read_csv(ROOT / "outputs/01_imported/alcohol_raw.csv")
    audit = derive(raw, cap_quantity=True)
    uncapped = derive(raw, cap_quantity=False)
    audit["drinks_per_week_uncapped"] = uncapped.drinks_per_week
    audit["quantity_capped"] = (raw.ALQ130.gt(15) & raw.ALQ130.lt(100)).fillna(False)
    audit["alcohol_observed"] = audit.drinks_per_week.notna()
    clean = audit.loc[audit.alcohol_observed, ["SEQN", "cycle", "drinks_per_week"]].copy()

    out = ROOT / "outputs/02_harmonised"
    out.mkdir(parents=True, exist_ok=True)
    clean.to_csv(out / "alcohol_drinks_per_week.csv", index=False)
    audit.to_csv(out / "alcohol_derivation_audit.csv", index=False)
    audit.groupby("cycle").agg(
        imported_n=("SEQN", "size"), observed_n=("alcohol_observed", "sum"),
        nondrinker_n=("screened_nondrinker", "sum"), zero_n=("drinks_per_week", lambda x: x.eq(0).sum()),
        capped_quantity_n=("quantity_capped", "sum"), median_drinks_week=("drinks_per_week", "median"),
        p99_drinks_week=("drinks_per_week", lambda x: x.quantile(.99)),
    ).reset_index().to_csv(out / "alcohol_by_cycle.csv", index=False)
    pd.DataFrame([{
        "rule": "cap average drinks per drinking day at 15",
        "participants_affected": int(audit.quantity_capped.sum()),
        "maximum_capped_drinks_per_week": float(audit.drinks_per_week.max()),
        "maximum_uncapped_drinks_per_week": float(audit.drinks_per_week_uncapped.max()),
    }]).to_csv(out / "alcohol_cap_sensitivity.csv", index=False)
    print(f"Alcohol harmonisation: {len(clean):,}/{len(raw):,} scored; "
          f"zero={(clean.drinks_per_week == 0).mean():.1%}; "
          f"quantity cap affected {int(audit.quantity_capped.sum()):,}")


if __name__ == "__main__":
    main()
