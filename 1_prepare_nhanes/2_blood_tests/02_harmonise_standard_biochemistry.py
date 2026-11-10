# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/2_blood_tests/02_harmonise_standard_biochemistry.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Renames biochemistry variables to one name per marker and
#               applies the assay corrections: creatinine calibration (1999-2000,
#               2005-2006) and the 2017-2018 Roche-to-Beckman bridging equations;
#               nonpositive values set to missing.
#
# METHODS       Supp. Section Zero, "Assay calibration" and "Cleaning" (Table S0.3).
# =============================================================================
"""Apply frozen biochemistry renames, assay corrections and complete-casing."""

from pathlib import Path

import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
RAW = ROOT / "outputs" / "01_imported" / "standard_biochemistry_raw.csv"
OUTPUT = ROOT / "outputs" / "02_harmonised"
MARKERS = [
    "LBXSAL", "LBXSATSI", "LBXSASSI", "LBXSAPSI", "LBXSBU", "LBXSCA", "LBXSC3SI",
    "LBXSCR", "LBXSGTSI", "LBXSIR", "LBXSLDSI", "LBXSPH", "LBXSTB", "LBXSTP",
    "LBXSUA", "LBXSNASI", "LBXSKSI", "LBXSCLSI", "LBXSOSSI", "LBXSGB", "LBXSCH",
    "LBXSTR", "LBXSGL",
]
ALIASES = {"LBDSAPSI": "LBXSAPSI", "LBDSLDSI": "LBXSLDSI", "LBDSCR": "LBXSCR",
           "LBDSPH": "LBXSPH", "LBDSTB": "LBXSTB"}


def apply_rule(frame, cycle, marker, rule_id, function, equation, log):
    mask = frame.cycle.eq(cycle) & frame[marker].notna()
    before = frame.loc[mask, marker].copy()
    frame.loc[mask, marker] = function(before)
    after = frame.loc[mask, marker]
    log.append({"rule_id": rule_id, "cycle": cycle, "marker": marker, "equation": equation,
                "values_transformed": int(mask.sum()), "before_min": before.min(),
                "before_median": before.median(), "before_max": before.max(),
                "after_min": after.min(), "after_median": after.median(), "after_max": after.max()})


def main() -> None:
    frame = pd.read_csv(RAW)
    changes = []
    # Fill the canonical 2001-2002 slots from their historical LBD aliases.
    mask_2001 = frame.cycle.eq("2001-2002")
    for source, target in ALIASES.items():
        source_values = frame.loc[mask_2001, source]
        frame.loc[mask_2001, target] = source_values
        changes.append({"rule_id": f"rename_{source}", "cycle": "2001-2002", "marker": target,
                        "equation": f"{target}={source}", "values_transformed": int(source_values.notna().sum())})

    rules = [
        ("1999-2000", "LBXSCR", "BIO-CAL-01", lambda x: 1.013*x + .147, "1.013*x + 0.147"),
        ("2005-2006", "LBXSCR", "BIO-CAL-02", lambda x: .978*x - .016, "0.978*x - 0.016"),
        ("2017-2018", "LBXSAL", "BIO-CAL-03", lambda x: 1.044*x + .01128, "1.044*x + 0.01128"),
        ("2017-2018", "LBXSAPSI", "BIO-CAL-04", lambda x: 10**(1.001*np.log10(x)-.04294), "10**(1.001*log10(x)-0.04294)"),
        ("2017-2018", "LBXSASSI", "BIO-CAL-05", lambda x: 1.018*x + 3.762, "1.018*x + 3.762"),
        ("2017-2018", "LBXSATSI", "BIO-CAL-06", lambda x: 1.013*x + 2.688, "1.013*x + 2.688"),
        ("2017-2018", "LBXSBU", "BIO-CAL-07", lambda x: 1.001*x + .4488, "1.001*x + 0.4488"),
        ("2017-2018", "LBXSCR", "BIO-CAL-08", lambda x: 1.051*x - .06945, "1.051*x - 0.06945"),
        ("2017-2018", "LBXSGTSI", "BIO-CAL-09", lambda x: .8042*x + 2.363, "0.8042*x + 2.363"),
        ("2017-2018", "LBXSIR", "BIO-CAL-10", lambda x: .9776*x - 4.494, "0.9776*x - 4.494"),
        ("2017-2018", "LBXSLDSI", "BIO-CAL-11", lambda x: .8568*x + 2.062, "0.8568*x + 2.062"),
        ("2017-2018", "LBXSUA", "BIO-CAL-12", lambda x: .9323*x + .2326, "0.9323*x + 0.2326"),
        ("2017-2018", "LBXSTR", "BIO-CAL-13", lambda x: .9655*x - 7.020, "0.9655*x - 7.020"),
        ("2017-2018", "LBXSCH", "BIO-CAL-14", lambda x: 1.046*x - 2.203, "1.046*x - 2.203"),
    ]
    for cycle, marker, rule_id, function, equation in rules:
        apply_rule(frame, cycle, marker, rule_id, function, equation, changes)

    nonpositive_rows = []
    for marker in MARKERS:
        bad = frame[marker].notna() & frame[marker].le(0)
        nonpositive_rows.append({"marker": marker, "nonpositive_to_missing": int(bad.sum())})
        frame.loc[bad, marker] = np.nan

    before_complete = len(frame)
    clean = frame[["SEQN", "cycle"] + MARKERS].dropna(subset=MARKERS).reset_index(drop=True)
    if clean.SEQN.duplicated().any():
        raise ValueError("harmonised biochemistry contains duplicate SEQN")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    clean.to_csv(OUTPUT / "standard_biochemistry.csv", index=False)
    pd.DataFrame(changes).to_csv(OUTPUT / "standard_biochemistry_change_log.csv", index=False)
    pd.DataFrame(nonpositive_rows).to_csv(OUTPUT / "standard_biochemistry_nonpositive_log.csv", index=False)
    by_cycle = clean.groupby("cycle").size().rename("complete_rows").reset_index()
    by_cycle.to_csv(OUTPUT / "standard_biochemistry_by_cycle.csv", index=False)
    miss = frame.groupby("cycle")[MARKERS].agg(lambda x: int(x.isna().sum())).reset_index()
    miss.to_csv(OUTPUT / "standard_biochemistry_missingness.csv", index=False)
    print(f"harmonised biochemistry: {len(clean):,} complete rows from {before_complete:,} imported participants")


if __name__ == "__main__":
    main()

