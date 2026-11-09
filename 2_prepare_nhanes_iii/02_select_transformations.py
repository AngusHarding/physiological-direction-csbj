"""Select transformations within NHANES III, and compare against continuous NHANES.

SPECIFICATION: SPECIFICATION.md Decision 7 and section 6a layers 3-4.

Outcome-blind: no mortality variable is read.

The selection PROCEDURE is imported from the main pipeline rather than
reimplemented -- `select_full_cohort_transformations.evaluate_marker` -- so the
rule provably transfers unchanged. Only the DATA differ, and therefore only the
choices. Reimplementing it would risk a silent divergence in exactly the step
Decision 7 says must be identical.

Selection is unsupervised and runs on the NHANES III analysis cohort, not on its
young-healthy subset. That mirrors the main pipeline's correction of defect 1
(SUPPLEMENTARY section 0.4): a transformation must normalise the distribution it
is actually applied to.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from nhanes3_locations import CONFIG, MODERN_CONFIG, OUTPUTS, modern_transformation_rule  # public_code

MAIN = modern_transformation_rule()  # public_code: 1_prepare_nhanes/7_marker_panel/02_select_transformations.py

OUT = OUTPUTS  # public_code


def main() -> None:
    panel = pd.read_csv(OUT / "nhanes3_panel_raw.csv")
    markers = pd.read_csv(OUT / "nhanes3_panel_markers.csv").marker.tolist()
    print(f"cohort {len(panel):,} | markers {len(markers)}")
    print(f"rule imported from main pipeline: candidates={MAIN.CANDIDATES}, "
          f"skew<={MAIN.ABS_SKEW_LIMIT}, kurt<={MAIN.ABS_EXCESS_KURTOSIS_LIMIT}, "
          f"qq_r2>={MAIN.QQ_R2_MINIMUM}")

    ledgers, chosen = [], []
    for m in markers:
        values = panel[m].dropna().to_numpy(float)
        # Step 1 of the locked rule: keep identity if the raw scale already
        # satisfies all three adequacy criteria.
        identity = MAIN.diagnostics(values, "identity")
        retain = (identity["abs_skew"] <= MAIN.ABS_SKEW_LIMIT
                  and identity["abs_excess_kurtosis"] <= MAIN.ABS_EXCESS_KURTOSIS_LIMIT
                  and (1.0 - identity["qq_deviation"]) >= MAIN.QQ_R2_MINIMUM)
        ledger = MAIN.evaluate_marker(values, m, retain_identity=retain)
        ledgers.append(ledger)
        win = ledger.loc[ledger.selected].iloc[0]
        chosen.append({
            "marker": m, "selected_transform": win.candidate,
            "orientation": win.orientation, "basis": win.get("selection_basis", ""),
            "raw_skew": round(float(identity["abs_skew"]), 3),
            "transformed_skew": round(float(win.abs_skew), 3),
            "transformed_kurtosis": round(float(win.abs_excess_kurtosis), 3),
            "qq_r2": round(float(win.qq_r2), 4),
            "adequacy": win.adequacy,
        })
    pd.concat(ledgers).to_csv(OUT / "transformation_ledger.csv", index=False)
    sel = pd.DataFrame(chosen)
    sel.to_csv(OUT / "transformation_selected.csv", index=False)

    print("\n=== transformations selected in NHANES III ===")
    print(sel[["marker", "selected_transform", "raw_skew",
               "transformed_skew", "adequacy"]].to_string(index=False))
    print("\ncounts:", dict(sel.selected_transform.value_counts()))
    print("adequacy:", dict(sel.adequacy.value_counts()))

    # ---- Layer 4: side by side against continuous NHANES ------------------
    # Compared against the recipe the main analysis APPLIES (v2). The submitted
    # version read final_38_panel_spec.csv, a superseded recipe, and reported 13
    # of 34 in agreement; against the applied recipe it is 22 of 34. Corrected
    # 15 October 2026 (AD-009 amendment 1; found in revision note N3a, R013).
    # No estimate depends on this table.
    main_spec = pd.read_csv(
        MODERN_CONFIG / "final_38_transformations_v2.csv")  # public_code
    mapping = pd.read_csv(CONFIG / "nhanes3_marker_map.csv")  # public_code
    link = dict(zip(mapping.marker, mapping.main_panel_var))
    link.update({"ACR": "ACR", "eGFR": "eGFR", "WHtR": "WHtR"})
    main_tf = dict(zip(main_spec["marker"], main_spec["selected_transform"]))

    rows = []
    for r in sel.itertuples():
        mv = link.get(r.marker)
        rows.append({
            "marker": r.marker, "main_panel_var": mv,
            "nhanes3_transform": r.selected_transform,
            "main_transform": main_tf.get(mv, "not in main panel"),
            "agree": r.selected_transform == main_tf.get(mv),
        })
    cmp = pd.DataFrame(rows)
    cmp.to_csv(OUT / "transformation_comparison.csv", index=False)
    matched = cmp[cmp.main_transform.ne("not in main panel")]
    print("\n=== transformation agreement with continuous NHANES ===")
    print(f"comparable markers : {len(matched)}")
    print(f"same choice        : {int(matched.agree.sum())} "
          f"({100*matched.agree.mean():.0f}%)")
    diff = matched[~matched.agree]
    if len(diff):
        print("\ndiffering choices (a result of re-derivation, not a defect):")
        print(diff[["marker", "nhanes3_transform", "main_transform"]].to_string(index=False))
    print(f"\nwritten -> {OUT}")


if __name__ == "__main__":
    main()
