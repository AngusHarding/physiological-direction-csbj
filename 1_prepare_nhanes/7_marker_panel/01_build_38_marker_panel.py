"""Assemble and lock the untransformed final 38-biomarker panel.

Beginner's map
--------------
Every preceding chunk produced a participant table for one laboratory or
physiological component. This script starts with the eligible cohort and left-
joins those tables by ``SEQN``. A left join never invents participants: it keeps
the eligible rows and adds values where a component has a matching participant.

The historical 29,053-person cohort was locked before six auxiliary columns
were removed from scoring. We reproduce that exact completeness boundary, then
save the final 38 separately from screening/redundancy metadata.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
OUT = ROOT / "outputs/05_panel"

SOURCES = {
    "standard_biochemistry": ROOT / "outputs/02_harmonised/standard_biochemistry.csv",
    "cardiovascular": ROOT / "outputs/02_harmonised/cardiovascular.csv",
    "cbc_derived": ROOT / "outputs/03_derived/cbc_derived_features.csv",
    "hba1c": ROOT / "outputs/02_harmonised/glycated_haemoglobin.csv",
    "anthropometry_derived": ROOT / "outputs/03_derived/anthropometry_derived_features.csv",
    "cholesterol_hdl": ROOT / "outputs/02_harmonised/cholesterol_hdl.csv",
    "whole_blood_count": ROOT / "outputs/02_harmonised/whole_blood_count.csv",
    "folate": ROOT / "outputs/02_harmonised/folate.csv",
}


def egfr_ckdepi_2021(creatinine, age, sex):
    """CKD-EPI 2021 race-free creatinine eGFR; sex is NHANES 1=male, 2=female."""
    creatinine=np.asarray(creatinine,float); female=np.asarray(sex)==2
    kappa=np.where(female,.7,.9); alpha=np.where(female,-.241,-.302); ratio=creatinine/kappa
    estimate=142*np.minimum(ratio,1)**alpha*np.maximum(ratio,1)**-1.2*0.9938**np.asarray(age,float)
    return np.where(female,estimate*1.012,estimate)


def main() -> None:
    spec = pd.read_csv(ROOT / "config/panel_merge_spec.csv").sort_values("order")
    eligible = pd.read_csv(ROOT / "outputs/04_cohorts/eligible_adult_cohort.csv")
    panel = eligible[["SEQN", "RIAGENDR", "RIDAGEYR"]].copy()
    flow = [{"step": "eligible_base", "eligible_rows": len(panel),
             "cumulative_complete_rows": len(panel), "variables_added": 0}]
    cumulative = []

    # Read each source once and add only the columns assigned to it in the spec.
    for source, path in SOURCES.items():
        columns = spec.loc[spec.source.eq(source), "marker"].tolist()
        frame = pd.read_csv(path, usecols=["SEQN"] + columns)
        if frame.SEQN.duplicated().any():
            raise ValueError(f"{source} contains duplicate SEQN")
        panel = panel.merge(frame, on="SEQN", how="left", validate="one_to_one")
        cumulative.extend(columns)
        flow.append({"step": f"merge_{source}", "eligible_rows": len(panel),
                     "cumulative_complete_rows": int(panel[cumulative].notna().all(axis=1).sum()),
                     "variables_added": len(columns)})

    # ACR needs the measured urine inputs, which are not final markers themselves.
    urine = pd.read_csv(ROOT / "outputs/02_harmonised/urine_albumin_creatinine.csv")
    panel = panel.merge(urine[["SEQN", "URXUMA", "URXUCR"]], on="SEQN", how="left", validate="one_to_one")
    panel["ACR"] = panel.URXUMA / panel.URXUCR * 100
    panel["eGFR"] = egfr_ckdepi_2021(panel.LBXSCR, panel.RIDAGEYR, panel.RIAGENDR)
    panel = panel.replace([np.inf, -np.inf], np.nan).drop(columns=["URXUMA", "URXUCR"])
    cumulative.extend(["ACR", "eGFR"])
    flow.append({"step": "derive_ACR_eGFR", "eligible_rows": len(panel),
                 "cumulative_complete_rows": int(panel[cumulative].notna().all(axis=1).sum()),
                 "variables_added": 2})

    completeness = spec.marker.tolist()
    missing = pd.DataFrame({"marker": completeness,
                            "missing_eligible_rows": [int(panel[m].isna().sum()) for m in completeness],
                            "nonmissing_eligible_rows": [int(panel[m].notna().sum()) for m in completeness]})
    # Sequential attrition makes the ordering explicit; it is an audit device,
    # not a claim that a participant has only one reason for missingness.
    remaining = pd.Series(True, index=panel.index); attrition=[]
    for marker in completeness:
        newly_excluded = remaining & panel[marker].isna()
        attrition.append({"marker": marker, "role": spec.set_index("marker").loc[marker,"role"],
                          "newly_excluded": int(newly_excluded.sum()),
                          "remaining_after_step": int((remaining & ~panel[marker].isna()).sum())})
        remaining &= panel[marker].notna()
    locked = panel.loc[remaining].reset_index(drop=True)

    # The panel specification, not merge order, owns the public column order.
    final_spec = pd.read_csv(ROOT / "config/final_38_panel_spec.csv")
    final_markers = final_spec.marker.tolist()
    auxiliary = spec.loc[~spec.role.eq("final_38"), "marker"].tolist()
    final = locked[["SEQN", "RIAGENDR", "RIDAGEYR"] + final_markers]
    metadata = locked[["SEQN"] + auxiliary]
    if len(final_markers) != 38 or final.isna().any().any():
        raise AssertionError("final panel is not 38-marker complete-case")

    OUT.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUT / "final_38_biomarker_panel_raw.csv", index=False)
    metadata.to_csv(OUT / "final_38_panel_screening_metadata.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "final_38_panel_merge_flow.csv", index=False)
    missing.to_csv(OUT / "final_38_panel_missingness.csv", index=False)
    pd.DataFrame(attrition).to_csv(OUT / "final_38_panel_attrition_by_marker.csv", index=False)
    final_spec.merge(spec[["marker", "source"]], on="marker", how="left", validate="one_to_one").merge(
        pd.read_csv(ROOT / "config/biomarker_registry.csv"), on="marker", how="left",
        validate="one_to_one").to_csv(OUT / "final_38_panel_source_map.csv", index=False)
    print(f"locked raw final-38 panel: {len(final):,} participants x {len(final_markers)} biomarkers")
    print(f"sex counts: M={(final.RIAGENDR==1).sum():,}; F={(final.RIAGENDR==2).sum():,}")


if __name__ == "__main__":
    main()
