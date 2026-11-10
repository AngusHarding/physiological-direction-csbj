# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/1_cohort/04_import_survey_design.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Reads the survey-design variables (strata, sampling units,
#               weights) from the ten DEMO files.
#
# METHODS       Used by the survey-design bootstraps ("Uncertainty":
#               "primary sampling units within strata"). Point estimates are
#               unweighted (Methods, "Experimental design").
# =============================================================================
"""Extract NHANES complex-survey design variables from all ten DEMO files.

The public analysis spans ten equal two-year cycles. WTMEC2YR is therefore
divided by ten to create a pooled 20-year examination weight. The original
weight is retained so reviewers can verify the calculation.
"""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

FIELDS = ["SEQN", "WTMEC2YR", "WTINT2YR", "SDMVSTRA", "SDMVPSU"]


def main() -> None:
    sources = pd.read_csv(ROOT / "config/chunk1_sources.csv")
    frames, audit = [], []
    for row in sources.itertuples(index=False):
        raw = read_nhanes_xpt(raw_data_root() / row.demo_path)
        absent = [field for field in FIELDS if field not in raw]
        if absent:
            raise KeyError(f"{row.cycle}: missing survey fields {absent}")
        block = raw[FIELDS].copy()
        block.insert(1, "cycle", row.cycle)
        block["WTMEC20YR"] = block.WTMEC2YR / len(sources)
        block["WTINT20YR"] = block.WTINT2YR / len(sources)
        # PSU codes repeat across strata; the pair is the cluster identifier.
        block["survey_cluster"] = block.SDMVSTRA.astype("Int64").astype(str) + "_" + block.SDMVPSU.astype("Int64").astype(str)
        validate_participant_key(block, f"survey design {row.cycle}")
        frames.append(block)
        audit.append({"cycle": row.cycle, "rows": len(block),
                      "positive_mec_weight_n": int(block.WTMEC2YR.gt(0).sum()),
                      "strata_n": block.SDMVSTRA.nunique(),
                      "clusters_n": block.survey_cluster.nunique(),
                      "sum_WTMEC2YR": block.WTMEC2YR.sum(),
                      "sum_WTMEC20YR": block.WTMEC20YR.sum()})
    result = pd.concat(frames, ignore_index=True)
    validate_participant_key(result, "pooled survey design")
    out = ROOT / "outputs/01_imported"; out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out / "survey_design_1999_2018.csv", index=False)
    pd.DataFrame(audit).to_csv(out / "survey_design_import_audit.csv", index=False)
    print(f"Survey design: {len(result):,} participants; {result.SDMVSTRA.nunique()} strata; {result.survey_cluster.nunique()} clusters")


if __name__ == "__main__":
    main()
