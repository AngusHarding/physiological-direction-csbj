"""Import full-sample cholesterol/HDL and fasting triglycerides separately."""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key


def main() -> None:
    sources = pd.read_csv(ROOT / "config/lipid_sources.csv")
    base = raw_data_root()
    destination = ROOT / "outputs/01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    cholesterol, triglycerides, summaries = [], [], []
    for source in sources.itertuples(index=False):
        tc = read_nhanes_xpt(base / Path(source.total_cholesterol_path))[["SEQN", "LBXTC"]]
        hdl = read_nhanes_xpt(base / Path(source.hdl_path))[["SEQN", source.hdl_source_variable]]
        validate_participant_key(tc, f"total cholesterol {source.cycle}")
        validate_participant_key(hdl, f"HDL {source.cycle}")
        merged = tc.merge(hdl, on="SEQN", how="inner", validate="one_to_one")
        merged.insert(1, "cycle", source.cycle)
        # Preserve each historical HDL name in its own slot. Harmonisation into
        # one canonical field happens in the next, explicitly reviewable stage.
        merged = merged.rename(columns={source.hdl_source_variable: "HDL_source_value"})
        merged.insert(3, "HDL_source_variable", source.hdl_source_variable)
        cholesterol.append(merged)

        trig_raw = read_nhanes_xpt(base / Path(source.fasting_triglyceride_path))
        trig = trig_raw[["SEQN", "LBXTR"]].copy()
        trig.insert(1, "cycle", source.cycle)
        validate_participant_key(trig, f"fasting triglyceride {source.cycle}")
        triglycerides.append(trig)
        summaries.append({"cycle": source.cycle, "tc_rows": len(tc), "hdl_rows": len(hdl),
                          "paired_rows": len(merged), "fasting_triglyceride_rows": len(trig)})
    tc_hdl = pd.concat(cholesterol, ignore_index=True)
    trig = pd.concat(triglycerides, ignore_index=True)
    validate_participant_key(tc_hdl, "stacked cholesterol/HDL")
    validate_participant_key(trig, "stacked fasting triglycerides")
    tc_hdl.to_csv(destination / "cholesterol_hdl_raw.csv", index=False)
    trig.to_csv(destination / "fasting_triglycerides_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "lipid_import_summary.csv", index=False)
    print(f"lipid import: {len(tc_hdl):,} TC/HDL records; {len(trig):,} fasting records")


if __name__ == "__main__":
    main()
