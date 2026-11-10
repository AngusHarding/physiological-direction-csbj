# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/4_lifestyle/05_import_physical_activity.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Reads the physical-activity questions, keeping which
#               questionnaire version each cycle used.
#
# METHODS       Data import for Table 1 and Table 4.
# =============================================================================
"""Import physical-activity source responses while preserving instrument era."""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

GATES = ["PAQ605", "PAQ620", "PAQ635", "PAQ650", "PAQ665"]
DOSE_COLUMNS = ["PAQ610", "PAD615", "PAQ625", "PAD630", "PAQ640", "PAD645",
                "PAQ655", "PAD660", "PAQ670", "PAD675"]
ALL = ["SEQN", "cycle", "instrument_era", "PAQ180"] + GATES + DOSE_COLUMNS


def main() -> None:
    sources = pd.read_csv(ROOT / "config/physical_activity_sources.csv")
    out = ROOT / "outputs/01_imported"
    out.mkdir(parents=True, exist_ok=True)
    frames, summaries = [], []
    for row in sources.itertuples(index=False):
        raw = read_nhanes_xpt(raw_data_root() / row.source_file)
        expected = ["SEQN", "PAQ180"] if row.instrument_era == "PAQ180" else ["SEQN"] + GATES + DOSE_COLUMNS
        absent = [c for c in expected if c not in raw]
        if absent:
            raise KeyError(f"{row.cycle}: PAQ lacks {absent}")
        block = raw[expected].copy()
        block.insert(1, "cycle", row.cycle)
        block.insert(2, "instrument_era", row.instrument_era)
        block = block.reindex(columns=ALL)
        validate_participant_key(block, f"PAQ {row.cycle}")
        frames.append(block)
        summaries.append({"cycle": row.cycle, "instrument_era": row.instrument_era,
                          "source_file": row.source_file, "source_rows": len(raw), "output_rows": len(block)})
    combined = pd.concat(frames, ignore_index=True)
    validate_participant_key(combined, "stacked PAQ")
    combined.to_csv(out / "physical_activity_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "physical_activity_import_summary.csv", index=False)
    print(f"Physical-activity import: {len(combined):,} participants across {combined.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()

