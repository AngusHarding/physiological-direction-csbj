"""Import alcohol questionnaire variables without deriving consumption."""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

OLD_COLUMNS = ["SEQN", "ALQ101", "ALQ110", "ALQ120Q", "ALQ120U", "ALQ130"]
NEW_COLUMNS = ["SEQN", "ALQ111", "ALQ121", "ALQ130"]
ALL_COLUMNS = ["SEQN", "cycle", "instrument_era", "ALQ101", "ALQ110", "ALQ120Q",
               "ALQ120U", "ALQ111", "ALQ121", "ALQ130"]


def main() -> None:
    sources = pd.read_csv(ROOT / "config/alcohol_sources.csv")
    out = ROOT / "outputs/01_imported"
    out.mkdir(parents=True, exist_ok=True)
    frames, summaries = [], []
    for row in sources.itertuples(index=False):
        raw = read_nhanes_xpt(raw_data_root() / row.source_file)
        expected = OLD_COLUMNS if row.instrument_era == "count_unit" else NEW_COLUMNS
        # ALQ101 is absent in the first two cycles by questionnaire design.
        required = [c for c in expected if not (c == "ALQ101" and c not in raw)]
        absent = [c for c in required if c not in raw]
        if absent:
            raise KeyError(f"{row.cycle}: ALQ lacks {absent}")
        block = raw.reindex(columns=expected).copy()
        block.insert(1, "cycle", row.cycle)
        block.insert(2, "instrument_era", row.instrument_era)
        block = block.reindex(columns=ALL_COLUMNS)
        validate_participant_key(block, f"ALQ {row.cycle}")
        frames.append(block)
        summaries.append({"cycle": row.cycle, "instrument_era": row.instrument_era,
                          "source_file": row.source_file, "source_rows": len(raw),
                          "output_rows": len(block)})
    combined = pd.concat(frames, ignore_index=True)
    validate_participant_key(combined, "stacked ALQ")
    combined.to_csv(out / "alcohol_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "alcohol_import_summary.csv", index=False)
    print(f"Alcohol import: {len(combined):,} participants across {combined.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()

