"""Import the two smoking-status questions from all ten NHANES cycles.

This stage preserves questionnaire codes. It does not classify participants or
drop unanswered records; those scientific decisions belong to harmonisation.
"""

from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key


def main() -> None:
    sources = pd.read_csv(ROOT / "config/smoking_sources.csv")
    destination = ROOT / "outputs/01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    frames, summaries = [], []

    for row in sources.itertuples(index=False):
        raw = read_nhanes_xpt(raw_data_root() / row.source_file)
        required = ["SEQN", "SMQ020", "SMQ040"]
        absent = [column for column in required if column not in raw]
        if absent:
            raise KeyError(f"{row.cycle}: SMQ lacks {absent}")
        block = raw[required].copy()
        block.insert(1, "cycle", row.cycle)
        validate_participant_key(block, f"SMQ {row.cycle}")
        frames.append(block)
        summaries.append({
            "cycle": row.cycle,
            "source_file": row.source_file,
            "source_rows": len(raw),
            "output_rows": len(block),
            "SMQ020_answered_n": int(block.SMQ020.isin([1, 2]).sum()),
            "SMQ040_answered_n": int(block.SMQ040.isin([1, 2, 3]).sum()),
        })

    output = pd.concat(frames, ignore_index=True)
    validate_participant_key(output, "stacked SMQ")
    output.to_csv(destination / "smoking_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "smoking_import_summary.csv", index=False)
    print(f"Smoking import: {len(output):,} unique participants across {output.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()

