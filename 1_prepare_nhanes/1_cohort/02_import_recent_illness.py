"""Import source-coded recent-illness questionnaire variables across ten cycles.

This mirrors import_demographics.py but reads the Health Status Questionnaire
(HSQ).  At this stage we preserve the original NHANES codes.  In particular,
missing/refused/unknown responses are not converted to 'No'.
"""

from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent       # src/01_import
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))

from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key


def main() -> None:
    # Read the ten-cycle file recipe and create the output folder if necessary.
    sources = pd.read_csv(ROOT / "config" / "chunk1_sources.csv")
    destination = ROOT / "outputs" / "01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    frames = []
    summaries = []
    data_root = raw_data_root()

    for row in sources.itertuples(index=False):
        path = data_root / Path(row.hsq_path)
        raw = read_nhanes_xpt(path)
        # HSQ500/510/520 are the three recent-illness questions used by the
        # historical eligible-cohort definition.
        required = ["SEQN", "HSQ500", "HSQ510", "HSQ520"]
        absent = [column for column in required if column not in raw.columns]
        if absent:
            raise KeyError(f"{row.cycle}: HSQ lacks {absent}")
        clean = raw[required].copy()
        clean.insert(1, "cycle", row.cycle)
        validate_participant_key(clean, f"HSQ {row.cycle}")
        frames.append(clean)
        summaries.append(
            {
                "cycle": row.cycle,
                "source_file": row.hsq_path,
                "source_rows": len(raw),
                "output_rows": len(clean),
                # notna() tests each cell. all(axis=1) asks whether all three
                # question cells are present for that participant/row.
                "complete_required_responses": int(clean[["HSQ500", "HSQ510", "HSQ520"]].notna().all(axis=1).sum()),
                "missing_seqn": int(clean.SEQN.isna().sum()),
                "duplicate_seqn_rows": int(clean.duplicated("SEQN", keep=False).sum()),
            }
        )

    # Stack ten two-year cycles into one participant table.
    output = pd.concat(frames, ignore_index=True)
    validate_participant_key(output, "stacked HSQ")
    output.to_csv(destination / "recent_illness.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "recent_illness_import_summary.csv", index=False)
    print(f"recent illness: {len(output):,} unique participants across {output.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()
