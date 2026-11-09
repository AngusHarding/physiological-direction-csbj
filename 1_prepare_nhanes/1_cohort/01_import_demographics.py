"""Import minimally interpreted DEMO variables across ten NHANES cycles.

Run this file directly.  It reads one DEMO XPT file per survey cycle, keeps the
four variables needed for Chunk 1, adds the cycle label, stacks the ten tables,
and saves one CSV.  It does not yet decide who is eligible for the study.
"""

from pathlib import Path
import sys

import pandas as pd

# HERE is this script's folder (src/01_import).
# ROOT is the data_preparation folder, two levels above the script's folder.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
# Python normally searches only standard/package locations.  Adding HERE lets
# this standalone script import our neighbouring nhanes_reader.py module.
sys.path.insert(0, str(HERE))

from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key


def main() -> None:
    # The source table is a small recipe: one DEMO and one HSQ path per cycle.
    sources = pd.read_csv(ROOT / "config" / "chunk1_sources.csv")
    destination = ROOT / "outputs" / "01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    # We collect each cycle table in frames, then stack them after the loop.
    # summaries contains one compact QC row per cycle.
    frames = []
    summaries = []
    data_root = raw_data_root()

    for row in sources.itertuples(index=False):
        # itertuples gives convenient attributes such as row.cycle and
        # row.demo_path for each row of the source recipe.
        path = data_root / Path(row.demo_path)
        raw = read_nhanes_xpt(path)
        # Keep source names at import so the relationship with NHANES is obvious:
        # age, sex and pregnancy status plus the participant identifier.
        required = ["SEQN", "RIDAGEYR", "RIAGENDR", "RIDEXPRG"]
        absent = [column for column in required if column not in raw.columns]
        if absent:
            raise KeyError(f"{row.cycle}: DEMO lacks {absent}")
        # .copy() creates an independent table containing only these columns.
        clean = raw[required].copy()
        # Insert the cycle beside SEQN so every row retains its survey origin.
        clean.insert(1, "cycle", row.cycle)
        validate_participant_key(clean, f"DEMO {row.cycle}")
        frames.append(clean)
        summaries.append(
            {
                "cycle": row.cycle,
                "source_file": row.demo_path,
                "source_rows": len(raw),
                "output_rows": len(clean),
                "missing_seqn": int(clean.SEQN.isna().sum()),
                "duplicate_seqn_rows": int(clean.duplicated("SEQN", keep=False).sum()),
            }
        )

    # concat places cycle tables one below another; ignore_index creates a fresh
    # 0..N-1 row index that is unrelated to the participant identifier.
    output = pd.concat(frames, ignore_index=True)
    validate_participant_key(output, "stacked DEMO")
    # index=False avoids writing pandas' internal row numbers as a fake column.
    output.to_csv(destination / "demographics.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "demographics_import_summary.csv", index=False)
    print(f"demographics: {len(output):,} unique participants across {output.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()
