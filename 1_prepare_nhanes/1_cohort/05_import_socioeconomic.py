"""Import minimally interpreted socioeconomic DEMO variables across ten cycles.

Run this file directly.  It reads the same DEMO XPT file per survey cycle that
`import_demographics.py` reads, keeps the three socioeconomic variables, adds
the cycle label, stacks the ten tables, and saves one CSV.  It interprets
nothing: NHANES codes are preserved exactly as published.

WHY THIS IS A SEPARATE FILE
  `demographics.csv` is consumed downstream and its content determines hashes
  that later modules verify.  Adding columns to it would change those hashes and
  invalidate work that does not depend on socioeconomic position at all.
  Writing a separate table adds the variables without disturbing anything.

VARIABLES
  RIDRETH1  race/ethnicity, as published (1-5)
  DMDEDUC2  education, adults 20+ (1-5; 7 refused, 9 don't know)
  INDFMPIR  family income-to-poverty ratio (0-5, capped at 5 by NHANES)

  Codes 7 and 9 are sentinel non-responses, not education levels.  They are NOT
  recoded here, because import records what NHANES published; any analysis that
  uses DMDEDUC2 must handle them explicitly.  A QC column counts them so the
  decision cannot be made by accident downstream.
"""

from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))

from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

SENTINEL_EDUCATION = {7.0, 9.0}


def main() -> None:
    sources = pd.read_csv(ROOT / "config" / "chunk1_sources.csv")
    destination = ROOT / "outputs" / "01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    frames = []
    summaries = []
    data_root = raw_data_root()

    for row in sources.itertuples(index=False):
        path = data_root / Path(row.demo_path)
        raw = read_nhanes_xpt(path)
        required = ["SEQN", "RIDRETH1", "DMDEDUC2", "INDFMPIR"]
        absent = [column for column in required if column not in raw.columns]
        if absent:
            raise KeyError(f"{row.cycle}: DEMO lacks {absent}")
        clean = raw[required].copy()
        clean.insert(1, "cycle", row.cycle)
        validate_participant_key(clean, f"DEMO socioeconomic {row.cycle}")
        frames.append(clean)
        summaries.append(
            {
                "cycle": row.cycle,
                "source_file": row.demo_path,
                "source_rows": len(raw),
                "output_rows": len(clean),
                "ethnicity_answered_n": int(clean.RIDRETH1.notna().sum()),
                "education_answered_n": int(
                    clean.DMDEDUC2.notna().sum()
                    - clean.DMDEDUC2.isin(SENTINEL_EDUCATION).sum()
                ),
                "education_refused_or_unknown_n": int(
                    clean.DMDEDUC2.isin(SENTINEL_EDUCATION).sum()
                ),
                "income_ratio_answered_n": int(clean.INDFMPIR.notna().sum()),
            }
        )

    output = pd.concat(frames, ignore_index=True)
    validate_participant_key(output, "stacked DEMO socioeconomic")
    output.to_csv(destination / "socioeconomic.csv", index=False)
    pd.DataFrame(summaries).to_csv(
        destination / "socioeconomic_import_summary.csv", index=False
    )
    print(
        f"socioeconomic: {len(output):,} unique participants across "
        f"{output.cycle.nunique()} cycles"
    )
    print(
        f"  ethnicity answered   : {output.RIDRETH1.notna().sum():,}\n"
        f"  education answered   : {int(output.DMDEDUC2.notna().sum() - output.DMDEDUC2.isin(SENTINEL_EDUCATION).sum()):,}"
        f"  (refused/unknown {int(output.DMDEDUC2.isin(SENTINEL_EDUCATION).sum()):,})\n"
        f"  income ratio answered: {output.INDFMPIR.notna().sum():,}"
    )


if __name__ == "__main__":
    main()
