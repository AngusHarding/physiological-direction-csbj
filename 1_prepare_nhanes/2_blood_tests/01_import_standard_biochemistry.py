# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/2_blood_tests/01_import_standard_biochemistry.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Reads the standard biochemistry files of the ten cycles,
#               restoring SAS stored zeros; no correction yet.
#
# METHODS       Supp. Section Zero, "Missing-value codes".
# =============================================================================
"""Import standard-biochemistry source columns without assay correction.

The output is deliberately 'raw-ish': it has restored SAS literal zeros and an
explicit cycle label, but keeps historical variable aliases. Harmonisation is a
separate script so a reviewer can inspect values before and after correction.
"""

from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

CANONICAL = [
    "LBXSAL", "LBXSATSI", "LBXSASSI", "LBXSAPSI", "LBXSBU", "LBXSCA", "LBXSC3SI",
    "LBXSCR", "LBXSGTSI", "LBXSIR", "LBXSLDSI", "LBXSPH", "LBXSTB", "LBXSTP",
    "LBXSUA", "LBXSNASI", "LBXSKSI", "LBXSCLSI", "LBXSOSSI", "LBXSGB", "LBXSCH",
    "LBXSTR", "LBXSGL",
]
ALIASES = ["LBDSAPSI", "LBDSLDSI", "LBDSCR", "LBDSPH", "LBDSTB"]


def main() -> None:
    sources = pd.read_csv(ROOT / "config" / "standard_biochemistry_sources.csv")
    destination = ROOT / "outputs" / "01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    frames, summaries, missingness = [], [], []
    data_root = raw_data_root()

    for row in sources.itertuples(index=False):
        raw = read_nhanes_xpt(data_root / Path(row.relative_path))
        available = [name for name in CANONICAL + ALIASES if name in raw.columns]
        # reindex creates absent alias/canonical columns as NaN, giving every
        # cycle the same wide schema without pretending names were harmonised.
        frame = raw[["SEQN"] + available].reindex(columns=["SEQN"] + CANONICAL + ALIASES)
        frame.insert(1, "cycle", row.cycle)
        validate_participant_key(frame, f"biochemistry {row.cycle}")
        frames.append(frame)
        summaries.append({"cycle": row.cycle, "source_file": row.relative_path,
                          "source_rows": len(raw), "output_rows": len(frame),
                          "source_variables_found": len(available)})
        for variable in available:
            missingness.append({"cycle": row.cycle, "source_variable": variable,
                                "nonmissing": int(raw[variable].notna().sum()),
                                "missing": int(raw[variable].isna().sum())})

    output = pd.concat(frames, ignore_index=True)
    validate_participant_key(output, "stacked raw biochemistry")
    output.to_csv(destination / "standard_biochemistry_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "standard_biochemistry_import_summary.csv", index=False)
    pd.DataFrame(missingness).to_csv(destination / "standard_biochemistry_raw_missingness.csv", index=False)
    print(f"raw biochemistry: {len(output):,} unique participants across {output.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()

