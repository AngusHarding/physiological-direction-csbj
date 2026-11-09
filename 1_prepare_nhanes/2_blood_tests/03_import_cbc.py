"""Import analyser-reported CBC measurements without screening or derivation."""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

MEASURED = ["LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT",
            "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXMCVSI", "LBXRDW",
            "LBXPLTSI", "LBXMPSI"]


def main() -> None:
    sources = pd.read_csv(ROOT / "config" / "cbc_sources.csv")
    destination = ROOT / "outputs" / "01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    frames, summaries, missingness = [], [], []
    data_root = raw_data_root()
    for source in sources.itertuples(index=False):
        raw = read_nhanes_xpt(data_root / Path(source.relative_path))
        absent = [column for column in MEASURED if column not in raw.columns]
        if absent:
            raise KeyError(f"{source.cycle}: expected CBC variables absent: {absent}")
        frame = raw[["SEQN"] + MEASURED].copy()
        frame.insert(1, "cycle", source.cycle)
        validate_participant_key(frame, f"CBC {source.cycle}")
        frames.append(frame)
        summaries.append({"cycle": source.cycle, "source_file": source.relative_path,
                          "source_rows": len(raw), "output_rows": len(frame)})
        for column in MEASURED:
            missingness.append({"cycle": source.cycle, "variable": column,
                                "nonmissing": int(frame[column].notna().sum()),
                                "missing": int(frame[column].isna().sum())})
    combined = pd.concat(frames, ignore_index=True)
    validate_participant_key(combined, "stacked raw CBC")
    combined.to_csv(destination / "cbc_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "cbc_import_summary.csv", index=False)
    pd.DataFrame(missingness).to_csv(destination / "cbc_raw_missingness.csv", index=False)
    print(f"raw CBC: {len(combined):,} unique participants across {combined.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()
