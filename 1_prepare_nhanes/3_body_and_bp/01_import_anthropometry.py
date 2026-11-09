"""Import measured weight, height and waist circumference without derivation."""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

MEASURED = ["BMXWT", "BMXHT", "BMXWAIST"]


def main() -> None:
    sources = pd.read_csv(ROOT / "config/anthropometry_sources.csv")
    destination = ROOT / "outputs/01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    frames, summaries, missingness = [], [], []
    base = raw_data_root()
    for source in sources.itertuples(index=False):
        raw = read_nhanes_xpt(base / Path(source.relative_path))
        frame = raw[["SEQN"] + MEASURED].copy()
        frame.insert(1, "cycle", source.cycle)
        validate_participant_key(frame, f"anthropometry {source.cycle}")
        frames.append(frame)
        summaries.append({"cycle": source.cycle, "source_file": source.relative_path,
                          "source_rows": len(raw), "output_rows": len(frame)})
        for variable in MEASURED:
            missingness.append({"cycle": source.cycle, "variable": variable,
                                "nonmissing": int(frame[variable].notna().sum()),
                                "missing": int(frame[variable].isna().sum())})
    combined = pd.concat(frames, ignore_index=True)
    validate_participant_key(combined, "stacked raw anthropometry")
    combined.to_csv(destination / "anthropometry_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "anthropometry_import_summary.csv", index=False)
    pd.DataFrame(missingness).to_csv(destination / "anthropometry_raw_missingness.csv", index=False)
    print(f"raw anthropometry: {len(combined):,} examination participants")


if __name__ == "__main__":
    main()
