# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/3_body_and_bp/04_import_blood_pressure.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Reads all four blood-pressure readings and the pulse.
#
# METHODS       Data import for Supp. Section Zero, "Blood pressure and pulse".
# =============================================================================
"""Import the ten NHANES blood-pressure examination files without deriving values.

This first stage deliberately preserves readings 1 through 4. Reading 1 is not
used in the final average, but retaining it makes that choice visible and lets a
reviewer compare the excluded settling measurement with the retained repeats.
"""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key

COLUMNS = [
    "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4",
    "BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4",
    "BPXPLS", "BPXPULS",
]


def main() -> None:
    sources = pd.read_csv(ROOT / "config" / "cardiovascular_sources.csv")
    output_dir = ROOT / "outputs" / "01_imported"
    output_dir.mkdir(parents=True, exist_ok=True)
    frames, summaries, missingness = [], [], []
    data_root = raw_data_root()

    for source in sources.itertuples(index=False):
        raw = read_nhanes_xpt(data_root / Path(source.relative_path))
        # reindex fails visibly only after recording absent columns as all-missing;
        # the explicit check gives a clearer explanation to a novice reviewer.
        absent = [name for name in COLUMNS if name not in raw.columns]
        if absent:
            raise KeyError(f"{source.cycle}: expected BP variables absent: {absent}")
        frame = raw[["SEQN"] + COLUMNS].copy()
        frame.insert(1, "cycle", source.cycle)
        validate_participant_key(frame, f"cardiovascular {source.cycle}")
        frames.append(frame)
        summaries.append({"cycle": source.cycle, "source_file": source.relative_path,
                          "source_rows": len(raw), "output_rows": len(frame)})
        for variable in COLUMNS:
            missingness.append({"cycle": source.cycle, "variable": variable,
                                "nonmissing": int(frame[variable].notna().sum()),
                                "missing": int(frame[variable].isna().sum())})

    combined = pd.concat(frames, ignore_index=True)
    validate_participant_key(combined, "stacked raw cardiovascular data")
    combined.to_csv(output_dir / "cardiovascular_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(output_dir / "cardiovascular_import_summary.csv", index=False)
    pd.DataFrame(missingness).to_csv(output_dir / "cardiovascular_raw_missingness.csv", index=False)
    print(f"raw cardiovascular: {len(combined):,} unique participants across {combined.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()
