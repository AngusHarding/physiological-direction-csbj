"""Import urine kidney measurements and the separate kidney-failure flag."""

from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key


def main() -> None:
    sources = pd.read_csv(ROOT / "config" / "kidney_sources.csv")
    base = raw_data_root()
    destination = ROOT / "outputs" / "01_imported"
    destination.mkdir(parents=True, exist_ok=True)
    urine_frames, flag_frames, summaries = [], [], []
    for source in sources.itertuples(index=False):
        raw = read_nhanes_xpt(base / Path(source.relative_path))
        if source.component == "urine":
            frame = raw[["SEQN", "URXUMA", "URXUCR"]].copy()
            frame.insert(1, "cycle", source.cycle)
            urine_frames.append(frame)
        else:
            # Only definite yes/no answers survive. Refused, don't-know and SAS
            # missing values remain unknown rather than being labelled healthy.
            frame = raw[["SEQN"]].copy()
            frame["cycle"] = source.cycle
            frame["kidney_fail"] = raw[source.source_variable].map({1.0: 1.0, 2.0: 0.0})
            flag_frames.append(frame[["SEQN", "cycle", "kidney_fail"]])
        validate_participant_key(frame, f"{source.component} {source.cycle}")
        summaries.append({"component": source.component, "cycle": source.cycle,
                          "source_file": source.relative_path, "source_rows": len(raw)})
    urine = pd.concat(urine_frames, ignore_index=True)
    flags_raw = pd.concat(flag_frames, ignore_index=True)
    validate_participant_key(urine, "stacked urine")
    validate_participant_key(flags_raw, "stacked raw kidney flags")
    urine.to_csv(destination / "kidney_urine_raw.csv", index=False)
    flags_raw.to_csv(destination / "kidney_failure_flag_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(destination / "kidney_import_summary.csv", index=False)
    print(f"kidney import: {len(urine):,} urine records; {len(flags_raw):,} questionnaire records")


if __name__ == "__main__":
    main()
