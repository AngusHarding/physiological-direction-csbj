"""Import chronic-disease questionnaire fields without interpreting responses.

Each questionnaire remains a separate participant table until it is merged by
SEQN.  This matters because an absent osteoporosis questionnaire is structural
non-coverage, not evidence that every participant answered "no".
"""

from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key


DOMAIN_FIELDS = {
    "mcq": ["MCQ010", "MCQ092", "MCQ160A", "MCQ160B", "MCQ160C",
            "MCQ160D", "MCQ160E", "MCQ160F", "MCQ160G", "MCQ160I",
            "MCD160M", "MCQ160M", "MCQ160K", "MCQ160L", "MCQ160N",
            "MCQ160O", "MCQ220"],
    "diq": ["DIQ010"],
    "kiq": ["KIQ020", "KIQ022"],
    "bpq": ["BPQ020"],
    "osq": ["OSQ060"],
}


def read_domain(cycle: str, domain: str, relative_path: str) -> pd.DataFrame:
    """Read available requested fields; fail if an expected file has no fields."""
    raw = read_nhanes_xpt(raw_data_root() / relative_path)
    fields = [field for field in DOMAIN_FIELDS[domain] if field in raw.columns]
    if not fields:
        raise KeyError(f"{cycle} {domain.upper()}: none of {DOMAIN_FIELDS[domain]} found")
    block = raw[["SEQN"] + fields].copy()
    validate_participant_key(block, f"{cycle} {domain.upper()}")
    return block


def main() -> None:
    sources = pd.read_csv(ROOT / "config/disease_sources.csv", keep_default_na=False)
    out = ROOT / "outputs/01_imported"
    out.mkdir(parents=True, exist_ok=True)
    frames, summaries = [], []

    for row in sources.itertuples(index=False):
        merged = None
        for domain in DOMAIN_FIELDS:
            relative_path = getattr(row, f"{domain}_file")
            if not relative_path:
                continue
            block = read_domain(row.cycle, domain, relative_path)
            merged = block if merged is None else merged.merge(
                block, on="SEQN", how="outer", validate="one_to_one"
            )
            summaries.append({
                "cycle": row.cycle, "domain": domain.upper(),
                "source_file": relative_path, "source_rows": len(block),
                "fields_imported": "|".join(c for c in block if c != "SEQN"),
            })
        merged.insert(1, "cycle", row.cycle)
        validate_participant_key(merged, f"{row.cycle} combined disease import")
        frames.append(merged)

    combined = pd.concat(frames, ignore_index=True, sort=False)
    validate_participant_key(combined, "stacked disease import")
    combined.to_csv(out / "self_reported_disease_raw.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "self_reported_disease_import_summary.csv", index=False)
    print(f"Disease import: {len(combined):,} participants across {combined.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()
