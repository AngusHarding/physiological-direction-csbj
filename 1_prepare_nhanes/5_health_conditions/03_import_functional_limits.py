"""Import the three raw NHANES functional-impairment reference screens."""

from hashlib import sha256
from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key  # noqa: E402

ITEMS = ["PFQ051", "PFQ054", "PFQ057"]


def fingerprint(path: Path) -> str:
    digest = sha256(path.read_bytes()).hexdigest().upper()
    return digest


def main() -> None:
    sources = pd.read_csv(ROOT / "config/functional_impairment_sources.csv")
    raw_root = raw_data_root()
    frames, manifest = [], []
    for row in sources.itertuples(index=False):
        path = raw_root / row.relative_path
        raw = read_nhanes_xpt(path)
        present = [item for item in ITEMS if item in raw.columns]
        frame = raw[["SEQN"] + present].copy()
        for item in ITEMS:
            if item not in frame:
                frame[item] = float("nan")  # public_code: numeric empty column, so pandas keeps the column type now and in future versions
        frame["cycle"] = row.cycle
        frames.append(frame[["SEQN", "cycle"] + ITEMS])
        manifest.append({"cycle": row.cycle, "relative_path": row.relative_path,
                         "sha256": fingerprint(path), "rows": len(raw),
                         "items_present": "|".join(present)})
    combined = pd.concat(frames, ignore_index=True)
    validate_participant_key(combined, "functional impairment")
    out = ROOT / "outputs/01_imported"
    out.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out / "functional_impairment.csv", index=False)
    pd.DataFrame(manifest).to_csv(ROOT / "outputs/provenance/functional_impairment_raw_manifest.csv", index=False)
    print(f"functional impairment: {len(combined):,} unique participants")


if __name__ == "__main__":
    main()
