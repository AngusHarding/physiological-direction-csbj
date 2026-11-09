"""Extract immutable FPED-derived tables from the frozen CRAN data package.

This is a package-replication boundary, not yet an official-USDA reconstruction.
The source tarball is retained and hashed; no live R code or package installation
is used. ``pyreadr`` reads the archived R data objects into pandas.
"""

from pathlib import Path
import hashlib
import tarfile
import pandas as pd
import pyreadr

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT, HEI_PACKAGES  # public_code: folders and write guard (shared/locations.py)
ARCHIVE_DIR = HEI_PACKAGES  # public_code: the CRAN archives, read only
DATA_ARCHIVE = ARCHIVE_DIR / "heiscore.data_0.0.1.tar.gz"
CODE_ARCHIVE = ARCHIVE_DIR / "heiscore_0.1.4.tar.gz"
EXPECTED = {
    "heiscore.data_0.0.1.tar.gz": "a325f6032b07bc6a312250eb7c7c537bf0a453d1ebf06bf66d77144e099712de",
    "heiscore_0.1.4.tar.gz": "c21a7f72bcd0fbad8bbae4d2ecc55f0f3bf8cca270c6d0cc0f90946220f99f60",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_extracted() -> Path:
    for archive in (DATA_ARCHIVE, CODE_ARCHIVE):
        observed = sha256(archive)
        if observed != EXPECTED[archive.name]:
            raise ValueError(f"Hash mismatch for {archive.name}: {observed}")
    destination = ROOT / "outputs/hei_packages_extracted"  # public_code: unpack into outputs
    destination.mkdir(parents=True, exist_ok=True)
    required = destination / "heiscore.data/data/HEI_scoring_standards.rda"
    if not required.exists():
        with tarfile.open(DATA_ARCHIVE, "r:gz") as archive:
            archive.extractall(destination)
        with tarfile.open(CODE_ARCHIVE, "r:gz") as archive:
            archive.extractall(destination)
    return destination


def read_single_object(path: Path) -> pd.DataFrame:
    objects = pyreadr.read_r(str(path))
    if len(objects) != 1:
        raise ValueError(f"Expected one R object in {path}, found {list(objects)}")
    return next(iter(objects.values()))


def main() -> None:
    extracted = ensure_extracted()
    sources = pd.read_csv(ROOT / "config/hei_sources.csv")
    out = ROOT / "outputs/01_imported"
    out.mkdir(parents=True, exist_ok=True)
    blocks, summaries = [], []
    data_dir = extracted / "heiscore.data/data"
    for row in sources.itertuples(index=False):
        block = read_single_object(data_dir / f"{row.package_object}.rda")
        block.insert(1, "cycle", row.cycle)
        block.insert(2, "cycle_code", row.cycle_code)
        if block.SEQN.isna().any() or block.SEQN.duplicated().any():
            raise ValueError(f"{row.cycle}: invalid SEQN")
        blocks.append(block)
        summaries.append({"cycle": row.cycle, "package_object": row.package_object,
                          "rows": len(block), "columns": len(block.columns),
                          "two_day_weight_observed_n": int(block.WTDR2D.notna().sum())})
    combined = pd.concat(blocks, ignore_index=True)
    if combined.SEQN.duplicated().any():
        raise ValueError("Stacked HEI package data contain duplicate SEQN")
    combined.to_csv(out / "hei_package_fped_inputs.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "hei_package_import_summary.csv", index=False)

    standards = read_single_object(data_dir / "HEI_scoring_standards.rda")
    # public_code: config/ is read only, so check the stored copy instead of rewriting it
    stored = (ROOT / "config/hei_2020_scoring_standards.csv").read_text(encoding="utf-8")
    if standards.to_csv(index=False, lineterminator="\n") != stored.replace("\r\n", "\n"):
        raise ValueError("config/hei_2020_scoring_standards.csv differs from the CRAN package's "
                         "HEI_scoring_standards; check which one is right before scoring HEI")
    pd.DataFrame([{"file": p.name, "sha256": sha256(p), "bytes": p.stat().st_size}
                  for p in (CODE_ARCHIVE, DATA_ARCHIVE, ARCHIVE_DIR / "heiscore_0.1.3.tar.gz")
                  if p.exists()]).to_csv(
        ROOT / "outputs/provenance/hei_package_manifest.csv", index=False)
    print(f"HEI package import: {len(combined):,} rows across {combined.cycle.nunique()} cycles")


if __name__ == "__main__":
    main()
