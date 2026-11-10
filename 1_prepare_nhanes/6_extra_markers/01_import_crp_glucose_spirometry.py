# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/6_extra_markers/01_import_crp_glucose_spirometry.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Reads the measurements outside the panel: C-reactive protein,
#               fasting glucose and insulin, and spirometry (with its quality
#               grade).
#
# METHODS       Results 7 candidates; CRP also feeds Phenotypic Age (Results 4).
# =============================================================================
"""Import the three candidate system measurements from raw NHANES files.

METHODS: see SPECIFICATION.md.

These measurements sit outside the 38-marker panel and outside the rebuilt
data-preparation pipeline, so they are imported here rather than there: they are
inputs to one diagnostic question, not to the panel, and adding them upstream
would place unused columns in tables that other modules hash.

Nothing is interpreted at import beyond what the survey itself defines. Source
variable names are preserved so the relationship with NHANES stays visible.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)

from nhanes_reader import raw_data_root, read_nhanes_xpt  # noqa: E402

OUTPUT = ROOT / "outputs/extra_markers"  # public_code: with the other outputs

# hsCRP. Two assay eras: LBXCRP through 2009-2010, LBXHSCRP from 2015-2016.
# 2011-2014 has no measurement at all and is simply absent from this table.
CRP_SOURCES = {
    "1999-2000": ("1999-2000/laboratory/LAB11.XPT", "LBXCRP", "LBXCRP"),
    "2001-2002": ("2001-2002/laboratory/L11_B.XPT", "LBXCRP", "LBXCRP"),
    "2003-2004": ("2003-2004/laboratory/L11_C.XPT", "LBXCRP", "LBXCRP"),
    "2005-2006": ("2005-2006/laboratory/CRP_D.XPT", "LBXCRP", "LBXCRP"),
    "2007-2008": ("2007-2008/laboratory/CRP_E.XPT", "LBXCRP", "LBXCRP"),
    "2009-2010": ("2009-2010/laboratory/CRP_F.XPT", "LBXCRP", "LBXCRP"),
    "2015-2016": ("2015-2016/laboratory/HSCRP_I.XPT", "LBXHSCRP", "LBXHSCRP"),
    "2017-2018": ("2017-2018/laboratory/HSCRP_J.XPT", "LBXHSCRP", "LBXHSCRP"),
}

# Fasting glucose and insulin. Insulin moves to its own file from 2013-2014.
GLUCOSE_SOURCES = {
    "1999-2000": ["1999-2000/laboratory/LAB10AM.XPT"],
    "2001-2002": ["2001-2002/laboratory/L10AM_B.XPT"],
    "2003-2004": ["2003-2004/laboratory/L10AM_C.XPT"],
    "2005-2006": ["2005-2006/laboratory/GLU_D.XPT"],
    "2007-2008": ["2007-2008/laboratory/GLU_E.XPT"],
    "2009-2010": ["2009-2010/laboratory/GLU_F.XPT"],
    "2011-2012": ["2011-2012/laboratory/GLU_G.XPT"],
    "2013-2014": ["2013-2014/laboratory/GLU_H.XPT", "2013-2014/laboratory/INS_H.XPT"],
    "2015-2016": ["2015-2016/laboratory/GLU_I.XPT", "2015-2016/laboratory/INS_I.XPT"],
    "2017-2018": ["2017-2018/laboratory/GLU_J.XPT", "2017-2018/laboratory/INS_J.XPT"],
}

SPIROMETRY_SOURCES = {
    "2007-2008": "2007-2008/examination/SPX_E.xpt",
    "2009-2010": "2009-2010/examination/SPX_F.xpt",
    "2011-2012": "2011-2012/examination/SPX_G.xpt",
}


def decode(series: pd.Series) -> pd.Series:
    """XPT character columns arrive as bytes; return plain strings."""
    if series.dtype == object:
        return series.apply(lambda v: v.decode().strip() if isinstance(v, bytes) else v)
    return series


def import_crp(root: Path) -> pd.DataFrame:
    frames = []
    for cycle, (path, variable, _) in CRP_SOURCES.items():
        raw = read_nhanes_xpt(root / path)
        if variable not in raw.columns:
            raise KeyError(f"{cycle}: {path} lacks {variable}")
        block = raw[["SEQN", variable]].copy()
        block.columns = ["SEQN", "crp_raw"]
        block.insert(1, "cycle", cycle)
        # The assay era is recorded, not corrected: 1999-2010 and 2015-2018 are
        # different measurements and are standardised within era downstream.
        block["assay_era"] = "LBXCRP_1999_2010" if variable == "LBXCRP" \
            else "LBXHSCRP_2015_2018"
        frames.append(block)
    return pd.concat(frames, ignore_index=True)


def import_glucose(root: Path) -> pd.DataFrame:
    frames = []
    for cycle, paths in GLUCOSE_SOURCES.items():
        merged = None
        for path in paths:
            raw = read_nhanes_xpt(root / path)
            keep = ["SEQN"] + [v for v in ("LBXGLU", "LBXIN", "PHAFSTHR", "WTSAF2YR")
                               if v in raw.columns]
            block = raw[keep].copy()
            merged = block if merged is None else merged.merge(
                block, on="SEQN", how="outer", suffixes=("", "_dup"))
        merged = merged[[c for c in merged.columns if not c.endswith("_dup")]]
        for variable in ("LBXGLU", "LBXIN", "PHAFSTHR", "WTSAF2YR"):
            if variable not in merged.columns:
                merged[variable] = float("nan")  # public_code: numeric empty column, so pandas keeps the column type now and in future versions
        merged.insert(1, "cycle", cycle)
        merged["fasting_hours_recorded"] = merged.PHAFSTHR.notna().any()
        frames.append(merged[["SEQN", "cycle", "LBXGLU", "LBXIN", "PHAFSTHR",
                              "WTSAF2YR", "fasting_hours_recorded"]])
    return pd.concat(frames, ignore_index=True)


def import_spirometry(root: Path) -> pd.DataFrame:
    frames = []
    for cycle, path in SPIROMETRY_SOURCES.items():
        raw = read_nhanes_xpt(root / path)
        block = raw[["SEQN", "SPXNFEV1", "SPXNFVC", "SPXNQEFF"]].copy()
        block["SPXNQEFF"] = decode(block.SPXNQEFF)
        block.insert(1, "cycle", cycle)
        frames.append(block)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    root = raw_data_root()

    crp = import_crp(root)
    crp.to_csv(OUTPUT / "candidate_crp_raw.csv", index=False)
    print(f"hsCRP        : {len(crp):,} rows, {crp.cycle.nunique()} cycles, "
          f"{crp.crp_raw.notna().sum():,} measured")
    print(f"               eras: {crp.groupby('assay_era').size().to_dict()}")
    print("               2011-2014 has no measurement and is absent by design")

    glucose = import_glucose(root)
    glucose.to_csv(OUTPUT / "candidate_glucose_raw.csv", index=False)
    print(f"glucose      : {len(glucose):,} rows, {glucose.cycle.nunique()} cycles, "
          f"{glucose.LBXGLU.notna().sum():,} measured")
    print(f"insulin      : {glucose.LBXIN.notna().sum():,} measured")
    print(f"fasting hours: recorded in "
          f"{glucose[glucose.PHAFSTHR.notna()].cycle.nunique()} of 10 cycles")

    spirometry = import_spirometry(root)
    spirometry.to_csv(OUTPUT / "candidate_spirometry_raw.csv", index=False)
    quality = spirometry.SPXNQEFF.value_counts(dropna=False).to_dict()
    print(f"spirometry   : {len(spirometry):,} rows, {spirometry.cycle.nunique()} cycles, "
          f"{spirometry.SPXNFEV1.notna().sum():,} with FEV1")
    print(f"               quality grades: {quality}")
    print(f"\nwritten -> {OUTPUT}")


if __name__ == "__main__":
    main()
