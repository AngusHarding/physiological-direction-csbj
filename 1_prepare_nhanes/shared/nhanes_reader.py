# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/shared/nhanes_reader.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Reads NHANES .XPT files and restores the zeros that the SAS
#               transport format stores as tiny numbers.
#
# METHODS       Supp. Section Zero, "Missing-value codes".
# =============================================================================
"""Shared, minimally interpretive NHANES XPT reader for data preparation.

Beginner's map
--------------
An NHANES ``.XPT`` file is a SAS transport file: a table stored in a format
commonly used for public health datasets.  ``pandas`` reads that table into a
``DataFrame``.  A DataFrame is like an Excel worksheet held in Python: rows are
participants/records and columns are variables.

This module deliberately does only three small jobs:

1. find the local folder containing the raw NHANES files;
2. read one XPT file without changing its scientific meaning;
3. check that ``SEQN`` is a valid unique participant identifier.

Study exclusions do *not* belong here.  Keeping the reader simple lets us reuse
it without accidentally applying a cohort rule twice.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


# pandas sometimes represents a literal SAS zero as approximately 5.4e-79.
# Any finite number this close to zero is restored to the literal value 0.0.
SAS_ZERO_THRESHOLD = 1e-50


def raw_data_root() -> Path:
    """Return the local raw-data root without embedding it in public code."""
    # Environment variables let each user keep raw data wherever they choose.
    # The public repository therefore does not reveal or depend on Angus's path.
    value = os.environ.get("NHANES_RAW_ROOT")
    if not value:
        raise RuntimeError(
            "NHANES_RAW_ROOT is not set. See config/raw_data_paths.example.yaml."
        )
    # Path() gives us safe path operations. resolve() produces an absolute path.
    root = Path(value).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"NHANES raw-data root does not exist: {root}")
    return root


def read_nhanes_xpt(path: Path) -> pd.DataFrame:
    """Read an XPT file and restore literal zeros misdecoded by pandas.

    True SAS missing values remain NaN. Numeric values with absolute magnitude
    below 1e-50 include exact zero and the observed XPT zero representation;
    these are restored to 0.0, matching the certified historical reader.
    """
    # One call turns the complete SAS table into a pandas DataFrame.
    frame = pd.read_sas(path)
    # Only numeric columns can contain the unusual tiny-number representation.
    numeric = frame.select_dtypes("number")
    # mask(condition, replacement) replaces cells where condition is True.
    frame[numeric.columns] = numeric.mask(numeric.abs() < SAS_ZERO_THRESHOLD, 0.0)
    return frame


def validate_participant_key(frame: pd.DataFrame, label: str) -> None:
    """Fail when a materialised participant table has missing or duplicate SEQN."""
    if "SEQN" not in frame.columns:
        raise KeyError(f"{label}: SEQN column is absent")
    # isna() marks missing cells; duplicated(..., keep=False) marks every member
    # of a duplicate group. sum() counts True values because True behaves as 1.
    missing = int(frame["SEQN"].isna().sum())
    duplicates = int(frame.duplicated("SEQN", keep=False).sum())
    if missing or duplicates:
        raise ValueError(
            f"{label}: invalid participant key (missing={missing}, "
            f"rows_in_duplicate_groups={duplicates})"
        )
