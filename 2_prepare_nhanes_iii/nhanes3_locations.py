"""Where the NHANES III preparation scripts find their inputs and put their outputs.

Every script in 2_prepare_nhanes_iii/ imports this file first. It does the same
two jobs as 1_prepare_nhanes/shared/locations.py, and shares its write guard:

1. It names each folder once, so moving the data means editing this file only.

2. It switches on the write guard and points it at this folder's outputs/:
   a script may write only there. Any other write inside this repository or
   the data folders stops the script with "write guard: refused to write".

Run 1_prepare_nhanes/ first: steps 2 and 4 read two of its tables, and step 2
reuses its transformation code.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# This folder and the modern preparation it builds on
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent                     # 2_prepare_nhanes_iii/
CONFIG = HERE / "config"                                   # hand-written tables, read only
OUTPUTS = HERE / "outputs"                                 # everything the scripts make
PREPARE_NHANES = HERE.parent / "1_prepare_nhanes"
MODERN_CONFIG = PREPARE_NHANES / "config"

sys.path.insert(0, str(PREPARE_NHANES / "shared"))
import locations as _guard  # noqa: E402  (switches on the shared write guard)

# ---------------------------------------------------------------------------
# Data location: set in data_locations.py, at the top of this code (already
# importable: locations.py above put that folder on the path).
# ---------------------------------------------------------------------------
import data_locations  # noqa: E402

# NHANES III public files: adult.dat/.sas, exam.dat/.sas, lab.dat/.sas and the
# NCHS 2019 linked mortality file NHANES_III_MORT_2019_PUBLIC.dat.
RAW_NHANES3 = Path(data_locations.NHANES3_RAW)

_guard.write_only_in(OUTPUTS)
_guard.protect(RAW_NHANES3)
OUTPUTS.mkdir(exist_ok=True)


def modern_transformation_rule():
    """The modern transformation script (1_prepare_nhanes/7_marker_panel/02), loaded as a module.

    Loading it only defines its rule (candidates, limits, `diagnostics`,
    `evaluate_marker`); its main() does not run.
    """
    path = PREPARE_NHANES / "7_marker_panel" / "02_select_transformations.py"
    spec = importlib.util.spec_from_file_location("modern_select_transformations", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
