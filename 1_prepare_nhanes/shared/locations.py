"""Where the preparation scripts find their inputs and put their outputs.

Every script in 1_prepare_nhanes/ imports this file first. It does two jobs.

1. It names each folder once. The scripts build every path from these names,
   so moving the data means editing this file and nothing else.

2. It switches on a write guard. From then on the script may create or change
   files only inside `1_prepare_nhanes/outputs/` (and the system's temporary
   folder). An attempt to write anywhere in this code's repository or in the
   data folders below stops the script with "write guard: refused to write".
   This keeps the raw data, and any other copy of the prepared data, safe from
   a mistake in a path.

Setting the environment variable NHANES_ALLOW_WRITES to exactly "1" switches
the guard off. A normal run never needs that.
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# This folder
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]      # 1_prepare_nhanes/
CONFIG = ROOT / "config"                         # hand-written tables, read only
OUTPUTS = ROOT / "outputs"                       # everything the scripts make

# ---------------------------------------------------------------------------
# Data locations: set in data_locations.py, at the top of this code, the one
# file to edit before a first run.
# ---------------------------------------------------------------------------
sys.path.insert(0, str(ROOT.parent))
import data_locations  # noqa: E402

# NHANES 1999-2018 public files, one folder per cycle (1999-2000/, ...).
RAW_NHANES = Path(data_locations.NHANES_RAW)

# NCHS 2019 public-use linked mortality files (NHANES_<cycle>_MORT_2019_PUBLIC.dat).
MORTALITY = Path(data_locations.NHANES_MORTALITY)

# The CRAN archives heiscore_0.1.4.tar.gz and heiscore.data_0.0.1.tar.gz.
HEI_PACKAGES = Path(data_locations.HEI_PACKAGES)

# The NHANES reader and the mortality import read these two variables.
os.environ["NHANES_RAW_ROOT"] = str(RAW_NHANES)
os.environ["NHANES_MORTALITY_ROOT"] = str(MORTALITY)

OUTPUTS.mkdir(exist_ok=True)
(OUTPUTS / "provenance").mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Write guard
# ---------------------------------------------------------------------------
def _repository(start: Path) -> Path:
    """The git repository holding this code (or this folder's parent if none)."""
    for folder in [start, *start.parents]:
        if (folder / ".git").exists():
            return folder
    return start.parent


_PROTECTED = [_repository(ROOT), RAW_NHANES, MORTALITY, HEI_PACKAGES]
_ALLOWED = [OUTPUTS, Path(tempfile.gettempdir())]
_WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_APPEND | os.O_CREAT | os.O_TRUNC


def _inside(path: Path, folder: Path) -> bool:
    try:
        path.relative_to(folder.resolve())
        return True
    except ValueError:
        return False


def _check(target) -> None:
    if not isinstance(target, (str, bytes, os.PathLike)):
        return                                   # a file number, not a path
    try:
        path = Path(os.fsdecode(target)).resolve()
    except (OSError, ValueError):
        return
    if any(_inside(path, folder) for folder in _ALLOWED):
        return
    if any(_inside(path, folder) for folder in _PROTECTED):
        allowed = " or ".join(str(folder) for folder in _ALLOWED[:-1])   # the last is the temporary folder
        raise PermissionError(f"write guard: refused to write {path}. "
                              f"Preparation scripts write only inside {allowed}.")


def write_only_in(folder: Path) -> None:
    """From now on, let scripts write only inside `folder` (and the temporary folder),
    instead of 1_prepare_nhanes/outputs/. 2_prepare_nhanes_iii uses this for its own outputs/."""
    _ALLOWED[:-1] = [Path(folder)]


def protect(folder: Path) -> None:
    """Add a folder that scripts must never write in (for example the NHANES III raw files)."""
    _PROTECTED.append(Path(folder))


def _audit(event: str, args: tuple) -> None:
    if event == "open":
        target, mode, flags = args
        if (isinstance(mode, str) and any(c in mode for c in "wax+")) or \
                (mode is None and isinstance(flags, int) and flags & _WRITE_FLAGS):
            _check(target)
    elif event in ("os.remove", "os.rmdir", "os.mkdir", "os.truncate", "shutil.rmtree", "os.chmod"):
        _check(args[0])
    elif event in ("os.rename", "shutil.move"):
        _check(args[0])
        _check(args[1])
    elif event == "shutil.copyfile":
        _check(args[1])                          # copying from a protected folder is a read


if os.environ.get("NHANES_ALLOW_WRITES") != "1" and not getattr(sys, "_prepare_write_guard", False):
    sys.addaudithook(_audit)                     # cannot be removed once added, by design
    sys._prepare_write_guard = True
