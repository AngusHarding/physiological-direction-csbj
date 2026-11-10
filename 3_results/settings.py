# =============================================================================
# READ ME FIRST: settings.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Holds every fixed choice the analysis depends on (folders, age
#               floor, accident code, survey cycles, age splines, seeds, draw
#               counts), and a guard that stops any script writing outside
#               outputs/. It calculates nothing.
#
# METHODS       "Participants and endpoints": the age floor of 52 ("Why 52") and
#               the accident code behind the principal endpoint.
#               "Statistical analysis": the age spline used for sex differences
#               and the 500 bootstrap draws behind Results 1 and 2.
#
# USED BY       Every Results script.
#
# WHERE TO LOOK
#   AGE_MIN = 52           Methods, "Mortality cohorts": "participants aged 52-79"
#   ACCIDENT_CODE = 4      Methods, "Endpoints": "other than an identified accident"
#   CYCLES                 the ten survey cycles, 1999-2000 to 2017-2018
#   ERAS[...]["age_spline"] Methods, "Sex differences": "an age spline"
#   ERAS[...]["period"]    "survey cycle or phase"
#   BOOTSTRAP_DRAWS = 500  the participant bootstrap of Results 1 and 2
# =============================================================================
"""Settings shared by every script: where the data are, and the fixed choices of the analysis.

Read this file first. Every folder and every constant the pipeline uses is
defined here, once. If a number in the paper depends on a choice (an age
limit, a random seed, a number of bootstrap draws), the choice is below.

This file defines no estimate; it holds the choices the estimates depend on.

Folders
-------
PREP       Prepared NHANES 1999-2018 tables: the outputs of ../1_prepare_nhanes/,
           which turns the raw NHANES files into one row per participant.
PREP_CONFIG  The hand-written tables of ../1_prepare_nhanes/config/.
III_PREP   Prepared NHANES III (1988-1994) tables: the outputs of ../2_prepare_nhanes_iii/,
           including the questionnaire and examination variables and the linked
           mortality file that the Results use.
EXTRACTS   Three prepared candidate-measurement tables (CRP, fasting glucose,
           spirometry), made by 1_prepare_nhanes/6_extra_markers/.
LOCKED_TEN The development ranking that fixes the ten-marker panel used by
           Results 7: the ranking of the reviewed Results 4 run, kept in inputs/
           (a full run of s4_panel.py remakes it as outputs/s4/locked_ten_modern.csv).
REPO       The git repository that holds this code. Output files record their
           paths relative to it, and the write guard protects it.
           (LOCAL_REPO is the same folder, under the name runs.py uses.)
RENDERER   The drawing code for the accepted Figures 1 and 2 (render_figures.py,
           this folder); it draws results the pipeline hands it and reads no
           estimates of its own.
OUTPUTS    Where this pipeline writes. One subfolder per Results section.

The write guard
---------------
Importing this file switches on a guard. Inside REPO and inside LOCAL_REPO, a
script may create, change or delete files only inside OUTPUTS (and in the
system's temporary folder); any other write there stops the script with an
error. So rerunning the pipeline cannot damage the input data, the archived
analysis or the pipeline's own code, wherever the data are read from. Writes
outside both repositories (for example matplotlib's font cache in the user's
home folder) are not the pipeline's business and are allowed.
Setting NHANES_ALLOW_WRITES=1 (exactly "1") switches the guard off; there is
no reason to do so for a normal run.

Inputs read
-----------
The same hook records every data file opened for reading inside REPO or
LOCAL_REPO (Python source files excepted). runs.py hashes those files into
each section's run record, so a stored bootstrap can be matched to the exact
data that produced it. Worker processes keep their own list, so while a
section runs (runs.start sets NHANES_READS_LOG) every process, parent and
workers alike, also appends each newly read file to reads_<process id>.txt in
that folder; runs.all_reads() gathers them (R120-F14).
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Folders
# ---------------------------------------------------------------------------
PIPELINE = Path(__file__).resolve().parent                 # 3_results/
PUBLIC = PIPELINE.parent                                   # public_code: the three stages
LOCAL_REPO = next((f for f in PUBLIC.parents if (f / ".git").exists()), PUBLIC).resolve()  # public_code
REPO = LOCAL_REPO                                          # public_code: one repository holds code and prepared data
PREP = PUBLIC / "1_prepare_nhanes" / "outputs"             # public_code
PREP_CONFIG = PUBLIC / "1_prepare_nhanes" / "config"       # public_code
III_PREP = PUBLIC / "2_prepare_nhanes_iii" / "outputs"     # public_code
EXTRACTS = PREP / "extra_markers"                          # public_code
INPUTS = PIPELINE / "inputs"
LOCKED_TEN = INPUTS / "full_development_ranking.csv"       # public_code
RENDERER = PIPELINE / "render_figures.py"                  # public_code
OUTPUTS = PIPELINE / "outputs"

# ---------------------------------------------------------------------------
# Cohorts and endpoint
# ---------------------------------------------------------------------------
AGE_MIN = 52                  # analysis cohorts: examined at age 52 or older
ACCIDENT_CODE = 4             # NCHS leading-cause code for accidents
# The ten continuous-NHANES survey cycles, in time order.
CYCLES = ("1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
          "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018")

# ---------------------------------------------------------------------------
# The two survey eras
# ---------------------------------------------------------------------------
# Each era was analysed with its own conventions in the reviewed analysis, and
# the pipeline keeps them so that it reproduces the reviewed numbers exactly.
#
#   period     the column that marks when a participant was examined: the
#              two-year survey cycle (1999-2018) or the survey phase (NHANES
#              III ran in two phases, 1988-91 and 1991-94).
#   age_spline the flexible age term used when comparing women and men.
#   cox        which Cox-model software fits the hazard ratios (see models.py).
#   seed       the random seed for the participant bootstrap.
ERAS = {
    "modern": {
        "label": "NHANES 1999–2018",
        "age": "RIDAGEYR",
        "period": "cycle",
        "age_spline": "bs(RIDAGEYR, df=4, degree=3, include_intercept=False)",
        "cox": "statsmodels, stratified by cycle",
        "seed": 20260803,
    },
    "iii": {
        "label": "NHANES III",
        "age": "age",
        "period": "phase",
        "age_spline": 'cr(age, df=4, constraints="center")',
        "cox": "lifelines, unstratified",
        "seed": 20261015,
    },
}

BOOTSTRAP_DRAWS = 500         # participant-bootstrap draws behind the Results 1 and 2 intervals
# (Results 3-7 resample differently and set their own draw counts and seeds in their scripts.)

# ---------------------------------------------------------------------------
# Write guard (see the module docstring)
# ---------------------------------------------------------------------------
_TEMP = Path(tempfile.gettempdir()).resolve()
_WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_APPEND | os.O_CREAT | os.O_TRUNC


def _inside(path: Path, folder: Path) -> bool:
    try:
        path.relative_to(folder)
        return True
    except ValueError:
        return False


_PROTECTED = (LOCAL_REPO,)           # public_code: this repository (code and prepared data)
READS: set[Path] = set()             # data files opened for reading (see "Inputs read" above)


def _resolve(target) -> Path | None:
    if not isinstance(target, (str, bytes, os.PathLike)):
        return None                             # a file number, not a path
    try:
        return Path(os.fsdecode(target)).resolve()
    except (OSError, ValueError):
        return None


_GUARD_ON = os.environ.get("NHANES_ALLOW_WRITES") != "1"


def _check(target) -> None:
    if not _GUARD_ON:
        return
    path = _resolve(target)
    if path is None or _inside(path, OUTPUTS) or _inside(path, _TEMP):
        return
    if any(_inside(path, root) for root in _PROTECTED):
        raise PermissionError(f"write guard: refused to write {path}. Inside the repositories the "
                              f"pipeline writes only inside {OUTPUTS}.")


READS_ENV = "NHANES_READS_LOG"     # folder for each process's list of files read (set by runs.start)


def _record_read(target) -> None:
    path = _resolve(target)
    if path is None or path.suffix in (".py", ".pyc") or _inside(path, OUTPUTS) or path in READS:
        return
    if any(_inside(path, root) for root in _PROTECTED) and path.is_file():
        READS.add(path)
        folder = os.environ.get(READS_ENV)
        if folder:                              # shared with the parent: workers' reads reach the run record
            fd = os.open(os.path.join(folder, f"reads_{os.getpid()}.txt"), os.O_WRONLY | os.O_APPEND | os.O_CREAT)
            try:
                os.write(fd, (str(path) + "\n").encode("utf-8"))
            finally:
                os.close(fd)


def _audit(event: str, args: tuple) -> None:
    if event == "open":
        target, mode, flags = args
        if (isinstance(mode, str) and any(c in mode for c in "wax+")) or \
                (mode is None and isinstance(flags, int) and flags & _WRITE_FLAGS):
            _check(target)
        else:
            _record_read(target)
    elif event in ("os.remove", "os.rmdir", "os.mkdir", "os.truncate", "shutil.rmtree"):
        _check(args[0])
    elif event in ("os.rename", "shutil.move"):
        _check(args[0])
        _check(args[1])


if not getattr(sys, "_nhanes_write_guard", False):
    sys.addaudithook(_audit)                    # cannot be removed once added, by design
    sys._nhanes_write_guard = True


def output_folder(section: str) -> Path:
    """The folder a section writes into, created if needed."""
    folder = OUTPUTS / section
    folder.mkdir(parents=True, exist_ok=True)
    return folder
