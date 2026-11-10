# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/run_all.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Runs every NHANES 1999-2018 preparation script in folder and
#               file order, each in its own process.
#
# METHODS       None; it calculates nothing itself.
# =============================================================================
"""Run every preparation script for NHANES 1999-2018, in order.

    python -B run_all.py              # all 46 scripts
    python -B run_all.py 7_marker_panel 8_adverse_distance   # only these folders

Folders run in number order (1_cohort, 2_blood_tests, ...) and the scripts in
each folder in number order, each in its own Python process, exactly as if typed
one by one from this folder. The first script that fails stops the run.
"""
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> None:
    folders = sorted(f for f in HERE.iterdir() if f.is_dir() and f.name[0].isdigit())
    if len(sys.argv) > 1:
        unknown = set(sys.argv[1:]) - {f.name for f in folders}
        if unknown:
            sys.exit(f"no such folder: {', '.join(sorted(unknown))}")
        folders = [f for f in folders if f.name in sys.argv[1:]]
    scripts = [s for f in folders for s in sorted(f.glob("[0-9][0-9]_*.py"))]
    start = time.time()
    for n, script in enumerate(scripts, 1):
        name = script.relative_to(HERE).as_posix()
        print(f"\n[{n}/{len(scripts)}] {name}", flush=True)
        if subprocess.run([sys.executable, "-B", name], cwd=HERE).returncode != 0:
            sys.exit(f"stopped: {name} failed")
    print(f"\nall {len(scripts)} scripts done in {(time.time() - start) / 60:.1f} min")


if __name__ == "__main__":
    main()
