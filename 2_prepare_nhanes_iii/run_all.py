# =============================================================================
# READ ME FIRST: 2_prepare_nhanes_iii/run_all.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Runs the five NHANES III preparation steps in order.
#
# METHODS       None; it calculates nothing itself.
# =============================================================================
"""Run the five NHANES III preparation steps, in order.

    python -B run_all.py

Each step runs in its own Python process, exactly as if typed one by one from
this folder. The first step that fails stops the run. Run 1_prepare_nhanes/ first.
"""
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> None:
    scripts = sorted(HERE.glob("[0-9][0-9]_*.py"))
    start = time.time()
    for n, script in enumerate(scripts, 1):
        print(f"\n[{n}/{len(scripts)}] {script.name}", flush=True)
        if subprocess.run([sys.executable, "-B", script.name], cwd=HERE).returncode != 0:
            sys.exit(f"stopped: {script.name} failed")
    print(f"\nall {len(scripts)} steps done in {(time.time() - start) / 60:.1f} min")


if __name__ == "__main__":
    main()
