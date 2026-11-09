"""Run every Results section, then check the outputs against the paper.

This script estimates nothing itself. It runs the seven section scripts in
paper order with agreed arguments, then runs check_against_paper.py for all
seven sections together.

Two modes:

    python -B make_all.py              # SHORT (default): every point estimate, table and figure,
                                       # plus a few bootstrap draws and evaluations per section that are
                                       # compared one by one with the archive. About 7-10 minutes.
    python -B make_all.py --full       # FULL: every bootstrap and all 472 panel-selection evaluations,
                                       # so every interval is recomputed. Many hours; see FULL below.
    python -B make_all.py --dry-run    # print the commands either mode would run, and run nothing
    python -B make_all.py --help       # this help; runs nothing

Why the short mode is the default. The full bootstraps take hours. The short
mode uses the same code and the same random numbers for the draws it does
run; if those draws match the archived draws one by one and the point
estimates match, the code is doing what the reviewed analysis did. It does
not by itself recompute the intervals: the check reports them as pending
(or as verified from the archived draws, where those exist).

Before anything starts, the runner asks each section script for its list of
options and stops if an option it needs is missing, so a mistake in the plan
fails in seconds rather than after an hour of computing. A section that
fails stops the run.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time

import settings as S

# One script per Results section, in paper order, with the arguments for each mode.
# SHORT: draws/evaluations compared with the archive (expected/<section>_draws.csv).
#   s4 runs the first two five-versus-five and two NHANES III half-split evaluations, plus
#   F09, the NHANES III held-out fold whose full-panel gain is negative (it tests the rule
#   that leaves the retained-gain ratio missing there).
# FULL: every bootstrap draw and every evaluation. Approximate times on 12 cores, from the
#   archived runs: s1 4 min; s2 about an hour; s3 40-90 min; s4 30 min plus the benchmark
#   bootstrap; s5 hours; s6 and s7 hours.
PLAN = {
    "s1_sex_gap.py":       {"short": ["--check-draws", "2"],               "full": ["--bootstrap"]},
    "s2_reference.py":     {"short": ["--check-draws", "2"],               "full": ["--bootstrap"]},
    "s3_healthy_range.py": {"short": ["--check-draws", "2"],               "full": ["--bootstrap"]},
    "s4_panel.py":         {"short": ["--evaluations", "2", "--only", "F09", "--check-draws", "2"],
                            "full": ["--bootstrap"]},
    "s5_transport.py":     {"short": ["--check-draws", "3"],               "full": ["--bootstrap"]},
    "s6_screen.py":        {"short": ["--check-draws", "2"],               "full": ["--bootstrap"]},
    "s7_functional.py":    {"short": ["--check-draws", "1"],               "full": ["--bootstrap"]},
}


def options_of(script: str) -> str:
    """The script's own --help text (argparse exits before any data are read)."""
    done = subprocess.run([sys.executable, "-B", str(S.PIPELINE / script), "--help"],
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    if done.returncode != 0:
        raise SystemExit(f"{script} --help failed:\n{done.stderr}")
    return done.stdout


def check_plan(commands: list[tuple[str, list[str]]], workers: int | None) -> None:
    """Stop before running anything if any section does not accept an option the plan uses."""
    problems = []
    for script, args in commands:
        help_text = options_of(script)
        for option in [a for a in args if a.startswith("--")]:
            if option not in help_text:
                problems.append(f"{script} does not accept {option}")
        if workers is not None and "--workers" not in help_text:
            problems.append(f"{script} does not accept --workers")
    if problems:
        raise SystemExit("make_all: the run plan does not match the section scripts:\n  " + "\n  ".join(problems))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run every Results section, then check against the paper.")
    parser.add_argument("--full", action="store_true",
                        help="recompute every bootstrap and all panel-selection evaluations (many hours)")
    parser.add_argument("--workers", type=int, default=None,
                        help="processes per section (default: each script's own default)")
    parser.add_argument("--dry-run", action="store_true", help="print the commands and stop")
    args = parser.parse_args()                    # unknown options stop here, before any work

    mode = "full" if args.full else "short"
    workers = ["--workers", str(args.workers)] if args.workers else []
    commands = [(script, plan[mode] + workers) for script, plan in PLAN.items()]
    check_plan(commands, args.workers)
    print(f"make_all: {mode.upper()} run. Commands, in order:")
    for script, extra in commands:
        print(f"  python -B {script} {' '.join(extra)}".rstrip())
    print("  python -B check_against_paper.py --all")
    if args.dry_run:
        return 0

    for script, extra in commands:
        print(f"\n=== {script} {' '.join(extra)} ===", flush=True)
        start = time.time()
        subprocess.run([sys.executable, "-B", str(S.PIPELINE / script), *extra], check=True)
        print(f"=== {script}: {(time.time() - start) / 60:.1f} min ===", flush=True)
    print("\n=== checking against the paper ===", flush=True)
    return subprocess.run([sys.executable, "-B", str(S.PIPELINE / "check_against_paper.py"), "--all"]).returncode


if __name__ == "__main__":
    sys.exit(main())
