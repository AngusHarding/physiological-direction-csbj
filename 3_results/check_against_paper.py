"""Check every number the pipeline made against the reviewed analysis and the paper.

This script estimates nothing and writes nothing except a cache of file
fingerprints (outputs/.file_hashes.json). For each row of
expected/expected_values.csv it asks two questions:

  1. Reproduced?  Does the pipeline's value (and its 95% interval) equal the
     reviewed value to within 1e-8 (scaled by the value's size)? Most reviewed
     values were copied at full precision from the archived analysis or the
     reviewed revision packages; some cohort counts and table cells were taken
     from the printed manuscript, and a few rows carry a looser tolerance. The
     source of each (column `source`, a file in the authors' working archive)
     and the reason for any tolerance are recorded with it.
  2. As printed?  Rounded the way the paper rounds it, does the pipeline's value
     read exactly as the manuscript prints it? Rows marked `paper_issue` are
     known slips in the manuscript; they are listed separately, not as failures.

Tables. Each table the pipeline writes is compared with the reviewed table as
text, cell by cell, labels and footnote marks included. The only differences
allowed are those listed, with a reason, in expected/table_format_exceptions.csv;
they are reported separately as formatting differences. A cell whose value is
right but whose 95% interval is not there yet (no full bootstrap has been run)
counts as pending, like any other interval.

Bootstrap draws. A full bootstrap can take hours. Each section's short check
runs a few draws (`--check-draws N`; for Results 4, a few evaluations) and
writes them to outputs/<section>/check_draws.csv, and its run record
(outputs/<section>/run.json) says how many rows it wrote. Those draws must
match the archived draws in expected/<section>_draws.csv one by one: every
draw present in the archive, each exactly once, every archived quantity
present. Quantities that have no archived draws are listed as unverified.
Matching draws show that this code, with these random numbers, gives the
reviewed draws; they do not recompute the intervals, which need every draw.

Intervals and summaries. An interval the pipeline has not computed yet is
counted as pending. It is counted as verified on the archived results instead
when the authors' verifier has applied the pipeline's own summary rule to the
complete archived draws and left a receipt (expected/interval_verification.csv, with
expected/interval_verification_sources.json), and that receipt still holds
(R120-F13):
  - its layout is right and each interval or summary appears once;
  - its bounds (or value) equal the reviewed numbers and read as printed, and
    it used as many draws (or evaluations) as a full run does (the section's
    run.json, `full_draws`); a receipt that fails here fails the section;
  - the pipeline modules and verifier scripts it used have the same syntax-tree
    fingerprints now, and the archived files it read the same SHA-256. If not
    (for example where the archive is not available, as in the public copy of
    this code), its intervals count as pending: a full run recomputes them.
Results 4's selection summaries (which need all 472 evaluations) are treated
the same way in a partial run. Neither counts as a reproduction by a full run.

    python -B check_against_paper.py            # every section that has outputs
    python -B check_against_paper.py s1 s4      # some sections
    python -B check_against_paper.py --all      # all seven sections; any missing section fails

The exit code is 0 when everything passes, 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

import runs
import settings as S

TOLERANCE = 1e-8
DRAW_TOLERANCE = 1e-6           # draw-by-draw check; exceptions in expected/draw_tolerances.csv
EXPECTED = S.PIPELINE / "expected"
SECTIONS = ("s1", "s2", "s3", "s4", "s5", "s6", "s7")
TABLES = {"s1": ["table1"], "s2": ["table2"], "s6": ["table4"]}      # tables compared cell by cell
RECEIPT_COLUMNS = ["section", "key", "era", "kind", "draws", "low", "high", "value", "verified_by"]


def as_printed(value: float, decimals: int) -> str:
    """Format like the paper: thousands separators, a true minus sign."""
    text = f"{value:,.{decimals}f}"
    return text.replace("-", "−")


def _optional_csv(name: str) -> pd.DataFrame:
    path = EXPECTED / name
    return pd.read_csv(path, dtype=str, keep_default_na=False) if path.is_file() else pd.DataFrame()


# ---------------------------------------------------------------------------
# Receipts of the interval rule applied to the archived draws
# ---------------------------------------------------------------------------
def load_receipts() -> tuple[pd.DataFrame, dict, list]:
    """(receipt rows, their sources by section, layout failures)."""
    path, sources_path = EXPECTED / "interval_verification.csv", EXPECTED / "interval_verification_sources.json"
    empty = pd.DataFrame(columns=RECEIPT_COLUMNS)
    if not path.is_file():
        return empty, {}, []
    rows = pd.read_csv(path, keep_default_na=False, na_values=[""],
                       dtype={"section": str, "key": str, "era": str, "kind": str, "verified_by": str})
    if list(rows.columns) != RECEIPT_COLUMNS:
        return empty, {}, [f"interval_verification.csv: columns {list(rows.columns)}, expected {RECEIPT_COLUMNS}"]
    failures = []
    if rows.duplicated(["section", "key", "era"]).any():
        failures.append("interval_verification.csv: some intervals or summaries appear more than once")
    if not rows.kind.isin(["interval", "summary"]).all() or not rows.section.isin(SECTIONS).all():
        failures.append("interval_verification.csv: unknown section or kind")
    try:
        sources = json.loads(sources_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        sources = {}
    return (empty if failures else rows), sources, failures


def out_of_date(entry: dict | None) -> list[str]:
    """Why a section's receipt no longer describes the current code and archived files (empty if it does)."""
    if not entry or not entry.get("code") or not entry.get("archive") or not entry.get("verifier"):
        return ["no record of the code and archived files it used"]
    reasons = []
    for group in ("verifier", "code"):
        for name, digest in entry[group].items():
            path = runs._resolve_relative(name)
            if not path.is_file():
                reasons.append(f"{name} is missing")
            elif runs._syntax_hash(path) != digest:
                reasons.append(f"{name} has changed")
    present = {name: runs._resolve_relative(name) for name in entry["archive"]}
    reasons += [f"{name} is missing" for name, path in present.items() if not path.is_file()]
    now = runs.file_hashes([path for path in present.values() if path.is_file()])
    reasons += [f"{name} has changed" for name, digest in entry["archive"].items() if name in now and now[name] != digest]
    return reasons


def check_receipts(section: str, receipts: pd.DataFrame, sources: dict, expected: pd.DataFrame,
                   record: dict) -> tuple[list, set, str | None]:
    """(failures, (section, key, era) verified on the archived results, note if out of date)."""
    mine = receipts[receipts.section.eq(section)]
    if mine.empty:
        return [], set(), None
    stale = out_of_date(sources.get(section))
    if stale:
        more = f" and {len(stale) - 1} more" if len(stale) > 1 else ""
        return [], set(), (f"{len(mine)} intervals or summaries were verified by the authors on the archived "
                           f"results (expected/interval_verification.csv), which cannot be confirmed here "
                           f"({stale[0]}{more}); they count as pending a full run")
    targets = expected.set_index(["key", "era"])
    full = record.get("full_draws") or {}
    failures = []
    for r in mine.itertuples():
        label = f"{section} {r.key} [{r.era}] archive verification"
        if (r.key, r.era) not in targets.index:
            failures.append(f"{label}: not a number in expected_values.csv")
            continue
        want = targets.loc[(r.key, r.era)]
        tolerance = want.tolerance if pd.notna(want.tolerance) else TOLERANCE
        if r.kind == "interval":
            parts = (("low", r.low, want.reviewed_low, want.printed_low),
                     ("high", r.high, want.reviewed_high, want.printed_high))
            digits = want.ci_decimals if pd.notna(want.ci_decimals) else want.decimals
        else:
            parts, digits = (("value", r.value, want.reviewed, want.printed),), want.decimals
        known_slip = pd.notna(want.paper_issue) and want.paper_issue != ""
        for part, got, reviewed, printed in parts:
            if pd.isna(got) or pd.isna(reviewed) and (pd.isna(printed) or printed == ""):
                failures.append(f"{label}: no {part} to compare")
                continue
            if pd.notna(reviewed) and not abs(got - reviewed) <= tolerance * max(1.0, abs(reviewed)):
                failures.append(f"{label}: {part} {got!r} differs from reviewed {reviewed!r}")
            if pd.notna(printed) and printed != "" and not known_slip and as_printed(got, int(digits)) != str(printed):
                failures.append(f"{label}: {part} prints as {as_printed(got, int(digits))}, the paper says {printed}")
        if full.get(r.kind) is None:
            failures.append(f"{label}: the run record does not say how many draws a full run uses")
        elif r.draws != full[r.kind]:
            failures.append(f"{label}: made from {r.draws:g} draws or evaluations; a full run uses {full[r.kind]}")
    verified = set() if failures else set(zip(mine.section, mine.key, mine.era))
    return failures, verified, None


# ---------------------------------------------------------------------------
def check_numbers(section: str, expected: pd.DataFrame, corrections: list, verified: set) -> tuple[list, dict]:
    """(failures, counts) for one section's quoted numbers."""
    made = pd.read_csv(S.OUTPUTS / section / "numbers.csv")
    if made.duplicated(["key", "era"]).any():
        return [f"{section}: numbers.csv has duplicated keys"], {}
    made = made.set_index(["key", "era"])
    partial = any(k == "partial_run" for k, _ in made.index)
    failures, counts = [], {"pending": 0, "verified_on_archive": 0, "partial": 0}
    for r in expected.itertuples():
        label = f"{section} {r.key} [{r.era}]"
        if (r.key, r.era) not in made.index:
            if partial and (section, r.key, r.era) in verified:   # summary code checked on archived results
                counts["verified_on_archive"] += 1
            elif partial:                                 # a check run covered only some evaluations
                counts["partial"] += 1
            else:
                failures.append(f"{label}: not produced")
            continue
        got = made.loc[(r.key, r.era)]
        tolerance = r.tolerance if pd.notna(getattr(r, "tolerance", np.nan)) else TOLERANCE
        has_interval = pd.notna(r.reviewed_low) or (pd.notna(r.printed_low) and r.printed_low != "")
        if has_interval and pd.isna(got["low"]):
            if (section, r.key, r.era) in verified:
                counts["verified_on_archive"] += 1
            else:
                counts["pending"] += 1
            parts, printed = (("value", r.reviewed),), (("value", r.printed),)
        else:
            parts = (("value", r.reviewed), ("low", r.reviewed_low), ("high", r.reviewed_high))
            printed = (("value", r.printed), ("low", r.printed_low), ("high", r.printed_high))
        for part, want in parts:
            if pd.notna(want) and not abs(got[part] - want) <= tolerance * max(1.0, abs(want)):
                failures.append(f"{label}: {part} {got[part]!r} differs from reviewed {want!r}")
        for part, want in printed:
            if pd.notna(want) and want != "":
                digits = r.decimals
                if part != "value" and pd.notna(getattr(r, "ci_decimals", np.nan)):
                    digits = r.ci_decimals                 # a few intervals are printed less precisely
                shown = as_printed(got[part], int(digits)) if pd.notna(got[part]) else "missing"
                if shown != str(want):
                    message = f"{label}: {part} prints as {shown}, the paper says {want}"
                    issue = getattr(r, "paper_issue", np.nan)
                    (corrections if pd.notna(issue) and issue != "" else failures).append(message)
    return failures, counts


def check_tables(section: str, exceptions: pd.DataFrame) -> tuple[list, list, int]:
    """(failures, formatting differences allowed by expected/table_format_exceptions.csv,
    number of cells whose value is right but whose interval is still pending)."""
    failures, formatting, pending = [], [], 0
    for stem in TABLES.get(section, []):
        made = pd.read_csv(S.OUTPUTS / section / f"{stem}.csv", dtype=str, encoding="utf-8-sig", keep_default_na=False)
        want = pd.read_csv(EXPECTED / f"{section}_{stem}.csv", dtype=str, encoding="utf-8-sig", keep_default_na=False)
        if made.shape != want.shape or list(made.columns) != list(want.columns):
            failures.append(f"{section} {stem}: layout differs from the reviewed table "
                            f"({made.shape} vs {want.shape}, or different column headings)")
            continue
        allowed = {}
        if not exceptions.empty:
            for e in exceptions[(exceptions.section == section) & (exceptions.table == stem)].itertuples():
                allowed[(int(e.row), e.column)] = (e.made, e.reviewed, e.reason)
        for i in range(len(want)):
            for col in want.columns:
                a, b = made.at[i, col], want.at[i, col]
                if a == b:
                    continue
                where = f"{section} {stem}: row {i + 1} ('{want.iat[i, 0]}') / {col}"
                if allowed.get((i + 1, col), (None, None))[:2] == (a, b):
                    formatting.append(f"{where}: {a!r} vs reviewed {b!r} ({allowed[(i + 1, col)][2]})")
                elif a and b.startswith(a + " (") and b.endswith(")") and "(" not in a:
                    pending += 1                       # right value; its interval needs a full run
                else:
                    failures.append(f"{where}: made {a!r}, reviewed {b!r}")
    return failures, formatting, pending


def check_draws(section: str, record: dict, tolerances: pd.DataFrame) -> tuple[list, int, list]:
    """(failures, draws compared, quantities without archived draws) for one section's short check."""
    info = record.get("check_draws")
    made_file, want_file = S.OUTPUTS / section / "check_draws.csv", EXPECTED / f"{section}_draws.csv"
    if not info:
        return ([f"{section}: check_draws.csv exists but the run record lists no check draws"]
                if made_file.is_file() else []), 0, []
    if not made_file.is_file():
        return [f"{section}: the run record lists {info['rows']} check draws but check_draws.csv is missing"], 0, []
    if not want_file.is_file():
        return [f"{section}: no archived draws to compare with ({want_file.name})"], 0, []
    keys = info["keys"]
    made, want = pd.read_csv(made_file), pd.read_csv(want_file)
    failures = []
    if len(made) != info["rows"]:
        failures.append(f"{section} draws: {len(made)} rows, the run record says {info['rows']}")
    if made.empty:
        failures.append(f"{section} draws: none written")
    for frame, name in ((made, "made"), (want, "archived")):
        if frame.duplicated(keys).any():
            failures.append(f"{section} draws: duplicated keys in the {name} draws")
    if failures:
        return failures, 0, []
    both = made.merge(want, on=keys, how="left", suffixes=("", "__reviewed"), indicator=True)
    unmatched = both[both["_merge"].ne("both")]
    if len(unmatched):
        failures.append(f"{section} draws: {len(unmatched)} made draws are not in the archive, "
                        f"e.g. {unmatched[keys].iloc[0].to_dict()}")
    both = both[both["_merge"].eq("both")]
    # Every archived quantity is compared, including ones the archive leaves missing for these
    # draws: a number where the reviewed analysis has none is a failure (cf. R120-F12).
    for column in [c for c in want.columns if c not in keys]:
        if column not in made.columns:
            failures.append(f"{section} draws: archived quantity {column} was not produced")
            continue
        a, b = both[column].to_numpy(float), both[column + "__reviewed"].to_numpy(float)
        tolerance = np.full(len(both), DRAW_TOLERANCE)
        if not tolerances.empty:
            for t in tolerances[tolerances.section == section].itertuples():
                if any(column.startswith(p.rstrip("*")) for p in t.columns.split(";")):
                    rows = both["era"].eq(t.era).to_numpy() if "era" in both else np.ones(len(both), bool)
                    tolerance[rows] = float(t.tolerance)
        same_missing = np.isnan(a) & np.isnan(b)
        gap = np.where(same_missing, 0.0, np.abs(a - b) / np.maximum(1.0, np.abs(b)))
        if not (gap <= tolerance).all():               # a missing value on one side only also fails here
            failures.append(f"{section} draws: {column} differs by up to {np.nanmax(gap):.2g}"
                            + (" (missing on one side)" if np.isnan(gap).any() else ""))
    unverified = [c for c in made.columns if c not in keys and c not in want.columns]
    return failures, len(both), unverified


# ---------------------------------------------------------------------------
def main() -> int:
    sys.stdout.reconfigure(encoding="utf-8")       # the paper's minus sign, on any console
    parser = argparse.ArgumentParser(description="Check the pipeline's outputs against the reviewed analysis.")
    parser.add_argument("sections", nargs="*", metavar="section",
                        help="sections to check, s1 ... s7 (default: every section that has outputs)")
    parser.add_argument("--all", action="store_true", help="check all seven sections; a missing one fails")
    args = parser.parse_args()
    unknown = [s for s in args.sections if s not in SECTIONS]
    if unknown:
        parser.error(f"unknown section(s): {', '.join(unknown)}")

    expected = pd.read_csv(EXPECTED / "expected_values.csv", keep_default_na=False,
                           na_values=[""], dtype={"printed": str, "printed_low": str, "printed_high": str})
    exceptions, tolerances = _optional_csv("table_format_exceptions.csv"), _optional_csv("draw_tolerances.csv")
    receipts, sources, receipt_failures = load_receipts()
    for line in receipt_failures:
        print("FAIL: " + line)
    if args.all:
        sections = list(SECTIONS)
    else:
        sections = args.sections or [s for s in SECTIONS if (S.OUTPUTS / s / "numbers.csv").is_file()]

    total_failures, corrections, formatting_all = len(receipt_failures), [], []
    for section in sections:
        rows = expected[expected.section.eq(section)]
        folder = S.OUTPUTS / section
        if not (folder / "numbers.csv").is_file():
            print(f"{section}: FAIL (no outputs; run the section script)")
            total_failures += 1
            continue
        try:
            record = json.loads((folder / "run.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            print(f"{section}: FAIL (no run record outputs/{section}/run.json; rerun the section script)")
            total_failures += 1
            continue
        receipt_failures, verified, stale = check_receipts(section, receipts, sources, rows, record)
        failures, counts = check_numbers(section, rows, corrections, verified)
        table_failures, formatting, table_pending = check_tables(section, exceptions)
        draw_failures, draws, unverified = check_draws(section, record, tolerances)
        failures += receipt_failures + table_failures + draw_failures
        formatting_all += formatting
        total_failures += len(failures)
        verdict = "PASS" if not failures else f"FAIL ({len(failures)})"
        notes = []
        if counts.get("partial"):
            notes.append(f"{counts['partial']} numbers not made by this partial run")
        if counts.get("verified_on_archive"):
            notes.append(f"{counts['verified_on_archive']} intervals or summaries verified on the archived results")
        if counts.get("pending"):
            notes.append(f"{counts['pending']} intervals pending a full run")
        if table_pending:
            notes.append(f"{table_pending} table cells' intervals pending a full run")
        if stale:
            notes.append(stale)
        if draws:
            notes.append(f"{draws} draws or evaluations match the archive")
        if unverified:
            notes.append(f"{len(unverified)} drawn quantities have no archived draws to compare")
        extra = f" [{'; '.join(notes)}]" if notes else ""
        print(f"{section}: {len(rows)} numbers, {len(TABLES.get(section, []))} table(s): {verdict}{extra}")
        for line in failures:
            print("   " + line)
    if formatting_all:
        print("\nAllowed formatting differences (expected/table_format_exceptions.csv):")
        for line in formatting_all:
            print("   " + line)
    if corrections:
        print("\nKnown errors in the manuscript (the pipeline agrees with the reviewed analysis):")
        for line in corrections:
            print("   " + line)
    return 1 if total_failures else 0


if __name__ == "__main__":
    sys.exit(main())
