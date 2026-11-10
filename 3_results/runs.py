# =============================================================================
# READ ME FIRST: runs.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Keeps a record of each run (code fingerprint, data
#               fingerprints, package versions, fit warnings) and decides
#               whether stored bootstrap draws may be reused: only if the data,
#               seed and code are unchanged. Comments and docstrings are not
#               part of the code fingerprint, so these READ ME blocks do not
#               invalidate stored draws.
#
# METHODS       None. Bookkeeping only; it calculates no estimate.
# =============================================================================
"""Run records and stored bootstrap draws: knowing which data and code made each output.

This file estimates nothing. It answers two practical questions.

1. What did the last run of a section do? Each section script calls
   `start(...)` first and `Run.finish()` last. Between them it writes
   outputs/<section>/run.json: the command, the time, the bootstrap draws it
   checked against the archive, a fingerprint (SHA-256) of every data file it
   read, and a fingerprint of the pipeline code. check_against_paper.py reads
   this record to know exactly which draws to compare; a check file left over
   from an earlier run is deleted by `start`, so it can never be checked by
   mistake.

2. May a stored bootstrap be reused? A full bootstrap takes minutes to hours,
   so a section reuses its stored draws (e.g. outputs/s1/distance_bootstrap.csv)
   when it can. Every stored file has a manifest (<file>.manifest.json) that
   says where its draws came from. Reuse needs, in every case:
     - every expected draw exactly once, no infinite value, no empty draw;
     - the file unchanged since its manifest was written (SHA-256);
     - the same seed, and the same data: every input in the manifest unchanged,
       and no input read by this run that the manifest does not list.
   Then it depends on the provenance the manifest records:
     - "computed": made by a full run of this pipeline. Reused only if the code
       fingerprint is the current one. If the code has changed, the script
       stops: matching a few check draws cannot show what changed code would
       give for the other draws, so stored draws are never relabelled as the
       current code's (R120-F4). Remaking them is a full run. Its manifest's
       input list is completed when that run finishes, so it names every file
       the making run read, not only those read before the draws were saved.
     - "archive_registered": draws shown, draw by draw and quantity by
       quantity, to equal a complete archived bootstrap of the reviewed
       analysis, by revision/notes/lean_pipeline_20261030/register_archive_draws.py.
       The manifest names the archive files and their SHA-256. run.json reports
       these as the archived bootstrap, not as draws made by the current code.
   When a run also makes check draws, they must agree with the stored draws
   (a consistency check, never a reason to reuse). Otherwise the script stops
   and says how to remake the draws; it never mixes draws from different data.

The code fingerprint hashes each pipeline module's syntax tree with docstrings
removed, so editing comments or docstrings does not invalidate stored draws,
while any change to what the code computes does.
"""
from __future__ import annotations

import ast
import hashlib
import json
import os
import platform
import sys
import time
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

import models
import settings as S

_HASH_MEMO = S.OUTPUTS / ".file_hashes.json"     # (size, modification time) -> hash, to avoid rehashing
PACKAGES = ("numpy", "pandas", "scipy", "statsmodels", "lifelines", "patsy", "matplotlib", "openpyxl",
            "scikit-learn")


# ---------------------------------------------------------------------------
# Fingerprints
# ---------------------------------------------------------------------------
def _relative(path: Path) -> str:
    for root, label in ((S.LOCAL_REPO, ""), (S.REPO, "NHANES_REPO:")):
        try:
            return label + path.relative_to(root).as_posix()
        except ValueError:
            continue
    return path.as_posix()


def file_hashes(paths) -> dict[str, str]:
    """SHA-256 of each file, reusing a stored hash when size and modification time are unchanged."""
    try:
        memo = json.loads(_HASH_MEMO.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        memo = {}
    out = {}
    for path in sorted(Path(p) for p in paths):
        stat = path.stat()
        key, stamp = str(path), [stat.st_size, stat.st_mtime_ns]
        if key in memo and memo[key][0] == stamp:
            digest = memo[key][1]
        else:
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for block in iter(lambda: fh.read(1 << 22), b""):
                    h.update(block)
            digest = h.hexdigest()
            memo[key] = [stamp, digest]
        out[_relative(path)] = digest
    _HASH_MEMO.write_text(json.dumps(memo, indent=0), encoding="utf-8")
    return out


def _resolve_relative(name: str) -> Path:
    if name.startswith("NHANES_REPO:"):
        return S.REPO / name[len("NHANES_REPO:"):]
    return S.LOCAL_REPO / name


def _syntax_hash(path: Path) -> str:
    tree = ast.parse(path.read_text(encoding="utf-8-sig"))
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and body \
                and isinstance(body[0], ast.Expr) and isinstance(getattr(body[0], "value", None), ast.Constant) \
                and isinstance(body[0].value.value, str):
            node.body = body[1:] or [ast.Pass()]
    return hashlib.sha256(ast.dump(tree, include_attributes=False).encode()).hexdigest()


def code_fingerprint() -> dict[str, str]:
    """Syntax-tree hash of every pipeline module loaded in this process, plus the shared renderer."""
    files = {Path(m.__file__).resolve() for m in list(sys.modules.values())
             if getattr(m, "__file__", None) and Path(m.__file__).resolve().parent == S.PIPELINE}
    if any(getattr(m, "__file__", None) and Path(m.__file__).resolve() == S.RENDERER
           for m in list(sys.modules.values())):
        files.add(S.RENDERER)
    return {_relative(f): _syntax_hash(f) for f in sorted(files)}


def _versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for name in PACKAGES:
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            out[name] = "not installed"
    return out


# ---------------------------------------------------------------------------
# The run record
# ---------------------------------------------------------------------------
class Run:
    def __init__(self, section: str):
        self.section = section
        self.folder = S.output_folder(section)
        self.record = {"section": section, "command": [Path(sys.argv[0]).name, *sys.argv[1:]],
                       "started": time.strftime("%Y-%m-%d %H:%M:%S"), "check_draws": None, "stored_draws": {}}
        self.made: list[Path] = []          # stored files made in this run; finish() completes their input lists
        for stale in ("run.json", "check_draws.csv", "check_draws_long.csv"):
            (self.folder / stale).unlink(missing_ok=True)
        # Fit warnings and failed convergence (models.note_fit) and files read (settings), from this
        # process and its workers; the environment variables are inherited by worker processes.
        self.diagnostics = self.folder / "process_logs"
        self.diagnostics.mkdir(exist_ok=True)
        for old in [*self.diagnostics.glob("*.jsonl"), *self.diagnostics.glob("reads_*.txt")]:
            old.unlink()
        os.environ[models.DIAGNOSTICS_ENV] = str(self.diagnostics)
        os.environ[S.READS_ENV] = str(self.diagnostics)

    # -- draws compared with the archive ------------------------------------------------------
    def write_check_draws(self, frame: pd.DataFrame, keys: list[str], name: str = "check_draws.csv",
                          index: bool = False) -> None:
        """Write the short-check draws and record which keys they cover (check_against_paper uses this)."""
        frame.to_csv(self.folder / name, index=index)
        written = pd.read_csv(self.folder / name)
        self.record["check_draws"] = {"file": name, "keys": keys, "rows": int(len(written))}

    # -- stored full bootstraps -------------------------------------------------------------
    def save_draws(self, frame: pd.DataFrame, path: Path, keys: list[str], seed) -> None:
        frame.to_csv(path, index=False)
        manifest = {"draws": path.name, "provenance": "computed", "file_sha256": _sha256(path), "keys": keys,
                    "rows": int(len(frame)), "seed": seed, "inputs": file_hashes(all_reads()),
                    "code": code_fingerprint(), "made": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "made_by": self.record["command"], "note": "made by a full bootstrap run"}
        _manifest(path).write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        self.record["stored_draws"][path.name] = "made in this run"
        self.made.append(path)

    def stored_draws(self, path: Path, keys: list[str], seed, expected_keys: pd.DataFrame,
                     fresh: pd.DataFrame | None = None, tolerance: float = 1e-6) -> pd.DataFrame | None:
        """The stored draws, if they can be vouched for (see the module docstring); None if absent.

        expected_keys: one row per draw the file must contain (columns = keys).
        fresh: this run's check draws, used to adopt draws made by other code.
        """
        if not path.is_file():
            self.record["stored_draws"][path.name] = "absent"
            return None
        frame = pd.read_csv(path)
        remake = f"Remake them with --bootstrap, or delete {path.name} to run without intervals."
        if frame.duplicated(keys).any():
            raise SystemExit(f"{path}: some draws appear twice. {remake}")
        have = set(map(tuple, frame[keys].astype(str).to_numpy()))
        want = set(map(tuple, expected_keys[keys].astype(str).to_numpy()))
        if have != want:
            raise SystemExit(f"{path}: has {len(have)} draws, expected {len(want)} "
                             f"({len(want - have)} missing, {len(have - want)} unexpected). {remake}")
        numeric = frame.drop(columns=keys).select_dtypes("number")
        if np.isinf(numeric.to_numpy(float)).any() or numeric.isna().all(axis=1).any():
            raise SystemExit(f"{path}: contains infinite values or empty draws. {remake}")

        manifest = _read_manifest(path)
        if manifest is None:
            raise SystemExit(f"{path}: has no manifest, so nothing shows which data and code made it. {remake}")
        self._same_origin(path, manifest, seed, remake)
        provenance = manifest.get("provenance")
        if provenance == "computed":
            if manifest.get("code") != code_fingerprint():
                raise SystemExit(
                    f"{path}: made by a different version of the code. Matching a few draws cannot show what "
                    f"this version gives for the rest, so the draws are not reused (R120-F4). {remake} "
                    f"(a full run). Draws equal to a complete archived bootstrap can instead be registered "
                    f"with revision/notes/lean_pipeline_20261030/register_archive_draws.py.")
            status = "reused: made by a full run of this code, same data and seed"
        elif provenance == "archive_registered":
            status = (f"archive-registered: equal to the complete archived bootstrap "
                      f"({manifest['archive']['summary']}); not recomputed by this code")
        else:
            raise SystemExit(f"{path}: its manifest records no known provenance ({provenance!r}). {remake}")
        if fresh is not None and not fresh.empty:
            worst = _largest_gap(fresh, frame, keys, path, remake)
            if not worst <= tolerance:
                raise SystemExit(f"{path}: this run's check draws differ from the stored draws "
                                 f"(largest scaled difference {worst:.2g}). {remake}")
            status += f"; this run's {len(fresh)} check draws agree with them (largest scaled difference {worst:.2g})"
        self.record["stored_draws"][path.name] = status
        return frame

    def _same_origin(self, path: Path, manifest: dict, seed, remake: str) -> None:
        """Stop unless the stored file, its seed and its data are those of this run."""
        if manifest.get("file_sha256") != _sha256(path):
            raise SystemExit(f"{path}: the file has changed since its manifest was written. {remake}")
        if manifest.get("seed") != seed:
            raise SystemExit(f"{path}: made with seed {manifest.get('seed')}, this run uses {seed}. {remake}")
        changed = [name for name, digest in manifest["inputs"].items()
                   if not _resolve_relative(name).is_file()
                   or file_hashes([_resolve_relative(name)]).get(name) != digest]
        new = sorted(set(file_hashes(all_reads())) - set(manifest["inputs"]))
        if changed or new:
            listed = (changed + new)[:3]
            raise SystemExit(f"{path}: made from different data ({len(changed)} input(s) changed, {len(new)} read "
                             f"now but not then: {', '.join(listed)}{' ...' if len(changed + new) > 3 else ''}). "
                             f"{remake}")

    # -- stored point results (recomputed rather than adopted when in doubt) -----------------
    def save_result(self, frame: pd.DataFrame, path: Path, seed) -> None:
        """Write a stored intermediate result (e.g. tuned penalties) with its manifest."""
        frame.to_csv(path, index=False)
        manifest = {"draws": path.name, "provenance": "computed", "file_sha256": _sha256(path),
                    "rows": int(len(frame)), "seed": seed, "inputs": file_hashes(all_reads()),
                    "code": code_fingerprint(), "made": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "made_by": self.record["command"], "note": "stored intermediate result"}
        _manifest(path).write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        self.record["stored_draws"][path.name] = "made in this run"
        self.made.append(path)

    def current_result(self, path: Path, seed) -> pd.DataFrame | None:
        """The stored result if its manifest matches this run's seed, data and code; otherwise None,
        and the caller recomputes it. Used where recomputing takes minutes, not hours."""
        manifest = _read_manifest(path)
        if not path.is_file() or manifest is None or manifest.get("provenance") != "computed" \
                or manifest.get("code") != code_fingerprint():
            self.record["stored_draws"][path.name] = "recomputed (absent, unrecorded or made by other code)"
            return None
        try:
            self._same_origin(path, manifest, seed, "")
        except SystemExit as why:
            self.record["stored_draws"][path.name] = f"recomputed ({why})"
            return None
        self.record["stored_draws"][path.name] = "reused (same data, seed and code)"
        return pd.read_csv(path)

    # -- the end of a run --------------------------------------------------------------------
    def fit_problems(self) -> dict:
        """Warnings and failed convergence recorded during the run, counted by kind and message."""
        entries = [json.loads(line) for f in sorted(self.diagnostics.glob("*.jsonl"))
                   for line in f.read_text(encoding="utf-8").splitlines() if line.strip()]
        counts: dict[str, int] = {}
        for e in entries:
            key = f"{e['kind']}: " + ("did not converge; " if not e["converged"] else "") + "; ".join(e["warnings"])
            counts[key] = counts.get(key, 0) + 1
        return {"fits_with_problems": len(entries), "by_message": counts, "first_examples": entries[:10]}

    def finish(self, **extra) -> None:
        self.record.update(extra)
        problems = self.fit_problems()
        self.record["fit_problems"] = problems
        if problems["fits_with_problems"]:
            print(f"  note: {problems['fits_with_problems']} model fits reported warnings or did not converge "
                  f"(listed in run.json under fit_problems)")
        inputs = file_hashes(all_reads())
        for path in self.made:                  # every file the run that made them read
            manifest = _read_manifest(path)
            manifest["inputs"] = inputs
            _manifest(path).write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        self.record.update(finished=time.strftime("%Y-%m-%d %H:%M:%S"), inputs=inputs,
                           code=code_fingerprint(), environment=_versions())
        (self.folder / "run.json").write_text(json.dumps(self.record, indent=1, default=str), encoding="utf-8")


def start(section: str) -> Run:
    return Run(section)


def all_reads() -> set[Path]:
    """Every data file read so far by this process and by the workers of this run (R120-F14)."""
    paths = set(S.READS)
    folder = os.environ.get(S.READS_ENV)
    if folder:
        for log in Path(folder).glob("reads_*.txt"):
            paths.update(Path(line) for line in log.read_text(encoding="utf-8").splitlines() if line)
    return {p for p in paths if p.is_file()}


def _sha256(path: Path) -> str:
    return next(iter(file_hashes([path]).values()))


def _largest_gap(fresh: pd.DataFrame, stored: pd.DataFrame, keys: list[str], path: Path, remake: str) -> float:
    """Largest scaled difference between this run's check draws and the same stored draws."""
    both = fresh.merge(stored, on=keys, how="left", suffixes=("", "__stored"), indicator=True)
    if not both["_merge"].eq("both").all():
        raise SystemExit(f"{path}: the check draws are not all in the stored file. {remake}")
    worst = 0.0
    for column in [c for c in fresh.columns if c not in keys]:
        if column not in stored.columns:
            raise SystemExit(f"{path}: has no column {column}. {remake}")
        a, b = both[column].to_numpy(float), both[column + "__stored"].to_numpy(float)
        gap = np.where(np.isnan(a) & np.isnan(b), 0.0, np.abs(a - b) / np.maximum(1.0, np.abs(b)))
        worst = max(worst, float(np.max(gap)) if np.isfinite(gap).all() else np.inf)
    return worst


def _manifest(path: Path) -> Path:
    return path.with_name(path.name + ".manifest.json")


def _read_manifest(path: Path) -> dict | None:
    try:
        return json.loads(_manifest(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def draw_keys(eras, draws, era_column: bool = True) -> pd.DataFrame:
    """The (era, draw) pairs a stored bootstrap must contain."""
    if not era_column:
        return pd.DataFrame({"draw": list(draws)})
    return pd.DataFrame([(e, d) for e in eras for d in draws], columns=["era", "draw"])
