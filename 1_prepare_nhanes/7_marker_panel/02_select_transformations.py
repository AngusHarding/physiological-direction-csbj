# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/7_marker_panel/02_select_transformations.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Chooses one shape transformation per marker from the five
#               candidates, using only the markers' distributions in the
#               29,053 complete-panel adults (women and men together), and
#               applies it. No outcome, age, sex difference or reference
#               membership is used.
#
# METHODS       "Panel and transformations": "one transformation per marker
#               (identity, square root, cube root, logarithm or reciprocal) was
#               chosen by a fixed rule on the pooled-sex adult panel, without
#               reference to outcomes." Supp. Section Zero, "Transformations".
#
# PSEUDOCODE
#   check the input is the certified raw panel (29,053 rows; SEQN, sex, age
#     and the 38 markers only, so no outcome column can be present)
#   for each marker:
#     if the raw values pass all three tests -> keep identity
#     else -> the first of sqrt, cbrt, log, reciprocal that passes all three
#     if none passes all three -> the one passing most; ties by rank sum,
#       then by that order
#   tests: |skew| <= 1, |excess kurtosis| <= 3, normal Q-Q R^2 >= 0.98
#   stability (descriptive only, never changes the choice): rerun within each
#     sex and in 100 sex-stratified bootstrap samples
# =============================================================================
"""Select and apply transformations using only the locked full cohort.

Beginner's map
--------------
For every biomarker, this script tries the five formulas declared in
TRANSFORMATION_PROTOCOL.md. It describes how close each result is to a
well-behaved, symmetric distribution, ranks the candidates, and applies the
winner. Sex and bootstrap reruns are diagnostics only: they never rewrite the
full-cohort choice.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
RAW_PANEL = ROOT / "outputs/05_panel/final_38_biomarker_panel_raw.csv"
SPEC = ROOT / "config/final_38_panel_spec.csv"
OUT = ROOT / "outputs/06_transform"
PROVENANCE = ROOT / "outputs/provenance"  # public_code: records go with the outputs

EXPECTED_SHA256 = "D9B2E795379571BFF781F2EE6503FA81B6AF5C04E59C47103519A0C25DC7F07A"
EXPECTED_ROWS = 29053
CANDIDATES = ("identity", "cbrt", "sqrt", "log", "reciprocal")
# Once identity fails its own gate, prefer the least aggressive adequate option.
PARSIMONY = ("sqrt", "cbrt", "log", "reciprocal")
ORIENTATION = {"identity": 1, "cbrt": 1, "sqrt": 1, "log": 1, "reciprocal": -1}
METRICS = ("abs_skew", "abs_excess_kurtosis", "qq_deviation", "tail_asymmetry")
ABS_SKEW_LIMIT = 1.0
ABS_EXCESS_KURTOSIS_LIMIT = 3.0
QQ_R2_MINIMUM = 0.98
BOOTSTRAPS = 100
SEED = 20260801


def file_sha256(path: Path) -> str:
    """Return a byte-for-byte fingerprint without changing the file."""
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def transform_values(values: np.ndarray, candidate: str) -> np.ndarray | None:
    """Apply one candidate; return None when its mathematical domain fails."""
    x = np.asarray(values, dtype=float)
    if not np.isfinite(x).all():
        return None
    if candidate == "identity":
        return x.copy()
    if candidate == "cbrt":
        return np.cbrt(x)
    if candidate == "sqrt":
        return np.sqrt(x) if np.all(x >= 0) else None
    if candidate == "log":
        return np.log(x) if np.all(x > 0) else None
    if candidate == "reciprocal":
        return 1.0 / x if np.all(x > 0) else None
    raise ValueError(f"unknown transformation: {candidate}")


# -- CALCULATION: the distribution tests -------------------------------------
# Supp. Section Zero: "Identity was kept when absolute bias-corrected skewness
#   was at most 1, absolute excess kurtosis at most 3 and normal Q-Q R^2 at
#   least 0.98." Skewness and kurtosis are scipy's bias-corrected estimates;
#   Q-Q R^2 is the squared correlation of the sorted values with normal
#   quantiles (scipy.stats.probplot).
def diagnostics(values: np.ndarray, candidate: str) -> dict:
    """Calculate the four predeclared distribution diagnostics."""
    transformed = transform_values(values, candidate)
    base = {"candidate": candidate, "valid": transformed is not None,
            "orientation": ORIENTATION[candidate]}
    if transformed is None or not np.isfinite(transformed).all() or np.ptp(transformed) == 0:
        return base | {m: np.inf for m in METRICS}

    q01, q50, q99 = np.quantile(transformed, [0.01, 0.50, 0.99])
    lower_gap, upper_gap = q50 - q01, q99 - q50
    tail = (abs(np.log(upper_gap / lower_gap))
            if lower_gap > 0 and upper_gap > 0 else np.inf)
    # scipy.probplot returns the expected normal order statistics and sorted data.
    theoretical, ordered = stats.probplot(transformed, dist="norm", fit=False)
    qq_r = np.corrcoef(theoretical, ordered)[0, 1]
    return base | {
        "abs_skew": abs(float(stats.skew(transformed, bias=False))),
        "abs_excess_kurtosis": abs(float(stats.kurtosis(transformed, fisher=True, bias=False))),
        "qq_deviation": float(1.0 - qq_r**2),
        "tail_asymmetry": float(tail),
        "transformed_min": float(np.min(transformed)),
        "transformed_max": float(np.max(transformed)),
    }


# -- CALCULATION: the decision rule ------------------------------------------
# Supp. Section Zero: "Otherwise the first transformation meeting all three
#   criteria was chosen, in the order square root, cube root, logarithm,
#   reciprocal; if none met all three, the one meeting most was chosen, with
#   ties broken by the rank sum of the three statistics and then by the same
#   order. No offsets were added, so a genuine zero ruled out the logarithm
#   and reciprocal."
def evaluate_marker(values: np.ndarray, marker: str, retain_identity: bool = False) -> pd.DataFrame:
    """Build the candidate ledger and apply the locked two-step decision rule."""
    rows = []
    for candidate in CANDIDATES:
        row = diagnostics(values, candidate)
        row.update(marker=marker, raw_min=float(np.min(values)),
                   raw_max=float(np.max(values)), zeros=int(np.sum(values == 0)),
                   unique_values=int(np.unique(values).size))
        rows.append(row)
    ledger = pd.DataFrame(rows)
    valid = ledger.valid
    ledger["skew_pass"] = ledger.abs_skew.le(ABS_SKEW_LIMIT)
    ledger["kurtosis_pass"] = ledger.abs_excess_kurtosis.le(ABS_EXCESS_KURTOSIS_LIMIT)
    ledger["qq_r2"] = 1.0 - ledger.qq_deviation
    ledger["qq_pass"] = ledger.qq_r2.ge(QQ_R2_MINIMUM)
    ledger["criteria_passed"] = ledger[["skew_pass", "kurtosis_pass", "qq_pass"]].sum(axis=1)
    ledger["adequacy"] = np.select(
        [~valid, ledger.criteria_passed.eq(3), ledger.criteria_passed.eq(2)],
        ["invalid", "acceptable", "borderline"], default="inadequate")
    # Ranks are retained for a deterministic best-available decision when no
    # nonlinear candidate reaches practical adequacy.
    for metric in ("abs_skew", "abs_excess_kurtosis", "qq_deviation"):
        ledger[f"rank_{metric}"] = np.nan
        ledger.loc[valid, f"rank_{metric}"] = ledger.loc[valid, metric].rank(method="average")
    ledger["rank_sum"] = ledger[["rank_abs_skew", "rank_abs_excess_kurtosis",
                                  "rank_qq_deviation"]].sum(axis=1, min_count=3)
    ledger["parsimony_rank"] = ledger.candidate.map(
        {"identity": 0, **{c: i + 1 for i, c in enumerate(PARSIMONY)}})
    ledger["selected"] = False
    if retain_identity:
        winner = ledger.index[ledger.candidate.eq("identity")][0]
        basis = "raw_identity_passed_all_three"
    else:
        nonlinear = ledger.loc[valid & ledger.candidate.ne("identity")]
        acceptable = nonlinear.loc[nonlinear.adequacy.eq("acceptable")]
        if not acceptable.empty:
            winner = acceptable.sort_values("parsimony_rank").index[0]
            basis = "least_aggressive_adequate"
        else:
            best_pass_count = nonlinear.criteria_passed.max()
            winner = nonlinear.loc[nonlinear.criteria_passed.eq(best_pass_count)].sort_values(
                ["rank_sum", "parsimony_rank"]).index[0]
            basis = "best_available_not_adequate"
    ledger.loc[winner, "selected"] = True
    ledger["selection_basis"] = basis
    return ledger


def choose(values: np.ndarray, marker: str, retain_identity: bool) -> str:
    """Return the selected candidate, used by stability reruns."""
    ledger = evaluate_marker(values, marker, retain_identity=retain_identity)
    return str(ledger.loc[ledger.selected, "candidate"].iloc[0])


# -- CHECK: outcome-blind input ----------------------------------------------
# Refuses any file except the certified raw panel: exact fingerprint, 29,053
# rows, and columns exactly SEQN, sex, age and the 38 markers. No mortality
# file is read anywhere in this script.
def validate_input(panel: pd.DataFrame, spec: pd.DataFrame) -> list[str]:
    """Refuse to analyse anything except the certified raw panel."""
    markers = spec.marker.tolist()
    expected_columns = ["SEQN", "RIAGENDR", "RIDAGEYR"] + markers
    if file_sha256(RAW_PANEL) != EXPECTED_SHA256:
        raise AssertionError("raw-panel SHA-256 does not match the certified input")
    if panel.shape != (EXPECTED_ROWS, 41) or list(panel.columns) != expected_columns:
        raise AssertionError("raw-panel dimensions or ordered columns changed")
    if panel.SEQN.duplicated().any() or panel.isna().any().any():
        raise AssertionError("raw panel must have unique SEQN and no missing cells")
    if set(panel.RIAGENDR.unique()) != {1.0, 2.0}:
        raise AssertionError("unexpected sex codes")
    return markers


def main() -> None:
    panel = pd.read_csv(RAW_PANEL)
    spec = pd.read_csv(SPEC)
    markers = validate_input(panel, spec)

    # Gate identity once on the locked full cohort. Resampling diagnostics do
    # not reopen that decision.
    retain_identity = {}
    for marker in markers:
        raw = diagnostics(panel[marker].to_numpy(), "identity")
        retain_identity[marker] = (raw["abs_skew"] <= ABS_SKEW_LIMIT
                                   and raw["abs_excess_kurtosis"] <= ABS_EXCESS_KURTOSIS_LIMIT
                                   and (1.0 - raw["qq_deviation"]) >= QQ_R2_MINIMUM)
    ledgers = [evaluate_marker(panel[m].to_numpy(), m, retain_identity[m]) for m in markers]
    ledger = pd.concat(ledgers, ignore_index=True)
    selected = ledger.loc[ledger.selected, ["marker", "candidate", "orientation", "selection_basis",
                                             "adequacy", "criteria_passed", "rank_sum", *METRICS,
                                             "qq_r2", "zeros", "unique_values"]].copy()
    selected = selected.rename(columns={"candidate": "selected_transform"})

    # Sex-specific winners are diagnostic. The primary winners above remain fixed.
    sex_rows = []
    for marker in markers:
        row = {"marker": marker}
        for sex, label in ((1.0, "male"), (2.0, "female")):
            row[f"{label}_winner"] = choose(panel.loc[panel.RIAGENDR.eq(sex), marker].to_numpy(),
                                              marker, retain_identity[marker])
        sex_rows.append(row)
    sex_stability = pd.DataFrame(sex_rows).merge(
        selected[["marker", "selected_transform"]], on="marker", validate="one_to_one")
    sex_stability["sex_disagreement"] = (
        sex_stability.male_winner.ne(sex_stability.female_winner))

    # Resample men and women separately, preserving their observed sample sizes.
    rng = np.random.default_rng(SEED)
    by_sex = {sex: np.flatnonzero(panel.RIAGENDR.to_numpy() == sex) for sex in (1.0, 2.0)}
    bootstrap_rows = []
    for replicate in range(1, BOOTSTRAPS + 1):
        sampled = np.concatenate([rng.choice(idx, size=len(idx), replace=True) for idx in by_sex.values()])
        for marker in markers:
            bootstrap_rows.append({"replicate": replicate, "marker": marker,
                                   "winner": choose(panel[marker].to_numpy()[sampled], marker,
                                                    retain_identity[marker])})
    bootstrap_long = pd.DataFrame(bootstrap_rows)
    bootstrap = (bootstrap_long.groupby(["marker", "winner"]).size().rename("wins")
                 .reset_index())
    bootstrap["selection_frequency"] = bootstrap.wins / BOOTSTRAPS
    primary_frequency = selected[["marker", "selected_transform"]].merge(
        bootstrap, left_on=["marker", "selected_transform"], right_on=["marker", "winner"],
        how="left", validate="one_to_one")
    primary_frequency["selection_frequency"] = primary_frequency.selection_frequency.fillna(0.0)
    selected = selected.merge(primary_frequency[["marker", "selection_frequency"]],
                              on="marker", validate="one_to_one")
    selected = selected.merge(sex_stability[["marker", "sex_disagreement"]],
                              on="marker", validate="one_to_one")
    selected["stability_flag"] = selected.sex_disagreement | selected.selection_frequency.lt(0.70)

    transformed = panel[["SEQN", "RIAGENDR", "RIDAGEYR"]].copy()
    choice = selected.set_index("marker").selected_transform
    for marker in markers:
        transformed[marker] = transform_values(panel[marker].to_numpy(), choice[marker])
    if transformed.isna().any().any() or not np.isfinite(transformed[markers].to_numpy()).all():
        raise AssertionError("selected transformations created missing or infinite values")

    historical = spec[["marker", "transform", "shift"]].rename(
        columns={"transform": "frozen_manuscript_transform", "shift": "frozen_manuscript_shift"})
    comparison = selected[["marker", "selected_transform"]].merge(
        historical, on="marker", validate="one_to_one")
    comparison["same_transform"] = comparison.selected_transform.eq(
        comparison.frozen_manuscript_transform)

    OUT.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(OUT / "transformation_candidate_ledger_v2.csv", index=False)
    selected.to_csv(OUT / "transformation_selection_v2.csv", index=False)
    # public_code: config/ is read only, so check the stored recipe instead of rewriting it
    recipe = selected[["marker", "selected_transform", "orientation", "selection_basis", "adequacy"]]
    stored = (ROOT / "config/final_38_transformations_v2.csv").read_text(encoding="utf-8")
    if recipe.to_csv(index=False, lineterminator="\n") != stored.replace("\r\n", "\n"):
        raise ValueError("config/final_38_transformations_v2.csv differs from the transformations "
                         "selected here; script 05 reads that file, so check which is right")
    sex_stability.to_csv(OUT / "transformation_sex_stability_v2.csv", index=False)
    bootstrap.to_csv(OUT / "transformation_bootstrap_stability_v2.csv", index=False)
    transformed.to_csv(OUT / "final_38_panel_transformed_unstandardised_v2.csv", index=False)
    comparison.to_csv(PROVENANCE / "transformation_manuscript_comparison_v2.csv", index=False)
    print(f"selected and applied transformations for {len(markers)} biomarkers")
    print(f"stable choices: {(~selected.stability_flag).sum()}/38")
    print(f"agreement with frozen manuscript: {comparison.same_transform.sum()}/38")


if __name__ == "__main__":
    main()
