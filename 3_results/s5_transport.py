# =============================================================================
# READ ME FIRST: s5_transport.py  (Results 5)
# -----------------------------------------------------------------------------
# RESULTS 5     "A frozen q score ranked deaths across survey eras."
#
# WHAT IT DOES  Learns the whole 32-marker scoring pipeline in one survey era,
#               freezes it, and applies it unchanged in the other era, where
#               only the score's association with death is estimated. Done in
#               both directions (forward: 1999-2018 -> III; reverse: III ->
#               1999-2018).
#
# METHODS       Supplementary Section Five, "Evaluation"; "Held-out
#               prediction" (concordance); "Uncertainty" (survey bootstrap).
#
# MAKES         Figure 5b,c numbers; Supp. Tables S5.2-S5.4; random
#               directions (Supp. Fig. S5.1); every Results 5 number.
#
# PSEUDOCODE
#   for each direction (source -> target):
#     SOURCE ONLY: transformation recipe (stored, per era);
#       young-healthy mean and SD by sex from source reference members;
#       direction u from Cox(female + 32 coordinates), strata = cycle or phase;
#       scale = SD of q in the source cohort
#     TARGET: coordinates from the SOURCE reference; q = (z . u) / source SD
#       association: Cox(female + q) -> log HR per source SD
#       ranking: concordance within sex and survey period
#     comparisons: target's own reference; a q trained and cross-fitted in the
#       target; 1,000 random directions oriented and scaled in the source;
#       post hoc: omit influential cases; a joint transformation recipe
#   intervals: 2,000 Rao-Wu survey-bootstrap draws in each era
# =============================================================================
"""Results 5: a frozen score ranked deaths across survey eras.

The question
------------
Does the score keep its prognostic ranking when every rule is learned in one
survey era and applied, unchanged, in the other?

How
---
Using the 32 markers measured comparably in both eras (inputs/common32_markers.csv),
the whole scoring pipeline is learned in the source era:
  1. each marker's shape transformation (the source era's stored recipe);
  2. each sex's young-healthy mean and SD (source reference members);
  3. the mortality direction (a Cox model of death on female + the 32
     coordinates, stratified by survey period);
  4. the score's scale (its SD in the source cohort).
It is then frozen and applied in the target era, where only the score's
association with death is estimated: a Cox model on female + score, and
concordance among participants of the same sex and survey period. No
baseline hazard is carried over, so this tests ranking, not calibrated risk.

  forward  NHANES 1999-2018 -> NHANES III (4,439 participants, 3,374 deaths)
  reverse  NHANES III -> NHANES 1999-2018 (11,497 participants, 2,494 deaths)

Comparisons: a score trained and cross-fitted within the target; the target's
own reference in place of the source's; 1,000 random scoring directions
(oriented and scaled in the source only); post hoc, omitting three target
participants with ALT recorded at 1 U/L, and a joint transformation recipe
informed by both eras' marker distributions (no outcomes).

Uncertainty: 2,000 Rao-Wu survey-bootstrap draws in each era (numbered from
0), as in the reviewed runs (E8b, E8c). Two kinds of interval:
  joint        the source reference, direction and scale are re-learned with
               the source draw, and the target association re-estimated with
               the target draw (frozen, own-reference and joint-recipe scores);
  conditional  the target-trained score is held at its point value and only its
               target association is re-estimated with the target draw.
Transformation choices are fixed throughout, so neither interval covers them.

Scope
-----
Population   both analysis cohorts on the 32 common markers (NHANES III rows with
             complete common markers: 4,439 participants).
Endpoint     accident-censored mortality.
Estimands    in the target era, the log hazard ratio per source-score SD (with female
             in the model) and concordance within sex and survey period.
Split        the source learns everything; the target only evaluates. The
             target-trained comparator is cross-fitted within the target.
Assumptions  each era's stored transformation recipe is fixed; point estimates do
             not use NHANES sample weights (the survey bootstrap respects the design
             through its resampling of sampling units).
Permitted    whether a score learned in one era ranks deaths in the other (ranking,
interpretation  not calibrated risk). The omitted-case and joint-recipe analyses were
             chosen after the primary results were seen; they qualify the primary
             result and do not replace it.

What it makes (in outputs/s5/)
------------------------------
  transport_estimates.csv  every point estimate, with intervals once the bootstrap has run
  random_directions.csv    the 1,000 random directions' target log hazard ratios and concordance
  check_draws.csv          survey-bootstrap draws for the draw-by-draw check (--check-draws N)
  numbers.csv              every Results 5 number quoted in the text

    python -B s5_transport.py                  # point estimates and random directions (minutes)
    python -B s5_transport.py --check-draws 3  # three survey-bootstrap draws, to compare with the archive
    python -B s5_transport.py --bootstrap      # all 2,000 survey-bootstrap draws (hours)
"""
from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
from sklearn.model_selection import GroupKFold   # noqa: E402

import data                               # noqa: E402
import models                             # noqa: E402
import runs                               # noqa: E402
import settings as S                      # noqa: E402

OUT = S.output_folder("s5")
ERAS = ("modern", "iii")
TESTS = {"forward": ("modern", "iii"), "reverse": ("iii", "modern")}     # (source, target)
RECIPE = {"modern": "modern_recipe", "iii": "nhanes3_recipe"}           # the source era's stored recipe
RANDOM_SEED, RANDOM_DIRECTIONS = 20261012, 1000
DRAWS = 2000


# ===========================================================================
# Data
# ===========================================================================
# -- COHORTS: both eras on the 32 common markers -----------------------------
# Results 5: "NHANES 1999-2018 to NHANES III (forward; 4,439 target
#   participants, 3,374 deaths) and NHANES III to NHANES 1999-2018 (reverse;
#   11,497 participants, 2,494 deaths)". The NHANES III cohort is smaller than
#   in Results 1-4 (4,474) because it needs the 32 common markers rather than
#   NHANES III's own 34.
def load() -> dict:
    """Both eras: the full panel (for references and survey weights) and the analysed rows."""
    keys = data.common32_markers().key.tolist()
    out = {}
    for era in ERAS:
        frame = data.load_common32(era)
        if era == "modern":
            people = data.load_modern().people
            people = people.assign(period=people.cycle)
            reference = frame.is_reference.to_numpy()
        else:
            complete = frame[frame.complete]
            mort = data.mortality_iii()
            people = complete[["SEQN", "sex", "age", "phase", "psu"]].merge(mort, on="SEQN", how="inner")
            people = people[people.age.ge(S.AGE_MIN) & people.permth_exm.gt(0)].reset_index(drop=True)
            people["female"] = people.sex.eq(2).astype(int)
            data.add_follow_up(people, "age")
            people["period"] = people.phase.astype(str)
            reference = (frame.is_reference & frame.complete).to_numpy()
        position = pd.Series(np.arange(len(frame)), index=frame.SEQN).loc[people.SEQN].to_numpy()
        out[era] = {"frame": frame, "people": people.reset_index(drop=True), "position": position,
                    "reference": reference, "sex": frame.sex.to_numpy(), "values": frame[keys].to_numpy(float),
                    "psu": frame.psu.to_numpy()[position]}
    assert (len(out["iii"]["people"]), out["iii"]["people"].event.sum()) == (4_439, 3_374)
    return out


def transformed(d: dict, recipe: str) -> np.ndarray:
    kinds = data.common32_markers()[recipe].tolist()
    return np.column_stack([data.transform(d["values"][:, j], k) for j, k in enumerate(kinds)])


# -- CALCULATION: coordinates from a reference -------------------------------
# (transformed value - reference mean) / reference SD, separately by sex: the
# same rule as the SS reference in Results 1-2. With survey-bootstrap weights,
# the mean and SD are weighted.
def scored(T_rows, sex_rows, T_ref, sex_ref, weights=None) -> np.ndarray:
    """(transformed value - reference mean) / reference SD, by sex. Weighted in bootstrap draws."""
    w = np.ones(len(T_ref)) if weights is None else weights
    out = np.empty_like(T_rows)
    for s in (1, 2):
        keep = (sex_ref == s) & (w > 0)
        ws, raw = w[keep], T_ref[keep]
        n = ws.sum()
        mu = (ws[:, None] * raw).sum(axis=0) / n
        sd = np.sqrt((ws[:, None] * (raw - mu) ** 2).sum(axis=0) / (n - 1))
        out[sex_rows == s] = (T_rows[sex_rows == s] - mu) / sd
    return out


# Cox model on attained age with entry, Efron ties and period strata (cycle in
# NHANES 1999-2018, phase in NHANES III, as Supp. Section Five states).
def cox_fit(people: pd.DataFrame, X: np.ndarray, weights=None, rows=None) -> np.ndarray:
    """Cox coefficients (attained age, Efron ties, survey-period strata), optionally weighted."""
    keep = np.ones(len(people), bool) if rows is None else rows
    if weights is not None:
        keep = keep & (weights > 0)
    p = people[keep]
    risk = models.RiskSets(p.entry_age.to_numpy(), p.exit_age.to_numpy(), p.event.to_numpy(),
                           p.period.to_numpy(str), None if weights is None else weights[keep])
    beta, ok = models.ridge_cox(risk, X[keep], np.zeros(X.shape[1]))
    if not ok:
        raise ArithmeticError("Cox fit did not converge")
    return beta


# -- CALCULATION: the target association -------------------------------------
# Results 5: "In the target we estimated only the q score's association with
#   mortality." "The hazard ratio per source SD was 1.25 forward ... and 1.76 in
#   reverse" [logHR, female in the model]. "Concordance, comparing participants
#   of the same sex and survey period, was 0.630 forward ... and 0.659 in
#   reverse" [C, models.concordance_within].
def evaluate(people: pd.DataFrame, score: np.ndarray, weights=None) -> dict:
    """Target association: log HR per source-score SD (with female in the model) and concordance."""
    beta = cox_fit(people, np.column_stack([people.female, score]), weights)
    c = models.concordance_within(people.entry_age.to_numpy(), people.exit_age.to_numpy(), people.event.to_numpy(),
                                  people.female.to_numpy(), people.period.to_numpy(str), score, weights)
    return {"logHR": float(beta[1]), "C": c}


# ===========================================================================
# One frozen-transport test
# ===========================================================================
# -- CALCULATION: everything learned in the source era -----------------------
# Results 5: "we fitted the complete scoring pipeline in one era: each marker's
#   transformation, the sex-specific young-healthy means and SDs, the Cox
#   mortality direction and the q score's scale. We then froze it and applied it
#   unchanged in the other era."
# How: reference mean and SD from source reference members; u = female-
#   adjusted Cox marker coefficients / their length (female left out of q);
#   sd = SD of q in the source cohort. No target data enter.
def learn_source(src: dict, T_src: np.ndarray, weights=None):
    """Steps 2-4 in the source: reference, direction, scale. Returns (coordinate scorer, u, sd)."""
    w_ref = None if weights is None else weights[src["reference"]]
    Zs = scored(T_src[src["position"]], src["sex"][src["position"]], T_src[src["reference"]],
                src["sex"][src["reference"]], w_ref)
    w_rows = None if weights is None else weights[src["position"]]
    beta = cox_fit(src["people"], np.column_stack([src["people"].female, Zs]), w_rows)
    u = beta[1:] / np.linalg.norm(beta[1:])
    q = Zs @ u
    if w_rows is None:
        sd = q.std(ddof=1)
    else:
        mu = np.dot(w_rows, q) / w_rows.sum()
        sd = np.sqrt(np.dot(w_rows, (q - mu) ** 2) / (w_rows.sum() - 1))
    return Zs, u, sd


# Target participants scored against the SOURCE reference (primary), or against
# the target's own reference (own=True): Results 5, "Ranking was similar when
# the target's own young-healthy reference replaced the source's (0.631 and
# 0.658)".
def target_coordinates(src, tgt, T_src, T_tgt, own=False, weights_src=None, weights_tgt=None):
    """Target coordinates with the source's reference (or, with own=True, the target's)."""
    ref, T_ref, w = (tgt, T_tgt, weights_tgt) if own else (src, T_src, weights_src)
    return scored(T_tgt[tgt["position"]], tgt["sex"][tgt["position"]], T_ref[ref["reference"]],
                  ref["sex"][ref["reference"]], None if w is None else w[ref["reference"]])


# -- COMPARATOR: a q score trained within the target -------------------------
# Results 5: "against 0.642 and 0.666 for q scores trained and cross-fitted
#   within each target era". Leave one cycle out in NHANES 1999-2018; ten folds
#   of sampling units in NHANES III.
def cross_fitted(tgt: dict, Z: np.ndarray) -> np.ndarray:
    """A score trained within the target: leave one cycle out (1999-2018) or ten sampling-unit folds (NHANES III)."""
    p = tgt["people"]
    q = np.zeros(len(Z))
    if p.period.nunique() > 2:
        splits = [(np.flatnonzero(p.period.ne(g)), np.flatnonzero(p.period.eq(g))) for g in sorted(p.period.unique())]
    else:
        splits = list(GroupKFold(10).split(Z, groups=tgt["psu"]))
    for train, test in splits:
        rows = np.zeros(len(Z), bool)
        rows[train] = True
        beta = cox_fit(p, np.column_stack([p.female, Z]), rows=rows)
        q[test] = Z[test] @ (beta[1:] / np.linalg.norm(beta[1:]))
    return q / q.std(ddof=1)


# -- CALCULATION: point estimates for one direction --------------------------
# Frozen score, own-reference score, target-trained score, and two post hoc
# checks: Results 5, "omitting them raised the hazard ratio to 1.45" (the three
# participants with ALT recorded at 1 U/L) and "A post hoc transformation
# recipe informed by both eras' marker distributions, without outcomes, gave
# 1.56". Random directions: "None of 1,000 random scoring directions, oriented
# and scaled using the source only, ranked deaths as well as the learned q
# score"; "Forward, 50 of the 1,000 random directions had a larger log hazard
# ratio". Each random direction's sign is set by its source Cox fit.
def point(d: dict, test: str, index: int) -> tuple[list, pd.DataFrame, np.ndarray]:
    """Point estimates for one test; also returns the target-trained score, held fixed in the draws."""
    s, t = TESTS[test]
    src, tgt = d[s], d[t]
    rows = []
    for label, recipe in (("frozen", RECIPE[s]), ("joint_recipe", "joint_recipe")):
        T_src, T_tgt = transformed(src, recipe), transformed(tgt, recipe)
        Zs, u, sd = learn_source(src, T_src)
        Zt = target_coordinates(src, tgt, T_src, T_tgt)
        rows.append({"test": test, "analysis": label, **evaluate(tgt["people"], Zt @ u / sd)})
        if label != "frozen":
            continue
        Zo = target_coordinates(src, tgt, T_src, T_tgt, own=True)
        rows.append({"test": test, "analysis": "own_reference", **evaluate(tgt["people"], Zo @ u / sd)})
        trained = cross_fitted(tgt, Zt)
        rows.append({"test": test, "analysis": "target_trained", **evaluate(tgt["people"], trained)})
        # Post hoc influence check (chosen after the primary results were seen). Forward: the
        # three highest frozen scores, all participants with ALT recorded at 1 U/L. Reverse: the
        # highest score in sampling unit 1999-2000|3|2, the most influential unit in the bootstrap.
        q = Zt @ u / sd
        if test == "forward":
            omit = np.argsort(-q)[:3]
            alt = tgt["values"][tgt["position"], data.common32_markers().key.tolist().index("alt")]
            assert (alt[omit] == 1).all(), "the three highest scores should be the ALT = 1 U/L participants"
        else:
            unit = np.flatnonzero(tgt["psu"] == "1999-2000|3|2")
            omit = unit[[np.argmax(q[unit])]]
        keep = np.ones(len(q), bool)
        keep[omit] = False
        rows.append({"test": test, "analysis": "omit_cases", "omitted": len(omit),
                     **evaluate(tgt["people"][keep].reset_index(drop=True), q[keep])})
        # 1,000 random directions, each oriented and scaled in the source only.
        rng = np.random.default_rng([RANDOM_SEED, index])
        U = rng.normal(size=(RANDOM_DIRECTIONS, Zs.shape[1]))
        U /= np.linalg.norm(U, axis=1)[:, None]
        scales = (Zs @ U.T).std(axis=0, ddof=1)
        random = []
        for j in range(RANDOM_DIRECTIONS):
            b = cox_fit(src["people"], np.column_stack([src["people"].female, Zs @ U[j] / scales[j]]))
            sign = 1 if b[1] >= 0 else -1
            random.append({"test": test, "direction": j, **evaluate(tgt["people"], sign * (Zt @ U[j]) / scales[j])})
        frozen = rows[0]
        random = pd.DataFrame(random)
        rows.append({"test": test, "analysis": "random_directions",
                     "logHR_above": int((random.logHR > frozen["logHR"]).sum()),
                     "C_above": int((random.C > frozen["C"]).sum())})
    return rows, random, trained


# ===========================================================================
# Survey bootstrap
# ===========================================================================
_W = {}


def _start(trained: dict) -> None:
    """Worker set-up: both eras, the replicate weights, and each test's fixed target-trained score."""
    d = load()
    weights = {}
    for stream, era in enumerate(ERAS):
        unit_weight, unit_of = models.survey_bootstrap_weights(d[era]["frame"], stream, DRAWS)
        weights[era] = (unit_weight, unit_of)
    _W.update(d=d, weights=weights, trained=trained)


# -- UNCERTAINTY: one survey-bootstrap draw ----------------------------------
# Results 5: "The intervals came from 2,000 survey-design bootstrap draws that
#   also refitted the source era's reference means and SDs, mortality direction
#   and score scale, with transformations held fixed." The target association
#   is re-estimated with the target's draw weights. The target-trained score is
#   held at its point value ("conditional"; Figure 5 legend: "for the
#   target-trained q score they hold its scoring fixed").
def one_draw(draw: int) -> dict:
    """One survey-bootstrap draw (0-based, as the archive numbers them).

    joint: the source reference, direction and scale are re-learned with the
    draw's weights (transformations fixed), and the target association is
    re-estimated with the target's weights.
    """
    d, weights = _W["d"], _W["weights"]
    w = {era: weights[era][0][draw, weights[era][1]] for era in ERAS}
    out = {"draw": draw}
    for k, (test, (s, t)) in enumerate(TESTS.items(), start=1):
        src, tgt = d[s], d[t]
        T_src, T_tgt = transformed(src, RECIPE[s]), transformed(tgt, RECIPE[s])
        _, u, sd = learn_source(src, T_src, w[s])
        wt = w[t][tgt["position"]]
        Zt = target_coordinates(src, tgt, T_src, T_tgt, weights_src=w[s])
        Zo = target_coordinates(src, tgt, T_src, T_tgt, own=True, weights_tgt=w[t])
        for mode, Z in (("primary", Zt), ("own_reference", Zo)):
            r = evaluate(tgt["people"], Z @ u / sd, wt)
            out[f"test{k}.{mode}.joint.logHR"], out[f"test{k}.{mode}.joint.C"] = r["logHR"], r["C"]
        # Conditional: the target-trained score at its point value, target association re-estimated.
        r = evaluate(tgt["people"], _W["trained"][test], wt)
        out[f"test{k}.target_trained.conditional.logHR"], out[f"test{k}.target_trained.conditional.C"] = r["logHR"], r["C"]
        # Post hoc joint recipe, same draws (archive: revision/notes/E8c_transport_posthoc)
        T_src, T_tgt = transformed(src, "joint_recipe"), transformed(tgt, "joint_recipe")
        _, u, sd = learn_source(src, T_src, w[s])
        r = evaluate(tgt["people"], target_coordinates(src, tgt, T_src, T_tgt, weights_src=w[s]) @ u / sd, wt)
        out[f"test{k}.joint_recipe.logHR"], out[f"test{k}.joint_recipe.C"] = r["logHR"], r["C"]
    return out


# Which drawn column gives each analysis's interval (joint, except the conditional target-trained score).
COLUMN = {"frozen": "primary.joint", "own_reference": "own_reference.joint", "joint_recipe": "joint_recipe",
          "target_trained": "target_trained.conditional"}


def attach_intervals(est: pd.DataFrame, boot: pd.DataFrame | None) -> pd.DataFrame:
    """Add 95% percentile intervals (numpy linear rule) of logHR and C to the point estimates."""
    est = est.copy()
    for quantity in ("logHR", "C"):
        est[f"{quantity}_low"] = est[f"{quantity}_high"] = np.nan
    if boot is None:
        return est
    for i, r in est.iterrows():
        k = 1 if r.test == "forward" else 2
        if r.analysis not in COLUMN:
            continue
        for quantity in ("logHR", "C"):
            low, high = models.percentile_interval(boot[f"test{k}.{COLUMN[r.analysis]}.{quantity}"])
            est.loc[i, f"{quantity}_low"], est.loc[i, f"{quantity}_high"] = low, high
    return est


def quoted_numbers(est: pd.DataFrame, d: dict) -> pd.DataFrame:
    """numbers.csv: every Results 5 number, hazard ratios exponentiated."""
    numbers = []
    for test, (s, t) in TESTS.items():
        p = d[t]["people"]
        numbers += [(f"{test}_participants", len(p)), (f"{test}_deaths", int(p.event.sum()))]
    for r in est.itertuples():
        if r.analysis == "random_directions":
            numbers += [(f"{r.test}_random_logHR_above", r.logHR_above), (f"{r.test}_random_C_above", r.C_above)]
            continue
        numbers += [(f"{r.test}_{r.analysis}_C", r.C, r.C_low, r.C_high),
                    (f"{r.test}_{r.analysis}_HR", np.exp(r.logHR), np.exp(r.logHR_low), np.exp(r.logHR_high))]
    numbers = [n if len(n) == 4 else (*n, np.nan, np.nan) for n in numbers]
    return pd.DataFrame(numbers, columns=["key", "value", "low", "high"]).assign(era="transport")


# ===========================================================================
# Figure 5
# ===========================================================================
def figure5(est: pd.DataFrame, d: dict) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    import displays as V

    fig = plt.figure(figsize=(185 * V.MM, 160 * V.MM), facecolor="white")
    a = fig.add_axes([.03, .69, .94, .25])
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.set_axis_off()
    boxes = [(.02, .38, "Learn in source era",
              "Transformations (32 markers)\nSex-specific young-reference means / SDs\nCox direction; source score SD", "#F0EEF7"),
             (.48, .2, "Apply in target era", "Apply frozen source objects\nwithout updating", "white"),
             (.76, .2, "Evaluate in target", "Score concordance\nRefit Cox association\nNo risk model transported", "white")]
    for x, w, title, body, fill in boxes:
        a.add_patch(FancyBboxPatch((x, .3), w, .6, boxstyle="round,pad=0.005,rounding_size=0.02",
                                   facecolor=fill, edgecolor=V.MUTED, lw=1))
        a.text(x + w / 2, .78, title, ha="center", fontsize=8, fontweight="bold")
        a.text(x + w / 2, .58, body, ha="center", va="center", fontsize=7, linespacing=1.5)
    for x0, x1 in ((.41, .47), (.69, .75)):
        a.annotate("", xy=(x1, .6), xytext=(x0, .6), arrowprops=dict(arrowstyle="-|>", color=V.PURPLE, lw=1.4))
    n_m, n_3 = d["modern"]["people"], d["iii"]["people"]
    a.text(.5, .1, f"Two separate source → target tests · NHANES 1999–2018: {len(n_m):,} participants / "
                   f"{int(n_m.event.sum()):,} deaths · NHANES III: {len(n_3):,} / {int(n_3.event.sum()):,} · ages 52–79",
           ha="center", fontsize=7, color=V.MUTED)

    e = est.set_index(["test", "analysis"])

    def draw(ax, y, value, low, high, colour, marker, filled):
        ax.hlines(y, low, high, color=colour, lw=1.5)
        ax.plot(value, y, marker=marker, ms=6, mew=1.3, mec=colour, mfc=colour if filled else "white")

    b = fig.add_axes([.22, .12, .3, .45])
    c = fig.add_axes([.72, .12, .26, .45])
    labels_b, labels_c = [], []
    for block, (test, title) in enumerate((("forward", "NHANES 1999–2018 → NHANES III"),
                                           ("reverse", "NHANES III → NHANES 1999–2018"))):
        top = 4.6 - block * 2.6
        r = e.loc[(test, "frozen")]
        draw(b, top, r.C, r.C_low, r.C_high, V.PURPLE, "o", True)
        r = e.loc[(test, "target_trained")]
        draw(b, top - 1, r.C, r.C_low, r.C_high, V.INK, "D", False)
        labels_b += [(top, "Frozen score"), (top - 1, "Trained in target")]
        b.text(.605, top + .55, title, fontsize=7, fontweight="bold")
        top = 5.0 - block * 3.0
        r = e.loc[(test, "frozen")]
        draw(c, top, np.exp(r.logHR), np.exp(r.logHR_low), np.exp(r.logHR_high), V.PURPLE, "o", True)
        r = e.loc[(test, "joint_recipe")]
        draw(c, top - .9, np.exp(r.logHR), np.exp(r.logHR_low), np.exp(r.logHR_high), V.PURPLE, "D", False)
        r = e.loc[(test, "omit_cases")]
        c.plot(np.exp(r.logHR), top - 1.8, marker="^", ms=6, mew=1.3, mec=V.PURPLE, mfc="white")
        omitted = "Three cases omitted" if test == "forward" else "One case omitted"
        labels_c += [(top, "Primary frozen score"), (top - .9, "Joint recipe\n(post hoc)"),
                     (top - 1.8, f"{omitted}\n(post hoc)")]
        short = title.replace("NHANES 1999–2018 →", "Modern →").replace("→ NHANES 1999–2018", "→ modern")
        c.text(2.05, top + .55, short, fontsize=7, fontweight="bold", ha="right")
    b.axhline(2.95, color=V.GRID, lw=1)
    c.axhline(2.5, color=V.GRID, lw=1)
    b.set(xlim=(.6, .7), ylim=(.4, 5.4), xlabel="Target concordance\n(within sex and survey period)")
    b.set_yticks([y for y, _ in labels_b])
    b.set_yticklabels([t for _, t in labels_b])
    c.set(xscale="log", xlim=(.95, 2.1), ylim=(.0, 5.8), xlabel="Hazard ratio per source SD\n(log scale)")
    c.xaxis.set_major_locator(FixedLocator([1, 1.25, 1.5, 2]))
    c.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    c.xaxis.set_minor_locator(NullLocator())
    c.axvline(1, color=V.MUTED, lw=.8, ls=(0, (3, 2)))
    c.set_yticks([y for y, _ in labels_c])
    c.set_yticklabels([t for _, t in labels_c])
    for ax in (b, c):
        ax.grid(axis="x", color=V.GRID, lw=.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0, pad=6)
    for key, x, y, title in (("a", .025, .965, "Source learning is frozen before target scoring"),
                             ("b", .025, .63, "Ranking in the target era"),
                             ("c", .53, .63, "Sensitivity of the hazard gradient")):
        fig.text(x, y, key, fontweight="bold", fontsize=10, va="center")
        fig.text(x + .033, y, title, fontsize=8.6, fontweight="bold", va="center")
    fig.text(.5, .015, "Bars: 95% intervals from 2,000 survey-design bootstrap draws (sampling units within strata); "
                       "trained-in-target bars hold that score fixed. Diamonds in c: post hoc joint recipe; "
                       "triangles: post hoc, point estimates only.",
             ha="center", fontsize=6.5, color=V.MUTED)
    V.save_figure(fig, OUT, "figure5")


# ===========================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check-draws", type=int, default=0, metavar="N")
    parser.add_argument("--bootstrap", action="store_true")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    run = runs.start("s5")

    d = load()
    rows, randoms, trained = [], [], {}
    for index, test in enumerate(TESTS):
        r, random, trained[test] = point(d, test, index)
        rows += r
        randoms.append(random)
        print(f"  {test}: point estimates and random directions done", flush=True)
    pd.concat(randoms).to_csv(OUT / "random_directions.csv", index=False)
    est = pd.DataFrame(rows)

    def survey_draws(n):
        with ProcessPoolExecutor(min(args.workers, n), initializer=_start, initargs=(trained,)) as pool:
            return pd.DataFrame(list(pool.map(one_draw, range(n))))

    stored, seed = OUT / "transport_bootstrap.csv", {"survey_weights": 20261011, "streams": list(ERAS)}
    fresh = None
    if args.check_draws:
        fresh = survey_draws(args.check_draws)
        run.write_check_draws(fresh, ["draw"])
    if args.bootstrap:
        boot = survey_draws(DRAWS)
        run.save_draws(boot, stored, ["draw"], seed)
    else:
        boot = run.stored_draws(stored, ["draw"], seed, runs.draw_keys(None, range(DRAWS), era_column=False), fresh)
    est = attach_intervals(est, boot)
    est.to_csv(OUT / "transport_estimates.csv", index=False)
    figure5(est, d)                                # interval bars appear once the bootstrap has run

    quoted_numbers(est, d).to_csv(OUT / "numbers.csv", index=False)
    run.finish(full_draws={"interval": DRAWS})
    print(f"done: outputs in {OUT}")


if __name__ == "__main__":
    main()
