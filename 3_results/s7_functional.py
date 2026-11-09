"""Results 7: lung function added prediction beyond the routine panel.

The question
------------
The routine panel contains no direct tests of function, such as lung function.
Does a measurement outside the panel (a direct functional test, or one of four
blood measures) improve held-out mortality prediction once the panel score is known?

How
---
For each candidate, in the participants who have it, three nested Cox models
are fitted on training folds and evaluated on the held-out fold:
  M0  baseline: sex, age, smoking (and height, for spirometry)
  M1  M0 + the panel mortality score q (learned without the held-out fold)
  M2  M1 + the candidate
Gains are held-out log partial likelihood differences summed over folds, per
1,000 evaluated deaths: the panel's gain (M1 - M0) and the candidate's (M2 - M1).
Folds: survey cycles in NHANES 1999-2018 (only cycles with at least 20
candidate deaths are evaluated; the others only train), ten folds of whole
sampling units in NHANES III. Main comparison, specified in advance: FEV1
beyond the full panel. A reduced baseline ("S2": sex, plus height for modern
spirometry) matches the original submission.

Candidates: high-sensitivity CRP, fasting glucose, fasting insulin, HOMA-IR and
FEV1 (quality grades A and B) in NHANES 1999-2018; FEV1, FVC, femur-neck bone
density, total bone mineral content, bioimpedance resistance and reactance and
tandem-stand time in NHANES III. Each is standardised (sign chosen so that
higher means worse), CRP within assay era.

Uncertainty: a bootstrap of 200 draws per era (numbered from 0): participants
within cycle x sex in NHANES 1999-2018, whole sampling units within strata in
NHANES III, each keeping their fold. Every fit, including the panel directions,
is repeated in each draw; 95% intervals are the 2.5th and 97.5th percentiles
(numpy's linear rule), left missing if more than 1% of draws fail, as in the
reviewed run (E5c, R038). The candidate measurements, their standardisation
and the locked ten markers are held fixed.

Scope
-----
Population   both analysis cohorts, restricted for each candidate to participants who
             have it (and height, for spirometry) and a known smoking status.
Endpoint     accident-censored mortality.
Estimands    the held-out gains of M1 over M0 (the panel) and M2 over M1 (the
             candidate), per 1,000 evaluated deaths, and the folds improved.
Split        directions and all three models are fitted on training folds; the
             evaluation uses the held-out fold, scaled with training statistics.
Assumptions  candidate standardisation and signs are fixed (higher is worse); the
             locked ten markers come from Results 4 (settings.LOCKED_TEN).
Permitted    whether a candidate adds held-out information beyond the panel in the
interpretation  participants who have it. Candidates are measured in different samples,
             so their gains are not ranked against each other; FEV1 beyond the full
             panel is the one comparison specified in advance.

What it makes (in outputs/s7/)
------------------------------
  gains.csv                 every candidate x baseline x panel: gains per 1,000 deaths, folds improved
  functional_bootstrap.csv  the 2 x 200 bootstrap draws (with --bootstrap) and its manifest
  intervals.csv             the 95% interval of both gains for every comparison
  figure7.svg / .png        Figure 7 (interval bars once the bootstrap exists)
  numbers.csv               every Results 7 number quoted in the text

    python -B s7_functional.py                  # point estimates (a few minutes); reuses stored draws
    python -B s7_functional.py --check-draws 1  # also draw 0 per era, compared with the archive
    python -B s7_functional.py --bootstrap      # all 2 x 200 draws
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

OUT = S.output_folder("s7")
CANDIDATE_DATA = S.EXTRACTS     # prepared candidate extracts
# The locked ten-marker panel of Table 3 (Results 4), in rank order.
TEN_MARKERS = S.LOCKED_TEN
SEED, DRAWS = 20261016, 200
SPIROMETRY = ("FEV1", "FVC")


def standard(x):
    return (x - x.mean()) / x.std(ddof=1)


# ===========================================================================
# Data
# ===========================================================================
def load_modern() -> dict:
    c = data.load_modern()
    d = c.people.copy()
    d[c.markers] = c.z.to_numpy()
    d["age"] = d.RIDAGEYR
    h = S.PREP / "02_harmonised"
    d = d.merge(pd.read_csv(h / "anthropometry.csv", usecols=["SEQN", "BMXHT"]), on="SEQN", how="left")
    d = d.merge(pd.read_csv(h / "smoking_status.csv", usecols=["SEQN", "Smoking_status"]), on="SEQN", how="left")
    d["height"], d["smoking"] = d.BMXHT, d.Smoking_status
    d["fold"] = pd.factorize(d.cycle, sort=True)[0]
    crp = pd.read_csv(CANDIDATE_DATA / "candidate_crp_raw.csv")
    crp = crp[crp.crp_raw.gt(0)].copy()
    crp["c_hsCRP"] = np.log(crp.crp_raw).groupby(crp.assay_era).transform(standard)   # within assay era
    glu = pd.read_csv(CANDIDATE_DATA / "candidate_glucose_raw.csv")
    glu = glu[glu.WTSAF2YR.notna()].copy()                                           # fasting subsample
    glu["c_glucose"] = standard(glu.LBXGLU.where(glu.LBXGLU.gt(0)))
    glu["c_insulin"] = standard(np.log(glu.LBXIN.where(glu.LBXIN.gt(0))))
    glu["c_homa"] = standard(np.log((glu.LBXGLU * glu.LBXIN / 22.5).where(glu.LBXGLU.gt(0) & glu.LBXIN.gt(0))))
    sp = pd.read_csv(CANDIDATE_DATA / "candidate_spirometry_raw.csv")
    sp["c_FEV1"] = standard(-sp.SPXNFEV1.where(sp.SPXNQEFF.isin(["A", "B"]) & sp.SPXNFEV1.gt(0)))  # lower is worse
    for frame, cols in ((crp, ["c_hsCRP"]), (glu, ["c_glucose", "c_insulin", "c_homa"]), (sp, ["c_FEV1"])):
        d = d.merge(frame[["SEQN"] + cols], on="SEQN", how="left", validate="one_to_one")
    ten = pd.read_csv(TEN_MARKERS).sort_values("rank").marker.tolist()[:10]
    candidates = {"hsCRP": "c_hsCRP", "fasting glucose": "c_glucose", "fasting insulin": "c_insulin",
                  "HOMA-IR": "c_homa", "FEV1": "c_FEV1"}
    return finish("modern", d, candidates, {"full": c.markers, "ten": ten})


def load_iii() -> dict:
    c = data.load_iii()
    d = c.people.copy()
    d[c.markers] = c.z.to_numpy()
    variables = {"SPPFEV1": ("FEV1", True), "SPPFVC": ("FVC", True), "BDPFNBMD": ("Femur neck BMD", True),
                 "BDPTOBMC": ("Total bone mineral content", True), "PEP12A1": ("Bioimpedance resistance", False),
                 "PEP12B1": ("Bioimpedance reactance", True), "PFPTTAND": ("Tandem stand time", True)}
    exam = data.read_nhanes3("exam", ["SEQN", "BMPHT"] + list(variables))
    adult = data.read_nhanes3("adult", ["SEQN", "HAR1", "HAR3"])
    d = d.merge(exam, on="SEQN", how="left", validate="one_to_one").merge(adult, on="SEQN", how="left", validate="one_to_one")
    d["height"] = d.BMPHT
    d["smoking"] = np.select([d.HAR1.eq(2), d.HAR1.eq(1) & d.HAR3.eq(1), d.HAR1.eq(1) & d.HAR3.eq(2)],
                             ["Never", "Current", "Former"], "Missing")
    d["smoking"] = d.smoking.replace("Missing", np.nan)
    units = d.SDPSTRA6.astype(str) + "_" + d.SDPPSU6.astype(str)
    d["fold"] = -1
    for f, (_, test) in enumerate(GroupKFold(10).split(d, groups=units)):
        d.loc[d.index[test], "fold"] = f
    candidates = {}
    for var, (name, reverse) in variables.items():
        d["c_" + var] = standard((-1 if reverse else 1) * d[var])     # sign: higher is worse
        candidates[name] = "c_" + var
    return finish("iii", d, candidates, {"full": c.markers})


def finish(era, d, candidates, panels) -> dict:
    """Drop participants with unknown smoking; list each candidate's evaluable folds."""
    d = d[d.smoking.notna()].reset_index(drop=True)
    folds = {}
    for name, col in candidates.items():
        has = d[col].notna() & (d.height.notna() if name in SPIROMETRY else True)
        folds[name] = [int(f) for f, g in d[has].groupby("fold") if era == "iii" or g.event.sum() >= 20]
    return {"era": era, "frame": d, "candidates": candidates, "panels": panels, "folds": folds}


# ===========================================================================
# Evaluation
# ===========================================================================
def risk_sets(era, d):
    return models.RiskSets(d.entry_age.to_numpy(), d.exit_age.to_numpy(), d.event.to_numpy(),
                           d.cycle.to_numpy(str) if era == "modern" else np.zeros(len(d), dtype=str))


def cox_fit(era, d, X):
    beta, ok = models.ridge_cox(risk_sets(era, d), X, np.zeros(X.shape[1]))
    if not ok:
        raise ArithmeticError("Cox fit did not converge")
    return beta


def direction(era, train, markers):
    """The fold's mortality direction: modern statsmodels-equivalent stratified fit; NHANES III lifelines,
    as in the reviewed analysis (see models.py)."""
    frame = train[["entry_age", "exit_age", "event", "female"] + markers + (["cycle"] if era == "modern" else [])]
    fit = models.cox(era, frame.reset_index(drop=True), ["female"] + markers)
    b = np.array([fit[m] for m in markers])
    return b / np.linalg.norm(b)


def baseline(c, name, spec, era):
    cols = [c.female.to_numpy(float)]
    if spec == "S1":
        cols += [c.age.to_numpy(float), c.smoking.eq("Former").to_numpy(float), c.smoking.eq("Current").to_numpy(float)]
    if name in SPIROMETRY and (spec == "S1" or era == "modern"):
        cols.append(c.height.to_numpy(float))
    return np.column_stack(cols)


def evaluate(package: dict, d: pd.DataFrame) -> list:
    era = package["era"]
    q_by_fold = {}
    for panel, markers in package["panels"].items():
        values = d[markers].to_numpy(float)
        q_by_fold[panel] = {f: values @ direction(era, d[d.fold.ne(f)], markers)
                            for f in sorted(package["frame"].fold.unique())}
    rows = []
    for name, col in package["candidates"].items():
        has = d[col].notna() & (d.height.notna() if name in SPIROMETRY else True)
        c, ids = d[has].reset_index(drop=True), np.flatnonzero(has)
        candidate = c[col].to_numpy(float)
        folds = package["folds"][name]
        deaths = int(c[c.fold.isin(folds)].event.sum())
        for spec in ("S1", "S2"):
            base = baseline(c, name, spec, era)
            for panel in package["panels"]:
                panel_gain = candidate_gain = 0.0
                improved = 0
                for f in folds:
                    test = c.fold.eq(f).to_numpy()
                    train = ~test
                    keep = [j for j in range(base.shape[1]) if base[train, j].std(ddof=1) > 1e-12]   # e.g. no current smokers
                    raw = np.column_stack([base[:, keep], q_by_fold[panel][f][ids], candidate])
                    X = (raw - raw[train].mean(axis=0)) / raw[train].std(axis=0, ddof=1)   # training scale
                    ll = []
                    for width in (len(keep), len(keep) + 1, len(keep) + 2):          # M0, M1, M2
                        beta = cox_fit(era, c[train], X[train, :width])
                        ll.append(risk_sets(era, c[test]).evaluate(X[test, :width], beta, order=0))
                    panel_gain += ll[1] - ll[0]
                    candidate_gain += ll[2] - ll[1]
                    improved += ll[2] > ll[1]
                rows.append({"era": era, "candidate": name, "baseline": spec, "panel": panel, "n": len(c),
                             "deaths": int(c.event.sum()), "deaths_evaluated": deaths, "folds": len(folds),
                             "folds_improved": improved, "panel_gain_per_1000": 1000 * panel_gain / deaths,
                             "candidate_gain_per_1000": 1000 * candidate_gain / deaths})
    return rows


def resample(d: pd.DataFrame, era: str, draw: int) -> pd.DataFrame:
    """Bootstrap: participants within cycle x sex (1999-2018), whole sampling units within strata (NHANES III).
    Each participant keeps their original fold."""
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 0 if era == "modern" else 1, draw]))
    if era == "modern":
        keys = d.cycle.astype(str) + "|" + d.female.astype(str)
        idx = np.concatenate([rng.choice(np.flatnonzero(keys.eq(k)), (keys == k).sum(), replace=True)
                              for k in sorted(keys.unique())])
    else:
        idx = []
        for _, g in d.groupby("SDPSTRA6", sort=True):
            units = sorted(g.SDPPSU6.unique())
            for unit in rng.choice(units, len(units), replace=True):
                idx.extend(g.index[g.SDPPSU6.eq(unit)].tolist())
        idx = np.asarray(idx, int)
    return d.iloc[idx].reset_index(drop=True)


_LOADED = {}


def one_draw(job) -> dict:
    """One bootstrap draw: every fit repeated on the resampled participants; both gains per comparison."""
    era, draw = job
    if era not in _LOADED:                       # each worker process loads an era once
        _LOADED[era] = load_modern() if era == "modern" else load_iii()
    package = _LOADED[era]
    out = {"era": era, "draw": draw}
    for r in evaluate(package, resample(package["frame"], era, draw)):
        stem = f"{r['candidate']}|{r['baseline']}|{r['panel']}"
        out[f"{stem}|candidate_gain"] = r["candidate_gain_per_1000"]
        out[f"{stem}|panel_gain"] = r["panel_gain_per_1000"]
    return out


def figure7(gains: pd.DataFrame, intervals: pd.DataFrame | None = None) -> None:
    """Figure 7: a, b each candidate's gain beyond the panel; c FEV1 beside the panel's own gain."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import displays as V

    g = gains[gains.baseline.eq("S1")]

    def bar(era, candidate, panel, quantity):
        r = g[g.era.eq(era) & g.candidate.eq(candidate) & g.panel.eq(panel)].iloc[0]
        low = high = np.nan
        if intervals is not None:
            hit = intervals[(intervals.era == era) & (intervals.key == f"{candidate}|S1|{panel}|{quantity}")]
            if len(hit):
                low, high = hit.iloc[0].low, hit.iloc[0].high
        return r[f"{quantity}_per_1000"], low, high, r

    fig = plt.figure(figsize=(185 * V.MM, 185 * V.MM), facecolor="white")
    order = {"modern": ["hsCRP", "fasting glucose", "fasting insulin", "HOMA-IR", "FEV1"],
             "iii": ["FEV1", "FVC", "Femur neck BMD", "Total bone mineral content", "Bioimpedance resistance",
                     "Bioimpedance reactance", "Tandem stand time"]}
    label = {"hsCRP": "hsCRP", "fasting glucose": "Fasting glucose", "fasting insulin": "Fasting insulin",
             "HOMA-IR": "HOMA-IR", "FEV1": "FEV$_1$", "FVC": "FVC", "Femur neck BMD": "Femoral-neck BMD",
             "Total bone mineral content": "Total bone mineral\ncontent", "Bioimpedance resistance": "Bioimpedance\nresistance",
             "Bioimpedance reactance": "Bioimpedance\nreactance", "Tandem stand time": "Tandem stand time"}
    for panel_letter, era, x0, title in (("a", "modern", .18, "NHANES 1999–2018"), ("b", "iii", .69, "NHANES III")):
        ax = fig.add_axes([x0, .53, .2, .33])
        names = order[era]
        for i, name in enumerate(names):
            y = len(names) - 1 - i
            primary = name == "FEV1"
            colour = V.PURPLE if primary else V.MUTED
            v, lo, hi, r = bar(era, name, "full", "candidate_gain")
            ax.hlines(y + (.12 if era == "modern" else 0), lo, hi, color=colour, lw=1.5)
            ax.plot(v, y + (.12 if era == "modern" else 0), "o", ms=6, mew=1.3, mec=colour,
                    mfc=colour if primary else "white")
            ax.text(1.02, y, f"{r.deaths_evaluated:,} / {r.n:,}", transform=ax.get_yaxis_transform(),
                    fontsize=6.5, color=V.MUTED, va="center")
            if era == "modern":
                v, lo, hi, _ = bar(era, name, "ten", "candidate_gain")
                ax.hlines(y - .12, lo, hi, color=colour, lw=1.5)
                ax.plot(v, y - .12, "D", ms=5.5, mew=1.3, mec=colour, mfc="white")
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels([label[n] for n in names][::-1])
        ax.set(xlim=(-12, 72), ylim=(-.6, len(names) - .4), xticks=[0, 20, 40, 60])
        ax.axvline(0, color=V.MUTED, lw=.8, ls=(0, (3, 2)))
        ax.grid(axis="x", color=V.GRID, lw=.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0, pad=6)
        ax.text(1.02, 1.03, "Evaluated deaths /\nparticipants", transform=ax.transAxes, fontsize=6.5, color=V.MUTED)
        fig.text(x0 - .155, .9, panel_letter, fontweight="bold", fontsize=10)
        fig.text(x0 - .12, .9, title, fontweight="bold", fontsize=8.6)
    fig.text(.5, .47, "Added held-out log partial likelihood (per 1,000 evaluated deaths)", ha="center", fontsize=8)
    c = fig.add_axes([.32, .1, .52, .27])
    rows = []
    for block, (era, name) in enumerate((("modern", "NHANES 2007–2012"), ("iii", "NHANES III"))):
        top = 3.6 - block * 2.4
        for k, (quantity, text, colour) in enumerate((("panel_gain", "Routine panel over baseline", V.MUTED),
                                                      ("candidate_gain", "FEV$_1$ over routine panel", V.PURPLE))):
            v, lo, hi, _ = bar(era, "FEV1", "full", quantity)
            c.hlines(top - k, lo, hi, color=colour, lw=1.6)
            c.plot(v, top - k, "o", color=colour, ms=6)
            rows.append((top - k, text))
        c.text(282, top + .5, name, ha="right", fontsize=7)
    c.axhline(1.9, color=V.GRID, lw=1)
    c.set_yticks([y for y, _ in rows])
    c.set_yticklabels([t for _, t in rows])
    c.set(xlim=(-12, 285), ylim=(-.4, 4.3), xlabel="Held-out log partial-likelihood gain\n(per 1,000 evaluated deaths)")
    c.axvline(0, color=V.MUTED, lw=.8, ls=(0, (3, 2)))
    c.grid(axis="x", color=V.GRID, lw=.5)
    c.set_axisbelow(True)
    c.tick_params(axis="y", length=0, pad=6)
    fig.text(.025, .41, "c", fontweight="bold", fontsize=10)
    fig.text(.06, .41, "FEV$_1$ increment relative to the routine-panel gain", fontweight="bold", fontsize=8.6)
    fig.legend(handles=[Line2D([0], [0], marker="o", lw=0, color=V.PURPLE, label="Primary FEV$_1$, full panel"),
                        Line2D([0], [0], marker="o", lw=0, mec=V.MUTED, mfc="white", label="Other full-panel tests (secondary)"),
                        Line2D([0], [0], marker="D", lw=0, mec=V.MUTED, mfc="white", label="Fixed ten-marker panel (secondary)")],
               loc="upper center", bbox_to_anchor=(.5, .995), ncol=3, frameon=False, fontsize=7)
    note = "" if intervals is not None else " Interval bars appear once the bootstrap has run."
    fig.text(.5, .015, "Candidate-specific samples; no ranking across tests." + note,
             ha="center", fontsize=6.5, color=V.MUTED)
    V.save_figure(fig, OUT, "figure7")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check-draws", type=int, default=0, metavar="N",
                        help="run only draws 0..N-1 per era and compare them with the archive (minutes)")
    parser.add_argument("--bootstrap", action="store_true", help="run all 2 x 200 bootstrap draws")
    parser.add_argument("--workers", type=int, default=8, help="processes for the draws")
    args = parser.parse_args()
    run = runs.start("s7")
    rows = []
    for package in (load_modern(), load_iii()):
        rows += evaluate(package, package["frame"])
        print(f"  {package['era']}: done", flush=True)
    gains = pd.DataFrame(rows)
    gains.to_csv(OUT / "gains.csv", index=False)

    def bootstrap_draws(numbers):
        jobs = [(era, b) for era in ("modern", "iii") for b in numbers]
        with ProcessPoolExecutor(min(len(jobs), args.workers)) as pool:
            return pd.DataFrame(list(pool.map(one_draw, jobs)))

    stored = OUT / "functional_bootstrap.csv"
    fresh = None
    if args.check_draws:
        fresh = bootstrap_draws(range(args.check_draws))
        run.write_check_draws(fresh, ["era", "draw"])
    if args.bootstrap:
        boot = bootstrap_draws(range(DRAWS))
        run.save_draws(boot, stored, ["era", "draw"], SEED)
    else:
        boot = run.stored_draws(stored, ["era", "draw"], SEED, runs.draw_keys(("modern", "iii"), range(DRAWS)), fresh)
    intervals = models.interval_table(boot) if boot is not None else None
    if intervals is not None:
        intervals.to_csv(OUT / "intervals.csv", index=False)
    figure7(gains, intervals)

    numbers = []
    for r in gains.itertuples():
        stem = f"{r.candidate}|{r.baseline}|{r.panel}"
        numbers += [(f"{stem}|candidate_gain", r.era, r.candidate_gain_per_1000),
                    (f"{stem}|panel_gain", r.era, r.panel_gain_per_1000),
                    (f"{stem}|n", r.era, r.n), (f"{stem}|deaths_evaluated", r.era, r.deaths_evaluated),
                    (f"{stem}|folds_improved", r.era, r.folds_improved), (f"{stem}|folds", r.era, r.folds)]
    numbers = pd.DataFrame(numbers, columns=["key", "era", "value"])
    if intervals is not None:
        numbers = numbers.merge(intervals[["era", "key", "low", "high"]], on=["era", "key"], how="left")
    else:
        numbers = numbers.assign(low=np.nan, high=np.nan)
    numbers.to_csv(OUT / "numbers.csv", index=False)
    run.finish(full_draws={"interval": DRAWS})       # per era
    print(f"done: outputs in {OUT}")


if __name__ == "__main__":
    main()
