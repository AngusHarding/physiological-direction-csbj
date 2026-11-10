# =============================================================================
# READ ME FIRST: s3_healthy_range.py  (Results 3)
# -----------------------------------------------------------------------------
# RESULTS 3     "Healthy-range values predicted mortality beyond abnormality
#               flags."
#
# WHAT IT DOES  Places every reading in one of five zones and asks, in held-out
#               survey cycles, whether exact values add mortality prediction
#               beyond knowing which readings are abnormal, and which zone the
#               added information comes from. Also compares adverse distance
#               (DA) with counts of abnormal markers. No mortality score q is
#               used here.
#
# METHODS       "Held-out prediction" (log partial likelihood per 1,000 deaths;
#               concordance); "Uncertainty"; "Survival models".
#
# MAKES         Figure 3a-d; Supp. Tables S3.3-S3.6; every Results 3 number.
#
# PSEUDOCODE
#   zone of each reading (person x marker): below limit / low corridor /
#     young-healthy / high corridor / above limit; clinical limits win
#   for each of 20 divisions of the ten cycles into 5 + 5, both directions:
#     design (which flag and slope columns exist) learned on the training half
#     penalty for each model chosen by leave-one-cycle-out inside training
#     fit each model on training; evaluate the log likelihood on the other half
#   models (all with sex + age):  A  nothing else
#                                 T  + clinical flags
#                                 TZ + corridor flags
#                                 Full + exact value within every zone
#                                 Full minus one zone group's values
#                                 single scores: DA, two abnormal counts
#   gain = (log lik. model 2 - log lik. model 1) / evaluated deaths x 1,000,
#          both halves pooled, then averaged over the 20 divisions
#   intervals: 200 participant draws (within cycle x sex), penalties fixed
# =============================================================================
"""Results 3: healthy-range values predicted mortality beyond abnormality flags.

The question
------------
A screen that registered only abnormal readings would discard whatever
information normal readings carry. Do exact marker values improve mortality
prediction beyond knowing which readings are clinically abnormal, and if so,
which readings carry the information?

How
---
Each reading (participant x marker) falls in one of five zones:
  below    below the clinical lower limit
  low      clinically normal but below the same-sex young-healthy interval ("low corridor")
  young    within the young-healthy interval (2.5th to 97.5th percentile)
  high     clinically normal but above the young-healthy interval ("high corridor")
  above    above the clinical upper limit
Clinical abnormality wins where intervals overlap. Limits are in
inputs/clinical_and_young_healthy_ranges.csv (mostly ABIM laboratory ranges).

Nested Cox models, all with sex and examination age:
  A      nothing else
  T      + clinical flags (below / above, per marker)
  TZ     + corridor flags (all five zones known)
  Full   + exact value within every zone (one slope per marker and zone)
  TZ+CX  Full without the young-zone slopes   (so Full - TZ+CX = what young-zone values add)
  TZ+YX  Full without the corridor slopes
  TZ+YC  Full without the beyond-limit slopes
and three single scores: adverse distance (DA), and the proportion of markers
beyond either clinical limit, or beyond the limit on the adverse side.

Every model is fitted on five survey cycles and evaluated on the other five,
both ways round, for 20 fixed divisions of the ten cycles
(inputs/cycle_partitions.csv). Models with marker terms are ridge-penalised;
each penalty is chosen by leave-one-cycle-out within the training half.
Gain = held-out log partial likelihood difference per 1,000 deaths, averaged
over the 20 divisions.

Scope
-----
Population   the NHANES 1999-2018 analysis cohort (11,497 participants, 38 markers).
Endpoint     accident-censored mortality.
Estimands    held-out log partial-likelihood gains (per 1,000 deaths, averaged over
             the 20 divisions) for the nested model comparisons above; held-out
             concordance of the main models; the share of readings in each zone.
Split        five cycles train, five evaluate, both ways, for 20 fixed divisions;
             penalties chosen inside the training half only.
Uncertainty  200 participant draws (within cycle x sex; one resample serves every
             division), with each division's tuned penalty held at its point value.
             The intervals therefore exclude penalty tuning and the choice of divisions.
Assumptions  the clinical and young-healthy limits, the transformations and the
             reference are fixed; ridge penalties shrink each zone's slopes.
Permitted    where in the range of each marker held-out mortality information sits,
interpretation  in likelihood terms. Not a causal statement, and not by itself a
             measure of clinical usefulness.

What it makes (in outputs/s3/)
------------------------------
  figure3.svg / .png      Figure 3
  zone_shares.csv         share of readings in each zone
  fits.csv                every held-out fit (division, half, model)
  gains.csv               the gains (per 1,000 deaths), with 95% intervals
  gain_bootstrap.csv      the 200 bootstrap draws behind the intervals (--bootstrap)
  numbers.csv             every Results 3 number quoted in the text

    python -B s3_healthy_range.py               # point estimates (penalty tuning: about 3 minutes on 12 cores)
    python -B s3_healthy_range.py --check-draws 2   # two bootstrap draws, to compare with the archive
    python -B s3_healthy_range.py --bootstrap   # also the 200-draw bootstrap (about 1.5 hours)
"""
from __future__ import annotations

import argparse
import os
import warnings
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
from statsmodels.duration.hazard_regression import PHReg   # noqa: E402

import data                               # noqa: E402
import displays as V                      # noqa: E402
import measures as M                      # noqa: E402
import models                             # noqa: E402
import runs                               # noqa: E402
import settings as S                      # noqa: E402

OUT = S.output_folder("s3")
INPUTS = S.PIPELINE / "inputs"
BELOW, LOW, YOUNG, HIGH, ABOVE = range(5)
ZONES = ("below", "low", "young", "high", "above")
GROUP = {BELOW: "X", LOW: "C", YOUNG: "Y", HIGH: "C", ABOVE: "X"}   # X beyond limits, C corridors, Y young
INWARD = {BELOW: LOW, LOW: YOUNG, ABOVE: HIGH, HIGH: YOUNG}         # where an unseen zone's readings go
SLOPES = {"TZ": "", "Full": "YCX", "TZ+CX": "CX", "TZ+YX": "YX", "TZ+YC": "YC"}   # which slope groups
SCORES = ("DA", "cut_either", "cut_adverse")
MODELS = ("A", "T", *SLOPES, *SCORES)
PENALISED = {"T", *SLOPES}
LAMBDA_GRID = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0)    # penalty per training death
BOOTSTRAP_SEED, BOOTSTRAP_DRAWS = 20261015, 200


# ===========================================================================
# Data: zones and scores (no outcome is used here)
# ===========================================================================
# -- INPUT: zones and scores (no outcome is used here) -----------------------
# Results 3: "In NHANES 1999-2018 (11,497 participants; 2,494 deaths; 38
#   markers), we placed each reading in one of five zones". The single scores
#   are DA and "the proportion of markers beyond either clinical limit, and the
#   proportion beyond the adverse-side limit".
class Data:
    def __init__(self):
        c = data.load_modern()
        p = c.people
        self.markers = c.markers
        self.cycle = p.cycle.to_numpy(str)
        self.sex = p.RIAGENDR.to_numpy(int)
        self.female = p.female.to_numpy(float)
        self.age = p.RIDAGEYR.to_numpy(float)
        self.entry, self.exit = p.entry_age.to_numpy(float), p.exit_age.to_numpy(float)
        self.event = p.event.to_numpy(int)
        self.coord = c.z.to_numpy(float)                 # SS coordinates, for the slopes
        raw = pd.read_csv(S.PREP / "07_master/analysis_master_raw.csv", usecols=["SEQN", *c.markers])
        self.raw = raw.set_index("SEQN").loc[p.SEQN, c.markers].to_numpy(float)
        self.zone = self.classify(pd.read_csv(INPUTS / "clinical_and_young_healthy_ranges.csv"))
        side = pd.read_csv(S.PREP_CONFIG / "fi_scoring_spec.csv").set_index("marker").loc[c.markers] \
            .adverse_direction.to_numpy(str)
        abnormal = (self.zone == BELOW) | (self.zone == ABOVE)
        adverse = np.where(side == "high_only", self.zone == ABOVE,
                           np.where(side == "low_only", self.zone == BELOW, abnormal))
        self.score = {"DA": M.adverse_distance(c.x), "cut_either": abnormal.mean(axis=1),
                      "cut_adverse": adverse.mean(axis=1)}
        parts = pd.read_csv(INPUTS / "cycle_partitions.csv")
        self.partitions = [(r.cycles_A.split("|"), r.cycles_B.split("|")) for r in parts.itertuples()]

    # -- CALCULATION: the five zones ---------------------------------------------
    # Results 3: "below or above the clinical limit; in a low or high corridor,
    #   clinically normal but outside the same-sex young-healthy interval; or
    #   within the young-healthy interval (Fig. 3a)". Figure 3 legend: "Clinical
    #   abnormality takes precedence where intervals overlap".
    # How: start everyone in the young zone; below the young-healthy 2.5th
    #   percentile -> low corridor; above the 97.5th -> high corridor; then the
    #   clinical limits are applied last, so they overwrite the corridors
    #   (precedence). Limits by marker and sex:
    #   inputs/clinical_and_young_healthy_ranges.csv (Supp. Tables S3.1-S3.2).
    def classify(self, ranges: pd.DataFrame) -> np.ndarray:
        ranges = ranges.set_index(["marker", "sex_code"])
        zone = np.full(self.raw.shape, YOUNG, dtype=np.int8)
        for j, marker in enumerate(self.markers):
            for s in (1, 2):
                r = ranges.loc[(marker, s)]
                rows = self.sex == s
                v = self.raw[rows, j]
                z = np.full(v.shape, YOUNG, dtype=np.int8)
                z[v < r.young_healthy_p2_5] = LOW
                z[v > r.young_healthy_p97_5] = HIGH
                # An "inclusive" limit is itself normal: only values beyond it are flagged.
                if pd.notna(r.clinical_healthy_low):
                    z[(v < r.clinical_healthy_low) if r.clinical_lower_inclusive else
                      (v <= r.clinical_healthy_low)] = BELOW
                if pd.notna(r.clinical_healthy_high):
                    z[(v > r.clinical_healthy_high) if r.clinical_upper_inclusive else
                      (v >= r.clinical_healthy_high)] = ABOVE
                zone[rows, j] = z
        return zone


# ===========================================================================
# The design learned from one training sample
# ===========================================================================
# -- CALCULATION: the model columns, learned from training rows only ---------
# Flags = one 0/1 column per non-young zone seen in training. Slopes = the
# exact value (coordinate minus the zone median) inside a zone, kept where the
# zone has at least 50 readings with 10 distinct values. Evaluation rows reuse
# the training design (Sample(..., design=train.design) in fit_half), so
# nothing about the evaluation half shapes the model.
class Design:
    """Which flag and slope columns exist, and their anchors, learned from training rows only.

    A flag is a 0/1 column for a non-young zone that occurs in training. A slope
    is a column (coordinate - zone median) inside one zone, kept only where the
    zone has at least 50 readings with at least 10 distinct values. A zone never
    seen in training is mapped to the nearest seen zone towards young.
    """

    def __init__(self, d: Data, rows: np.ndarray):
        zone, coord = d.zone[rows], d.coord[rows]
        k = zone.shape[1]
        self.observed = np.zeros((k, 5), bool)
        self.flags, self.slopes = [], []
        for m in range(k):
            for z in range(5):
                inside = zone[:, m] == z
                self.observed[m, z] = inside.any()
                if not inside.any():
                    continue
                if z != YOUNG:
                    self.flags.append((m, z))
                values = coord[inside, m]
                if inside.sum() >= 50 and np.unique(values).size >= 10:
                    self.slopes.append((m, z, float(np.median(values))))
        self.effective = np.tile(np.arange(5, dtype=np.int8), (k, 1))
        for m in range(k):
            for z in range(5):
                e = z
                while not self.observed[m, e]:
                    e = INWARD[e]
                self.effective[m, z] = e
        n_flags = len(self.flags)
        flag_zone = [z for _, z in self.flags]
        slope_zone = [z for _, z, _ in self.slopes]
        self.columns = {"A": [], "T": [i for i, z in enumerate(flag_zone) if z in (BELOW, ABOVE)]}
        for name, groups in SLOPES.items():
            self.columns[name] = list(range(n_flags)) + [n_flags + i for i, z in enumerate(slope_zone)
                                                         if GROUP[z] in groups]

    def matrix(self, d: Data, rows: np.ndarray) -> np.ndarray:
        zone, coord = d.zone[rows], d.coord[rows]
        X = np.zeros((len(rows), len(self.flags) + len(self.slopes)))
        for i, (m, z) in enumerate(self.flags):
            X[:, i] = self.effective[m, zone[:, m]] == z
        for i, (m, z, anchor) in enumerate(self.slopes):
            X[:, len(self.flags) + i] = np.where(zone[:, m] == z, coord[:, m] - anchor, 0.0)
        return X


class Sample:
    """One set of rows with the design, covariates and risk sets every model shares."""

    def __init__(self, d: Data, rows: np.ndarray, design: Design | None = None):
        self.rows = rows
        self.design = design or Design(d, rows)
        self.markers = self.design.matrix(d, rows)
        # sex and examination age (centred at 65 only for numerical conditioning)
        self.base = np.column_stack([d.female[rows], d.age[rows] - 65.0])
        self.entry, self.exit, self.event, self.strata = d.entry[rows], d.exit[rows], d.event[rows], d.cycle[rows]
        self.risk = models.RiskSets(self.entry, self.exit, self.event, self.strata)
        self.events = int(self.event.sum())

    def X(self, model: str) -> np.ndarray:
        return np.column_stack([self.base, self.markers[:, self.design.columns[model]]])


# Ridge penalty for each column: sex and age unpenalised, every marker column
# penalised by the same amount.
def penalty(X: np.ndarray, lam: float) -> np.ndarray:
    p = np.full(X.shape[1], lam)
    p[:2] = 0.0                                # sex and age are not penalised
    return p


# -- CALCULATION: choosing the penalty inside the training half ---------------
# Each training half holds five cycles; the penalty is chosen by
# leave-one-cycle-out within those five (best summed held-out log likelihood
# over the grid LAMBDA_GRID, extended if the best value is at an edge). The
# evaluation half is never used.
def tune(model: str, inner: list) -> float:
    """Choose the penalty by leave-one-cycle-out within the training half.

    Tries the grid from strongest to weakest (each fit starting from the last);
    if the best value is at an end of the grid, extends one decade at a time,
    twice at most. Returns the penalty per training death.
    """
    scores = {}

    def run(values):
        for lam in sorted(values, reverse=True):
            total = 0.0
            for fold in inner:
                X = fold["train"].X(model)
                beta, ok = models.ridge_cox(fold["train"].risk, X, penalty(X, lam * fold["train"].events),
                                            fold["warm"].get(model))
                if not ok:
                    total = -np.inf
                    break
                fold["warm"][model] = beta
                total += fold["test"].risk.evaluate(fold["test"].X(model), beta, order=0)
            scores[lam] = total
    run(LAMBDA_GRID)
    if not any(np.isfinite(v) for v in scores.values()):
        raise FloatingPointError(f"{model}: every penalty on the grid failed to converge in an inner fold; "
                                 f"no penalty can be chosen")
    for _ in range(2):
        best = max(scores, key=scores.get)
        if best == min(scores):
            run([best * 10 ** -0.5, best * 10 ** -1.0])
        elif best == max(scores):
            run([best * 10 ** 0.5, best * 10 ** 1.0])
        else:
            break
    return max(scores, key=scores.get)


# -- CALCULATION: held-out concordance ---------------------------------------
# Results 3: exact values "raised held-out concordance from 0.671 to 0.698".
# For each death, compare with every same-cycle participant still at risk;
# concordant if the person who died had the higher predicted risk.
def concordance(s: Sample, score: np.ndarray) -> tuple[float, float]:
    num = pairs = 0.0
    for g in np.unique(s.strata):
        idx = np.flatnonzero(s.strata == g)
        e, x, ev, h = s.entry[idx], s.exit[idx], s.event[idx], score[idx]
        dead = np.flatnonzero(ev == 1)
        comparable = (e[None, :] <= x[dead][:, None]) & (x[None, :] > x[dead][:, None])
        num += float(((h[dead][:, None] > h[None, :]) & comparable).sum()
                     + .5 * ((h[dead][:, None] == h[None, :]) & comparable).sum())
        pairs += float(comparable.sum())
    return num, pairs


# -- CALCULATION: one training half -> one evaluation half --------------------
# Results 3: "We compared nested Cox models on survey cycles withheld from
#   training: sex and age; then clinical abnormality flags; then corridor
#   flags; then exact values within each zone."
# How: build the design on training rows; evaluation rows use that design;
#   choose each penalty inside training (tune) unless given (bootstrap); fit
#   on training; log likelihood of the evaluation rows at the training
#   coefficients (PHReg(...).loglike(beta)); concordance.
def fit_half(d: Data, train_rows, eval_rows, lambdas: dict | None = None) -> dict:
    """Fit every model on the training rows and evaluate it on the evaluation rows.

    Without `lambdas`, penalties are tuned (point estimates); with them, the
    given penalties are reused (bootstrap).
    """
    train = Sample(d, train_rows)
    test = Sample(d, eval_rows, design=train.design)       # evaluation rows use the training design
    inner = []
    if lambdas is None:
        cycles = d.cycle[train_rows]
        for c in np.unique(cycles):
            fold = Sample(d, train_rows[cycles != c])
            inner.append({"train": fold, "test": Sample(d, train_rows[cycles == c], design=fold.design),
                          "warm": {}})
    out = {"events_eval": test.events}
    for model in MODELS:
        if model in SCORES:                                # one score, standardised on training rows
            s = d.score[model]
            mu, sd = s[train_rows].mean(), s[train_rows].std(ddof=1)
            Xtr = np.column_stack([train.base, (s[train_rows] - mu) / sd])
            Xte = np.column_stack([test.base, (s[eval_rows] - mu) / sd])
        else:
            Xtr, Xte = train.X(model), test.X(model)
        lam = 0.0
        if model in PENALISED:
            lam = lambdas[model] if lambdas is not None else tune(model, inner)
        beta, ok = models.ridge_cox(train.risk, Xtr, penalty(Xtr, lam * train.events))
        if not ok:
            raise ArithmeticError(f"{model} did not converge")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    # evaluates fixed coefficients; no fit happens here
            ll = float(PHReg(test.exit, Xte, status=test.event, entry=test.entry, strata=test.strata,
                             ties="efron").loglike(beta))
        num, pairs = concordance(test, Xte @ beta)
        out[model] = {"ll": ll, "penalty": lam, "conc_num": num, "conc_pairs": pairs}
    return out


# ===========================================================================
# Jobs: 20 divisions x 2 halves (x 200 draws for the bootstrap)
# ===========================================================================
_D = {}


def _start() -> None:
    _D["d"] = Data()


# -- UNCERTAINTY: one bootstrap draw -----------------------------------------
# Figure 3 legend: "Bars are 95% intervals from 200 participant-bootstrap
#   draws within cycle and sex, holding the selected penalties and
#   preprocessing fixed". One resample serves all 20 divisions.
def bootstrap_pool(d: Data, draw: int) -> np.ndarray:
    """Participants resampled with replacement within cycle x sex; one resample serves all divisions."""
    rng = np.random.default_rng([BOOTSTRAP_SEED, draw])
    parts = [rng.choice(idx, size=len(idx), replace=True)
             for c in np.unique(d.cycle) for s in (1, 2)
             for idx in [np.flatnonzero((d.cycle == c) & (d.sex == s))]]
    return np.sort(np.concatenate(parts))


# One job = one division, one direction (A trains or B trains), one draw.
# Divisions: inputs/cycle_partitions.csv (20 fixed divisions of the ten cycles).
def run_job(job: tuple) -> dict:
    d = _D["d"]
    division, half, draw, lambdas = job
    cycles_a, cycles_b = d.partitions[division]
    pool = np.arange(len(d.cycle)) if draw == 0 else bootstrap_pool(d, draw)
    in_a = np.isin(d.cycle[pool], cycles_a)
    train, test = (pool[in_a], pool[~in_a]) if half == "A" else (pool[~in_a], pool[in_a])
    return {"division": division + 1, "half": half, "draw": draw, **fit_half(d, train, test, lambdas)}


def flatten(results: list) -> pd.DataFrame:
    rows = []
    for r in results:
        for model in MODELS:
            rows.append({"draw": r["draw"], "division": r["division"], "half": r["half"],
                         "events_eval": r["events_eval"], "model": model, **r[model]})
    return pd.DataFrame(rows)


# -- CALCULATION: the gains reported in Results 3 ----------------------------
# Results 3: "Gains are held-out log partial likelihood per 1,000 deaths,
#   averaged over 20 divisions". "Clinical flags added 143 per 1,000 deaths to
#   sex and age" [flags_vs_sex_age = T - A]; "Corridor flags added 5.1"
#   [TZ - T]; "Exact values added 45.1 beyond all flags" [Full - TZ]; removing
#   young-zone values "lost 45.1" [Full - (Full without young slopes)];
#   "adverse distance (DA) added 122.9 per 1,000 deaths" [DA - A].
# How: within a division, sum log likelihoods over both halves and divide by
#   the evaluated deaths in both halves, x 1,000; main() then averages over
#   the 20 divisions and counts the divisions with a positive gain ("with a
#   gain in all 20 divisions").
def gains(fits: pd.DataFrame) -> pd.DataFrame:
    """Per draw and division: held-out log-likelihood differences per 1,000 deaths, both halves pooled."""
    out = []
    for (draw, division), block in fits.groupby(["draw", "division"]):
        ll = block.groupby("model").ll.sum()
        g = ll / block.drop_duplicates("half").events_eval.sum() * 1000
        c = block.groupby("model").conc_num.sum() / block.groupby("model").conc_pairs.sum()
        out.append({"draw": draw, "division": division,
                    "flags_vs_sex_age": g["T"] - g["A"],
                    "corridor_flags_vs_flags": g["TZ"] - g["T"],
                    "exact_values_vs_flags": g["Full"] - g["TZ"],
                    "drop_young": g["Full"] - g["TZ+CX"],
                    "drop_corridors": g["Full"] - g["TZ+YX"],
                    "drop_beyond_limits": g["Full"] - g["TZ+YC"],
                    **{f"{s}_vs_sex_age": g[s] - g["A"] for s in SCORES},
                    "DA_vs_cut_either": g["DA"] - g["cut_either"],
                    "DA_vs_cut_adverse": g["DA"] - g["cut_adverse"],
                    **{f"C_{m}": c[m] for m in ("TZ", "Full", "DA", "cut_either", "cut_adverse")}})
    return pd.DataFrame(out)


def fit_jobs(jobs: list, workers: int) -> list:
    with ProcessPoolExecutor(workers, initializer=_start) as pool:
        return list(pool.map(run_job, jobs, chunksize=1))


# ===========================================================================
# Figure 3
# ===========================================================================
# Draws Figure 3 from this run's numbers.
def figure3(shares: pd.Series, est: pd.DataFrame) -> None:
    """Figure 3: a, where readings lie; b, what exact values add; c, loss when detail is
    removed; d, adverse distance versus counts of abnormal markers. Gains per 1,000 deaths."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    YOUNG_C, CORRIDOR_C, BEYOND_C = V.GREEN, "#C9A227", "#B5452F"
    fig = plt.figure(figsize=(185 * V.MM, 185 * V.MM), facecolor="white")
    fig.text(.5, .985, "NHANES 1999–2018 · ages 52–79 · 11,497 participants · 2,494 endpoint deaths",
             ha="center", va="top", fontsize=7.5, color=V.MUTED)

    # a: the five zones (illustrative widths) and the share of readings in each group
    a = fig.add_axes([.055, .56, .42, .33])
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.set_axis_off()
    widths, left = [.16, .2, .28, .2, .16], 0.0
    fills = ["#F0D5D2", "#F4EAD2", "#CBE9E0", "#F4EAD2", "#F0D5D2"]
    names = ["Below\nclinical\nlimit", "Low\ncorridor", "Young-\nhealthy\nzone", "High\ncorridor", "Above\nclinical\nlimit"]
    for w, fill, name in zip(widths, fills, names):
        a.add_patch(Rectangle((left, .62), w, .2, facecolor=fill, edgecolor="white", lw=1))
        a.text(left + w / 2, .72, name, ha="center", va="center", fontsize=7)
        left += w
    a.text(0, .9, "436,886 readings = 11,497 people × 38 markers", fontsize=7,
           color=V.MUTED)
    a.plot([.16, .16, .84, .84], [.585, .56, .56, .585], color=V.MUTED, lw=.8)
    a.text(.5, .52, "Clinically normal range", ha="center", va="top", fontsize=7.5)
    for y, colour, text in ((.38, YOUNG_C, f"Young-healthy: {100 * shares['young']:.1f}% of readings"),
                            (.29, CORRIDOR_C, f"Corridors: {100 * (shares['low'] + shares['high']):.1f}% of readings"),
                            (.20, BEYOND_C, f"Beyond limits: {100 * (shares['below'] + shares['above']):.1f}% of readings")):
        a.add_patch(Rectangle((.01, y - .02), .025, .04, color=colour))
        a.text(.06, y, text, va="center", fontsize=8)
    a.text(0, .1, "Illustrative widths; some corridors may be empty.", fontsize=6.5, color=V.MUTED)

    def forest(ax, rows, xlabel, separator=None):
        e = est.set_index("key")
        for i, (key, label, colour, paired) in enumerate(rows):
            y = len(rows) - 1 - i
            r = e.loc[key]
            ax.hlines(y, r.low, r.high, color=colour, lw=1.6)
            ax.plot(r.estimate, y, marker="D" if paired else "o", ms=6, mew=1.3, mec=colour,
                    mfc="white" if paired else colour)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([label for _, label, _, _ in rows][::-1])
        ax.set(xlim=(-15, 170), ylim=(-.6, len(rows) - .4), xlabel=xlabel)
        ax.set_xticks([0, 50, 100, 150])
        ax.axvline(0, color=V.MUTED, lw=.8, ls=(0, (3, 2)))
        if separator is not None:
            ax.axhline(separator, color=V.GRID, lw=1)
        ax.grid(axis="x", color=V.GRID, lw=.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0, pad=6)

    forest(fig.add_axes([.75, .60, .235, .28]),
           [("flags_vs_sex_age", "Clinical flags vs\nsex + age", V.INK, False),
            ("corridor_flags_vs_flags", "Corridor flags vs\nclinical flags", CORRIDOR_C, False),
            ("exact_values_vs_flags", "Exact within-zone values\nvs all flags", YOUNG_C, False)],
           "Prediction gain\n(per 1,000 deaths)")
    forest(fig.add_axes([.235, .16, .235, .29]),
           [("drop_young", "Young-healthy\ndetail omitted", YOUNG_C, False),
            ("drop_corridors", "Corridor\ndetail omitted", CORRIDOR_C, False),
            ("drop_beyond_limits", "Beyond-limit\ndetail omitted", BEYOND_C, False)],
           "Prediction loss\n(per 1,000 deaths)")
    d = fig.add_axes([.75, .16, .235, .29])
    forest(d, [("DA_vs_sex_age", "Adverse distance (DA)", YOUNG_C, False),
               ("cut_either_vs_sex_age", "Binary: either side", V.MUTED, False),
               ("cut_adverse_vs_sex_age", "Binary: adverse side", V.MUTED, False),
               ("DA_vs_cut_either", "DA − binary\n(either side)", YOUNG_C, True),
               ("DA_vs_cut_adverse", "DA − binary\n(adverse side)", YOUNG_C, True)],
           "Gain or paired difference\n(per 1,000 deaths)", separator=1.5)
    d.text(.98, .98, "First three: gain over sex + age", transform=d.transAxes, ha="right", va="top",
           fontsize=6.5, color=V.MUTED)
    for key, x, y, title in (("a", .025, .935, "Where the readings lie"), ("b", .555, .935, "What exact values add"),
                             ("c", .025, .51, "Loss when detail is removed"),
                             ("d", .555, .51, "Distance versus counts of abnormal markers")):
        fig.text(x, y, key, fontweight="bold", fontsize=10, va="center")
        fig.text(x + .033, y, title, fontsize=8.6, fontweight="bold", va="center")
    fig.text(.5, .025, "Held-out Cox log partial likelihood; all forest panels use the same scale.",
             ha="center", fontsize=7, color=V.MUTED)
    V.save_figure(fig, OUT, "figure3")


# 95% interval = 2.5th-97.5th percentiles of the 200 draws, each draw already
# averaged over the 20 divisions.
def intervals(boot: pd.DataFrame) -> dict:
    """95% percentile interval (numpy linear rule) of each quantity over the bootstrap draws.

    Each draw's value is already the mean over the 20 cycle divisions (draw_gains), as in
    the reviewed run (R077), so the interval describes resampling of participants with the
    divisions and the tuned penalties held fixed.
    """
    return {key: models.percentile_interval(boot[key]) for key in boot.columns if key != "draw"}


# ===========================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--bootstrap", action="store_true", help="redo the 200 bootstrap draws (hours)")
    parser.add_argument("--check-draws", type=int, default=0, metavar="N",
                        help="run only bootstrap draws 1..N, to compare with the archive (minutes)")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    run = runs.start("s3")

    d = Data()
    shares = pd.Series({ZONES[z]: float((d.zone == z).mean()) for z in range(5)})
    abnormal = ((d.zone == BELOW) | (d.zone == ABOVE)).sum(axis=1)
    shares.rename("share").to_csv(OUT / "zone_shares.csv")

    # The point fits (with their tuned penalties) take about ten minutes; they are reused only
    # when their manifest shows the same data, seed and code (runs.current_result).
    fits_file = OUT / "fits.csv"
    point = run.current_result(fits_file, BOOTSTRAP_SEED)
    if point is None:
        print("point estimates: 20 divisions x 2 halves, penalties tuned")
        point = flatten(fit_jobs([(k, h, 0, None) for k in range(20) for h in "AB"], args.workers))
        run.save_result(point, fits_file, BOOTSTRAP_SEED)
    per_division = gains(point)

    def draw_gains(numbers):
        """Bootstrap draws with each division's tuned penalties held fixed; mean over divisions."""
        lambdas = {(r.division - 1, r.half): {} for r in point.itertuples()}
        for r in point[point.model.isin(PENALISED)].itertuples():
            lambdas[(r.division - 1, r.half)][r.model] = r.penalty
        jobs = [(k, h, b, lambdas[(k, h)]) for b in numbers for k in range(20) for h in "AB"]
        return gains(flatten(fit_jobs(jobs, args.workers))).groupby("draw").mean(numeric_only=True) \
            .drop(columns="division").reset_index()

    stored = OUT / "gain_bootstrap.csv"
    fresh = None
    if args.check_draws:
        fresh = draw_gains(range(1, args.check_draws + 1))
        run.write_check_draws(fresh, ["draw"])
        print(f"{args.check_draws} draws written to check_draws.csv")
    if args.bootstrap:
        boot = draw_gains(range(1, BOOTSTRAP_DRAWS + 1))
        run.save_draws(boot, stored, ["draw"], BOOTSTRAP_SEED)
    else:
        boot = run.stored_draws(stored, ["draw"], BOOTSTRAP_SEED,
                                runs.draw_keys(None, range(1, BOOTSTRAP_DRAWS + 1), era_column=False), fresh)

    rows = []
    bounds = intervals(boot) if boot is not None else {}
    for key in [c for c in per_division.columns if c not in ("draw", "division")]:
        v = per_division[key]
        low, high = bounds.get(key, (np.nan, np.nan))
        rows.append({"key": key, "estimate": v.mean(), "low": low, "high": high,
                     "divisions_positive": int((v > 0).sum())})
    est = pd.DataFrame(rows)
    est.to_csv(OUT / "gains.csv", index=False)
    figure3(shares, est)                           # interval bars appear once the bootstrap has run

    numbers = [("readings", d.zone.size, np.nan, np.nan), ("share_young_pct", 100 * shares["young"], np.nan, np.nan),
               ("share_beyond_limits_pct", 100 * (shares["below"] + shares["above"]), np.nan, np.nan),
               ("participants_no_abnormal", int((abnormal == 0).sum()), np.nan, np.nan)]
    for r in est.itertuples():
        numbers.append((r.key, r.estimate, r.low, r.high))
        if not r.key.startswith("C_"):
            numbers.append((f"{r.key}_divisions_positive", r.divisions_positive, np.nan, np.nan))
    pd.DataFrame(numbers, columns=["key", "value", "low", "high"]).assign(era="modern") \
        .to_csv(OUT / "numbers.csv", index=False)
    run.finish(full_draws={"interval": BOOTSTRAP_DRAWS})
    print(f"done: outputs in {OUT}")


if __name__ == "__main__":
    main()
