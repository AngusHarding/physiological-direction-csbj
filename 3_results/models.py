# =============================================================================
# READ ME FIRST: models.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  The statistical machinery, each piece written once: Cox models,
#               the women-minus-men contrast, a fast Cox solver (with an
#               optional ridge penalty), the cross-fitted mortality score q,
#               bootstraps, intervals, concordance, weighted quantiles and
#               Kaplan-Meier risk. It chooses no analysis; the s1-s7 scripts do.
#
# METHODS       "Statistical analysis": "Survival models", "Sex differences",
#               "Held-out prediction", "Screening risks", "Uncertainty";
#               "Physiological distance and the mortality score":
#               "Mortality score".
#
# USED BY       cox, sex_contrast, bootstraps: Results 1 and 2.
#               mortality_score (q, cross-fitted): Results 2.
#               RiskSets / ridge_cox: Results 3 (penalised); Results 2 (Table 2),
#                 5, 6 and 7 with the penalty set to zero (an ordinary Cox fit).
#               survey_bootstrap_weights, concordance_within: Results 5.
#               kaplan_meier_risk, weighted_quantile: Results 6.
#               interval_table: Results 6 and 7.
#
# PSEUDOCODE
#   cox(era, terms):
#     NHANES 1999-2018: statsmodels PHReg, time = attained age, entry = entry
#       age, Efron ties, a separate baseline hazard per survey cycle
#     NHANES III: lifelines, the same time scale and entry, unstratified
#   sex_contrast: OLS of (value standardised to cohort SD) on
#     female + age spline + survey period; return the female coefficient
#   mortality_score: for each fold: fit Cox(female + all coordinates) on the
#     training rows; u = marker coefficients / their length;
#     q(test rows) = coordinates . u   (female coefficient not used)
#   bootstrap: resample rows with replacement within cycle x sex
#     (sex x phase in NHANES III); 95% interval = 2.5th, 97.5th percentiles
# =============================================================================
"""The statistical models, each written once.

Cox regression (each era's reviewed program), the women-versus-men contrast,
ridge-penalised Cox with a vectorised likelihood (Results 3-7), the
cross-fitted mortality score q, the participant and survey-design bootstraps,
percentile intervals, concordance within sex and period, weighted quantiles
and Kaplan-Meier risk. This file chooses no analysis; the section scripts do.

Two Cox programs, on purpose
----------------------------
The reviewed analysis fitted NHANES 1999-2018 with statsmodels (stratified by
survey cycle) and NHANES III with lifelines (unstratified). The two programs
differ in one detail: whether someone who enters the study at exactly the age
another participant dies counts as being at risk at that death. That shifts
estimates in the third decimal place, so the pipeline keeps each era's program
in order to reproduce the reviewed numbers exactly (decision AD-012).

Fit problems are never hidden (R120-F10)
----------------------------------------
Statsmodels reports trouble only as warnings, and lifelines as warnings or an
error. The reviewed analysis kept these fits and recorded their warnings; it
did not discard them. The pipeline does the same: every warning or failed
convergence is written, with the size of the model, to the section's run
record (see `note_fit` and runs.py), and a fit whose coefficients are not
finite numbers stops the run, because no estimate can be built on it.
"""
from __future__ import annotations

import json
import os
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from patsy import dmatrix
from statsmodels.duration.hazard_regression import PHReg

import settings as S

DIAGNOSTICS_ENV = "NHANES_FIT_DIAGNOSTICS"     # folder for this run's fit-problem records (set by runs.start)


def note_fit(kind: str, caught, converged: bool = True, **context) -> None:
    """Record a fit's warnings or failed convergence in the run's diagnostics folder.

    Each process appends to its own file, so bootstrap workers' problems are kept
    too; runs.Run.finish gathers them into run.json. Clean fits write nothing.
    """
    messages = sorted({str(w.message) for w in caught})
    if converged and not messages:
        return
    folder = os.environ.get(DIAGNOSTICS_ENV)
    if not folder:
        return
    entry = {"kind": kind, "converged": bool(converged), "warnings": messages, **context}
    with open(os.path.join(folder, f"{os.getpid()}.jsonl"), "a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, default=str) + "\n")


def require_finite(coef, what: str) -> None:
    """Stop if a fit returned coefficients that are not finite numbers."""
    if not np.isfinite(np.asarray(coef, float)).all():
        raise FloatingPointError(f"{what}: the fit returned non-finite coefficients; no estimate can use it")


# -- CALCULATION: Cox model on attained age ----------------------------------
# Methods, "Survival models": "Cox models with attained age as the time scale,
#   delayed entry at examination age (interview age in NHANES III) [48],
#   Efron's approximation for ties [49], and survey-cycle-specific baseline
#   hazards in NHANES 1999-2018."
# How: NHANES 1999-2018 uses statsmodels PHReg with strata = survey cycle,
#   entry = entry age, ties = "efron". NHANES III uses lifelines (Efron ties by
#   default), unstratified, with entry_col = entry age. Returns log hazard
#   ratios (exp gives the hazard ratios in Fig. 1d and Fig. 2d).
def cox(era: str, frame: pd.DataFrame, terms: list[str]) -> dict:
    """Fit a Cox model of the endpoint on attained age; return {term: log hazard ratio}.

    `frame` needs entry_age, exit_age, event and the columns named in `terms`.
    Also returns standard errors under "se". Warnings are recorded (note_fit).
    """
    context = {"era": era, "terms": len(terms), "rows": len(frame), "events": int(frame.event.sum())}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if era == "modern":
            fit = PHReg(endog=frame.exit_age.to_numpy(float), exog=frame[terms],
                        status=frame.event.to_numpy(int), entry=frame.entry_age.to_numpy(float),
                        strata=frame.cycle.to_numpy(str), ties="efron").fit(disp=0)
            coef, se = np.asarray(fit.params, float), np.asarray(fit.bse, float)
        else:
            f = CoxPHFitter()
            f.fit(frame[["exit_age", "entry_age", "event"] + terms], duration_col="exit_age",
                  entry_col="entry_age", event_col="event")
            coef, se = f.params_[terms].to_numpy(float), f.standard_errors_[terms].to_numpy(float)
    note_fit("cox", caught, **context)
    require_finite(coef, f"Cox model ({context})")
    return {**dict(zip(terms, coef)), "se": dict(zip(terms, se))}


# -- CALCULATION: women-minus-men difference ---------------------------------
# Methods, "Sex differences": "women-minus-men differences in distances and
#   scores were estimated by linear regression on sex, an age spline and survey
#   cycle or phase."
# How: the value is first standardised to cohort SD; ordinary least squares on
#   female + age spline + period; the female coefficient is the difference in
#   SD. Splines and periods are in settings.ERAS.
# Results 1 (Fig. 1c) and Results 2 (Fig. 2c) report these numbers.
def sex_contrast(era: str, frame: pd.DataFrame, values: np.ndarray) -> float:
    """Women minus men in `values`, in cohort SDs, adjusted for age and period.

    Ordinary least squares of the standardised values on female, a flexible
    age curve (settings.ERAS[era]["age_spline"]) and the survey period.
    """
    e = S.ERAS[era]
    X = dmatrix(f"female + {e['age_spline']} + C({e['period']})", frame, return_type="dataframe")
    y = (values - values.mean()) / values.std(ddof=1)
    return float((np.linalg.pinv(X.to_numpy(float)) @ y)[list(X.columns).index("female")])


# ---------------------------------------------------------------------------
# Ridge-penalised Cox regression (Results 3)
# ---------------------------------------------------------------------------
# -- CALCULATION: the Cox partial likelihood ---------------------------------
# The same Cox likelihood as above (attained age, delayed entry, Efron ties,
# strata), coded so that it can be evaluated quickly many times. It is used to
# fit models (ridge_cox) and to score held-out data at fixed coefficients.
# Methods, "Held-out prediction": "evaluated by the evaluation-set Cox log
#   partial likelihood at training coefficients, reported as gain per 1,000
#   evaluated deaths". The "gain" is the difference of two such likelihoods,
#   both computed by evaluate(..., order=0) on the evaluation rows.
class RiskSets:
    """The Cox partial likelihood for one set of participants, ready for repeated evaluation.

    Efron ties, delayed entry (at risk while entry <= t <= exit, as statsmodels)
    and strata. The risk sets are coded once as a sparse +1/-1 matrix over the
    ordered death times, so the likelihood, its gradient and its Hessian for any
    coefficients are a few matrix products instead of a loop over deaths.
    """

    def __init__(self, entry, exit_, event, strata, weights=None):
        from scipy import sparse
        entry, exit_, event = map(np.asarray, (entry, exit_, event))
        n = len(entry)
        # Case weights (the survey bootstrap's replicate weights). With weights, each tied
        # death time contributes (sum of the dying's weights / number dying) per Efron step.
        # Rows with zero weight must be left out by the caller.
        self.w = np.ones(n) if weights is None else np.asarray(weights, float)
        lo, hi, death_k = np.full(n, -1), np.full(n, -2), np.full(n, -1)
        offset, counts = 0, []
        for s in np.unique(strata):
            idx = np.flatnonzero(strata == s)
            times = np.unique(exit_[idx][event[idx] == 1])
            if len(times) == 0:
                continue
            lo[idx] = offset + np.searchsorted(times, entry[idx], side="left")     # first death time at risk
            hi[idx] = offset + np.searchsorted(times, exit_[idx], side="right") - 1  # last death time at risk
            dead = idx[event[idx] == 1]
            death_k[dead] = offset + np.searchsorted(times, exit_[dead])
            counts.append(np.bincount(death_k[dead] - offset, minlength=len(times)))
            offset += len(times)
        self.K, self.n = offset, n
        valid = (lo >= 0) & (lo <= hi)
        rows = np.flatnonzero(valid)
        self.diff = sparse.csr_matrix((np.r_[np.ones(len(rows)), -np.ones(len(rows))],
                                       (np.r_[lo[rows], hi[rows] + 1], np.r_[rows, rows])), shape=(self.K + 1, n))
        self.lo, self.hi, self.valid = lo, hi, valid
        self.events = np.flatnonzero(event == 1)
        self.death_k = death_k[self.events]
        self.Dd = sparse.csr_matrix((np.ones(len(self.events)), (self.death_k, self.events)), shape=(self.K, n))
        d = np.concatenate(counts) if counts else np.zeros(0, int)
        self.kk = np.repeat(np.arange(self.K), d)                       # death time of each death
        self.frac = np.concatenate([np.arange(c) / c for c in d]) if len(d) else np.zeros(0)   # Efron fractions
        dead_weight = np.bincount(self.death_k, weights=self.w[self.events], minlength=self.K)
        self.fac = (dead_weight / np.maximum(d, 1))[self.kk]             # 1 for every step when unweighted

    def evaluate(self, X, beta, order: int = 2):
        """Log likelihood; with order >= 1 also the gradient; with order 2 also the Hessian."""
        eta = X @ beta
        shift = float(eta.max()) if len(eta) else 0.0
        r = self.w * np.exp(eta - shift)
        S0 = np.cumsum(self.diff @ r)[: self.K]
        S0D = self.Dd @ r
        c = S0[self.kk] - self.frac * S0D[self.kk]
        we = self.w[self.events]
        loglik = float((we * eta[self.events]).sum() - (self.fac * (np.log(c) + shift)).sum())
        if order == 0:
            return loglik
        rX = r[:, None] * X
        S1 = np.cumsum(self.diff @ rX, axis=0)[: self.K]
        m = (S1[self.kk] - self.frac[:, None] * (self.Dd @ rX)[self.kk]) / c[:, None]
        score = (we[:, None] * X[self.events]).sum(axis=0) - (self.fac[:, None] * m).sum(axis=0)
        if order == 1:
            return loglik, score
        inv = self.fac / c
        A = np.bincount(self.kk, weights=inv, minlength=self.K)
        B = np.bincount(self.kk, weights=self.frac * inv, minlength=self.K)
        cumA = np.r_[0.0, np.cumsum(A)]
        a = np.zeros(self.n)
        a[self.valid] = cumA[self.hi[self.valid] + 1] - cumA[self.lo[self.valid]]
        w = r * a
        w[self.events] -= r[self.events] * B[self.death_k]
        return loglik, score, -((X.T * w) @ X - (self.fac[:, None] * m).T @ m)


# -- CALCULATION: fitting a Cox model, optionally with a ridge penalty -------
# Maximises (Cox log likelihood - 0.5 x sum of penalty x beta^2) by Newton
# steps. With the penalty set to zero, as in Results 2 (Table 2), 5, 6 and 7,
# this is an ordinary unpenalised Cox fit. Only Results 3 uses a positive
# penalty, chosen within training data (s3_healthy_range.tune).
def ridge_cox(risk: RiskSets, X, penalty, start=None):
    """Maximise log likelihood - 0.5 * sum(penalty * beta^2) by Newton steps with step-halving.

    `penalty` gives each column its own weight (zero = unpenalised). Converged
    when the largest step is below 1e-9 and the largest penalised gradient below
    1e-6. Returns (beta, converged).
    """
    from scipy.linalg import cho_factor, cho_solve
    p = X.shape[1]
    beta = np.zeros(p) if start is None else start.copy()
    ll, g, H = risk.evaluate(X, beta)
    objective = ll - 0.5 * float((penalty * beta ** 2).sum())
    for _ in range(50):
        gradient = g - penalty * beta
        curvature = -H + np.diag(penalty)
        try:
            step = cho_solve(cho_factor(curvature), gradient)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(curvature, gradient, rcond=None)[0]
        t = 1.0
        while True:                                   # halve the step until the objective does not fall
            trial = beta + t * step
            new = risk.evaluate(X, trial, order=0) - 0.5 * float((penalty * trial ** 2).sum())
            if np.isfinite(new) and new >= objective - 1e-10 * max(1.0, abs(objective)):
                break
            t *= 0.5
            if t < 1e-10:
                trial = beta
                break
        moved = float(np.max(np.abs(trial - beta))) if p else 0.0
        beta = trial
        ll, g, H = risk.evaluate(X, beta)
        objective = ll - 0.5 * float((penalty * beta ** 2).sum())
        largest = float(np.max(np.abs(g - penalty * beta))) if p else 0.0
        if not (np.isfinite(objective) and np.isfinite(beta).all()):
            return beta, False
        if (moved < 1e-9 and largest < 1e-6) or (p and largest < 1e-9):
            return beta, True
        if moved == 0.0:
            return beta, False
    return beta, False


# ---------------------------------------------------------------------------
# The mortality score q (cross-fitted)
# ---------------------------------------------------------------------------
# -- CROSS-FITTING: which rows are held out ----------------------------------
# Results 2: each participant is scored "with a model fitted without their own
#   survey cycle (NHANES 1999-2018) or sampling-unit fold (NHANES III)".
# How: NHANES 1999-2018 holds out one survey cycle at a time (ten folds).
#   NHANES III makes ten folds of whole sampling units (stratum x PSU), so no
#   unit is split between training and test.
def folds(era: str, people: pd.DataFrame):
    """Training and test rows for cross-fitting, as (train, test) row-number arrays.

    NHANES 1999-2018: leave one survey cycle out (ten folds).
    NHANES III: ten folds of whole sampling units (stratum x primary sampling
    unit), so neighbours from the same sampling unit are never split between
    training and test.
    """
    if era == "modern":
        for cycle in S.CYCLES:
            test = people.cycle.eq(cycle).to_numpy()
            yield np.flatnonzero(~test), np.flatnonzero(test)
    else:
        from sklearn.model_selection import GroupKFold
        units = (people.SDPSTRA6.astype(str) + "_" + people.SDPPSU6.astype(str)).to_numpy()
        yield from GroupKFold(n_splits=10).split(people, groups=units)


# -- CALCULATION: the mortality score q --------------------------------------
# Methods, "Mortality score": "The score came from a Cox model of the principal
#   endpoint containing a female indicator and all marker coordinates, fitted on
#   both sexes together [47]. Its biomarker coefficients beta define the
#   mortality-associated direction u = beta/||beta||. The mortality score
#   q = u'z ... is a profile's signed position along that direction, and the
#   fitted female coefficient is omitted from it."
#   "Wherever q was used, its direction was fitted without the scored
#   participant's outcome."
# How: for each fold (see folds above): fit Cox(female + all z) on training
#   rows; u = marker coefficients / their length; q = z . u for the test rows.
#   This shared version is used by Results 2. Results 5, 6 and 7 fit q in
#   their own scripts, and Results 4 builds full and reduced q scores in
#   s4_panel.py (see each script's READ ME block). Results 3 uses no q score.
def mortality_score(era: str, people: pd.DataFrame, X: np.ndarray, units: list | None = None):
    """q: each participant's coordinates projected on a mortality direction learned without them.

    In each fold, a Cox model of the endpoint on female + all marker
    coordinates is fitted to the training rows; the marker coefficients,
    rescaled to length one, are the direction u. Test participants get
    q = coordinates . u. The sex coefficient is used in fitting but left out
    of q.

    Returns (q, units). Pass `units` (one direction per fold, from an earlier
    call) to project new coordinates on directions already learned instead of
    refitting.
    """
    q = np.full(len(X), np.nan)
    learned = []
    names = [f"m{j}" for j in range(X.shape[1])]
    for k, (train, test) in enumerate(folds(era, people)):
        if units is None:
            frame = people.iloc[train].reset_index(drop=True).copy()
            frame[names] = X[train]
            fit = cox(era, frame, ["female"] + names)
            beta = np.array([fit[n] for n in names])
            norm = np.linalg.norm(beta)
            if not norm > 1e-12:                   # no direction to project on (cf. s4 NEAR_ZERO_NORM)
                raise FloatingPointError(f"mortality direction in fold {k} has norm {norm}")
            u = beta / norm
        else:
            u = units[k]
        learned.append(u)
        q[test] = X[test] @ u
    return q, learned


# ---------------------------------------------------------------------------
# Participant bootstrap
# ---------------------------------------------------------------------------
# -- UNCERTAINTY: the participant bootstrap (Results 1 and 2) ----------------
# Methods, "Uncertainty": "analysis-specific bootstraps that resampled either
#   participants within sex and survey cycle, or primary sampling units within
#   strata [53], refitting the stated learned stages".
# How: NHANES 1999-2018 resamples participants within survey cycle x sex;
#   NHANES III within sex x survey phase. Each group keeps its size. Draws:
#   settings.BOOTSTRAP_DRAWS (500). What is refitted in each draw is set by the
#   calling script (Supp. Table S0.7).
def bootstrap_strata(era: str, people: pd.DataFrame) -> list[np.ndarray]:
    """Groups resampled separately, so every draw keeps each group's size.

    NHANES 1999-2018: survey cycle by sex (cycles in time order, men first).
    NHANES III: sex by survey phase.
    """
    if era == "modern":
        return [np.flatnonzero(people.cycle.eq(c).to_numpy() & people.female.eq(f).to_numpy())
                for c in S.CYCLES for f in (0, 1)]
    keys = people.female.astype(str) + "|" + people.phase.astype(str)
    return [np.flatnonzero(keys.to_numpy() == k) for k in sorted(keys.unique())]


def bootstrap_rows(era: str, strata: list[np.ndarray], draw: int) -> np.ndarray:
    """Row numbers for bootstrap draw `draw` (1, 2, ...): each group resampled with replacement.

    The random numbers depend only on the era's seed and the draw number, so
    any draw can be regenerated on its own.
    """
    rng = np.random.default_rng(np.random.SeedSequence([S.ERAS[era]["seed"], draw]))
    return np.concatenate([rng.choice(s, size=len(s), replace=True) for s in strata])


# -- UNCERTAINTY: the 95% interval -------------------------------------------
# The 2.5th and 97.5th percentiles of the bootstrap draws: the "95% interval"
# reported throughout the Results for bootstrap-based estimates.
def percentile_interval(values) -> tuple[float, float]:
    """The 2.5th and 97.5th percentiles (numpy's linear rule): the 95% bootstrap interval."""
    v = np.asarray(values, float)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def interval_table(draws: pd.DataFrame, max_missing: float = 0.01) -> pd.DataFrame:
    """The 95% percentile interval of every drawn quantity, per era (Results 6 and 7).

    `draws` has one row per (era, draw) and one column per quantity. A draw in which
    a quantity is missing is left out for that quantity, and the number of draws used
    is recorded. If more than `max_missing` of the draws are missing (1%, the rule of
    the reviewed Results 7 run), the interval is left missing rather than computed
    from a selected subset of draws.
    """
    rows = []
    for era, g in draws.groupby("era", sort=False):
        for key in [c for c in g.columns if c not in ("era", "draw")]:
            v = g[key].to_numpy(float)
            v = v[np.isfinite(v)]
            enough = len(v) and (len(g) - len(v)) <= max_missing * len(g)
            low, high = percentile_interval(v) if enough else (np.nan, np.nan)
            rows.append({"era": era, "key": key, "low": low, "high": high, "draws": len(v), "of": len(g)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Survey-design bootstrap (Results 5)
# ---------------------------------------------------------------------------
# -- UNCERTAINTY: the survey-design bootstrap (Results 5) --------------------
# Methods, "Uncertainty": "... or primary sampling units within strata [53]".
# How (Rao-Wu): in each stratum with n sampling units, draw n - 1 units with
# replacement; each participant's weight is n / (n - 1) x the number of times
# their unit was drawn. Results 5 uses 2,000 draws. (Results 6 has its own
# version, s6_screen.replicate_weights.)
def survey_bootstrap_weights(frame: pd.DataFrame, stream: int, draws: int = 2000, seed: int = 20261011):
    """Rao-Wu rescaled bootstrap: replicate weights that respect the survey design.

    In each sampling stratum with n primary sampling units, n - 1 units are
    drawn with replacement; a participant's weight is n / (n - 1) times the
    number of times their unit was drawn. Built on the whole panel (`frame`
    needs `stratum` and `psu`), so subgroups formed later share the draws.

    Returns (unit_weight, unit_of): unit_weight[draw, unit], and each frame
    row's unit. A row's weight in draw b is unit_weight[b, unit_of[row]].
    """
    rng = np.random.default_rng([seed, stream])
    names, blocks = [], []
    for _, units in frame.groupby("stratum").psu.unique().sort_index().items():
        units = sorted(units)
        n = len(units)
        drawn = rng.integers(0, n, size=(draws, n - 1))
        blocks.append(np.stack([(drawn == j).sum(axis=1) for j in range(n)], axis=1) * (n / (n - 1.0)))
        names.extend(units)
    unit_of = pd.Series(np.arange(len(names)), index=names).loc[frame.psu].to_numpy()
    return np.concatenate(blocks, axis=1), unit_of


# -- CALCULATION: concordance within sex and period --------------------------
# Methods, "Held-out prediction": "... and by risk-set concordance [51]."
# How: for each death at age t, compare the person who died with every
# same-sex, same-period participant still at risk at t; the pair is
# concordant if the person who died had the higher score (ties count one
# half). Used for the concordance results of Results 5.
def concordance_within(entry, exit_, event, sex, strata, score, weights=None) -> float:
    """Concordance comparing participants of the same sex and stratum (survey cycle or phase).

    For each death at age t, the comparators are same-sex, same-stratum
    participants who had entered by t and were still under observation after
    t; the pair is concordant if the one who died had the higher score (ties
    count one half). With weights, each pair counts w_death x w_comparator.
    """
    w = np.ones(len(exit_)) if weights is None else np.asarray(weights, float)
    numerator = denominator = 0.0
    for g in np.unique(strata):
        for s in np.unique(sex):
            ix = np.flatnonzero((strata == g) & (sex == s) & (w > 0))
            en, ex, d, v, ww = entry[ix], exit_[ix], event[ix], score[ix], w[ix]
            dead = np.flatnonzero(d == 1)
            comparable = (en[None, :] <= ex[dead, None]) & (ex[None, :] > ex[dead, None])
            pair = ww[dead, None] * ww[None, :] * comparable
            gap = v[dead, None] - v[None, :]
            numerator += float((pair * ((gap > 0) + .5 * (gap == 0))).sum())
            denominator += float(pair.sum())
    return numerator / denominator


# -- CALCULATION: percentiles with weights (Results 6) -----------------------
# Used for the screen's young-healthy medians, the 99th-percentile caps of the
# bounded distance and the within-sex 60th-percentile cut-offs, so that the
# same rule works on the survey-bootstrap weights. At equal weights it equals numpy's ordinary
# percentile.
def weighted_quantile(x: np.ndarray, w: np.ndarray, p: float) -> float:
    """Quantile with case weights; equals numpy's linear-interpolation quantile at unit weights.

    Values sorted, cumulative weights C, total N, position h = 1 + (N - 1) p;
    interpolate between the first values whose cumulative weight reaches
    floor(h) and ceil(h). Rows with zero weight are ignored.
    """
    keep = w > 0
    x, w = np.asarray(x, float)[keep], np.asarray(w, float)[keep]
    order = np.argsort(x, kind="mergesort")
    xs, C = x[order], np.cumsum(w[order])
    h = 1.0 + (C[-1] - 1.0) * p
    lo, hi = np.floor(h), np.ceil(h)

    def first_reaching(rank):
        reached = np.flatnonzero(C >= rank - 1e-9)
        return xs[reached[0]] if len(reached) else xs[-1]
    low, high = first_reaching(lo), first_reaching(hi)
    return float(low + (h - lo) * (high - low))


# -- CALCULATION: five-year risk of death (Results 6) ------------------------
# Methods, "Screening risks": "Kaplan-Meier estimates [52]."
# How: 1 - Kaplan-Meier survival at 60 months, optionally weighted.
def kaplan_meier_risk(months: np.ndarray, died: np.ndarray, weights: np.ndarray | None = None,
                      horizon: float = 60) -> float:
    """Probability of death within `horizon` months (1 - weighted Kaplan-Meier survival).

    Deaths at month 0 count; a death and a censoring at the same month are
    handled as simultaneous (the censored person is still at risk).
    """
    w = np.ones(len(months)) if weights is None else np.asarray(weights, float)
    keep = w > 0
    t, d, w = months[keep], died[keep], w[keep]
    survival = 1.0
    for time in np.unique(t[(d == 1) & (t <= horizon)]):
        survival *= max(0.0, 1 - w[(t == time) & (d == 1)].sum() / w[t >= time].sum())
    return 1 - survival
