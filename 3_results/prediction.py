# =============================================================================
# READ ME FIRST: prediction.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Held-out prediction in the matched NHANES 1999-2010 cohort:
#               learn a score on five survey cycles, predict the sixth, repeat
#               for all six, and pool the held-out predictions. Reports
#               concordance (C), log likelihood and five-year risk error.
#
# METHODS       "Statistical analysis": "Held-out prediction" ("evaluated by the
#               evaluation-set Cox log partial likelihood at training
#               coefficients ... and by risk-set concordance [51]. Both assess
#               relative-hazard prediction, not calibrated absolute risk.")
#
# USED BY       Results 2, Table 2 (s2_reference.py): is a sex term still needed?
#               Results 4, the matched Phenotypic Age comparison (s4_panel.py).
#
# PSEUDOCODE
#   for each held-out cycle (six cycles, 1999-2010):
#     train = the other five cycles
#     q weights = Cox(markers [+ female]) on train, rescaled to length one
#     evaluation model = Cox(age [+ female] + q) on train
#     held-out linear predictor and five-year risk for the held-out cycle
#   pool the six held-out cycles; then, for all, men and women:
#     C          = share of comparable pairs (same cycle) ranked correctly
#     risk error = mean predicted five-year risk - Kaplan-Meier observed risk
# =============================================================================
"""Held-out prediction: does a score predict death in survey cycles it never saw?

Used for Table 2 (does an explicit sex term improve prediction?) and the
PhenoAge comparison. The design is the same throughout:

  1. Hold out one survey cycle of the matched NHANES 1999-2010 cohort.
  2. On the other five cycles, learn the score's weights and an "evaluation"
     Cox model of the endpoint on age (+ sex) + the score. Table 2 (s2) also
     re-learns the young-healthy reference from those five cycles' reference
     members; the Phenotypic Age comparison (s4) uses coordinates from the
     prepared full-cohort reference, as the reviewed comparison did.
  3. Predict the held-out cycle: a linear predictor (for ranking) and a
     five-year risk of death (for calibration).
  4. Repeat for each of the six cycles, then pool the held-out predictions.

Metrics, all computed on held-out participants only:
  C            concordance: of comparable pairs (one died, the other was still
               alive and under observation at that age, same cycle), the share
               in which the one who died had the higher predicted risk.
  risk error   predicted minus observed five-year risk, in percentage points;
               predicted risk uses the training baseline hazard, observed risk
               is the Kaplan-Meier estimate. This checks calibration in the large,
               not risk in individuals.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg

import models


def training_risk_sets(people: pd.DataFrame) -> models.RiskSets:
    """The Cox risk sets of the training rows (attained age, cycle strata), built once per fold."""
    return models.RiskSets(people.entry_age.to_numpy(float), people.exit_age.to_numpy(float),
                           people.event.to_numpy(int), people.cycle.to_numpy(str))


def fit(risk: models.RiskSets, X: np.ndarray) -> np.ndarray:
    """Unpenalised Cox coefficients for columns X (Efron ties, delayed entry, cycle strata).

    Uses the vectorised solver in models.py rather than statsmodels: the same
    likelihood, about twenty times faster, which matters for 500 bootstrap draws.
    """
    beta, converged = models.ridge_cox(risk, X, np.zeros(X.shape[1]))
    if not converged:
        raise ArithmeticError("Cox fit did not converge")
    return beta


# -- CALCULATION: the q weights, learned on training cycles only -------------
# Results 2: "omitting sex both when the weights were learned and when the
#   q score was evaluated (Table 2)". with_sex = False drops the female column
#   from this fit; with_sex = True fits it but leaves it out of the score.
def learn_weights(people, X, train, risk, with_sex: bool) -> np.ndarray:
    """The mortality score q learned on training rows: unit-length marker weights.

    A sex term is fitted alongside the markers when `with_sex` is true; the
    score itself never contains it.
    """
    design = np.column_stack(([people.female.to_numpy()] if with_sex else []) + [X])
    beta = fit(risk, design[train])
    weights = beta[1:] if with_sex else beta
    return X @ (weights / np.linalg.norm(weights))


# -- CALCULATION: the evaluation model ---------------------------------------
# Cox model of the endpoint on age (+ female when with_sex) + the score, fitted
# on training rows; returns the linear predictor for every row. Table 2's
# "C, no sex term" and "C, with sex term" columns compare with_sex False/True.
# This model includes examination age as a covariate as well as the
# attained-age time scale. Table 2 legend: "Evaluation models were Cox models on
# attained age, stratified by survey cycle, containing examination age, the
# score and, where stated, a sex term."
def evaluation_model(people, score, train, risk, with_sex: bool) -> np.ndarray:
    """Linear predictor from a Cox model on age (+ sex) (+ the score, in training-SD units)."""
    terms = ([people.female.to_numpy()] if with_sex else []) + [people.RIDAGEYR.to_numpy()]
    if score is not None:
        terms.append(score / score[train].std(ddof=1))
    design = np.column_stack(terms)
    return design @ fit(risk, design[train])


def baseline_hazards(people: pd.DataFrame, eta: np.ndarray):
    """Cumulative baseline hazard per cycle (Efron/Breslow increments), from training rows."""
    out = []
    a, b = people.entry_age.to_numpy(), people.exit_age.to_numpy()
    event = people.event.to_numpy().astype(bool)
    cycle = people.cycle.to_numpy(str)
    for c in sorted(np.unique(cycle)):
        rows = cycle == c
        ac, bc, ec, w_eta = a[rows], b[rows], event[rows], eta[rows]
        shift = float(w_eta.max())
        w = np.exp(w_eta - shift)
        times = np.unique(bc[ec])
        increments = []
        for t in times:
            dead = (bc == t) & ec
            at_risk = (ac <= t) & (bc >= t)
            n = int(dead.sum())
            # Efron: with n tied deaths, the risk set shrinks by 1/n of the dying at each step.
            denominators = w[at_risk].sum() - np.arange(n) / n * w[dead].sum()
            increments.append(float((1 / denominators).sum()))
        out.append((times, np.r_[0., np.cumsum(increments)], shift, float(bc.max())))
    return out


# -- CALCULATION: predicted five-year risk -----------------------------------
# Each training cycle gives a baseline survival curve; predicted risk = 1 -
# average five-year survival from the participant's examination age. Results 2:
# "without a sex term the q score over-predicted average five-year risk by
# similar amounts in women and men (0.48 and 0.45 percentage points)".
def five_year_risk(train_people, train_eta, test_people, test_eta) -> np.ndarray:
    """Probability of an endpoint death within five years of examination.

    Each training cycle gives a baseline survival curve; the predicted risk is
    one minus the average, over training cycles, of each curve's five-year
    survival from the participant's age at examination.
    """
    a = test_people.entry_age.to_numpy()
    survival = []
    for times, cumulative, shift, last in baseline_hazards(train_people, train_eta):
        if (a + 5 > last).any():
            raise ArithmeticError("five-year horizon beyond a training cycle's follow-up")
        increase = cumulative[np.searchsorted(times, a + 5, side="right")] - \
            cumulative[np.searchsorted(times, a, side="right")]
        survival.append(np.exp(-np.exp(test_eta - shift) * increase))
    return 1 - np.mean(survival, axis=0)


# -- CALCULATION: concordance (C) --------------------------------------------
# Methods, "Held-out prediction": "risk-set concordance [51]". For each death
# at age t, every same-cycle participant still at risk at t is a comparator;
# the pair is concordant if the person who died had the higher predicted risk
# (ties count one half). Results 2: "change in concordance 0.0002".
def concordance_counts(people: pd.DataFrame, eta: np.ndarray) -> tuple[float, float]:
    """(concordant pairs, comparable pairs); ties in predicted risk count one half."""
    a, t = people.entry_age.to_numpy(), people.exit_age.to_numpy()
    deaths = np.flatnonzero(people.event.to_numpy())
    g = people.cycle.to_numpy(str)
    comparable = (a[None, :] <= t[deaths, None]) & (t[None, :] > t[deaths, None]) & (g[None, :] == g[deaths, None])
    gap = eta[deaths, None] - eta[None, :]
    return float((comparable & (gap > 0)).sum() + .5 * (comparable & (gap == 0)).sum()), float(comparable.sum())


# -- CALCULATION: held-out log likelihood ------------------------------------
# Methods, "Held-out prediction": "evaluation-set Cox log partial likelihood at
# training coefficients". The predictor is fixed; nothing is fitted here.
def log_likelihood(people: pd.DataFrame, eta: np.ndarray) -> float:
    """Cox partial log-likelihood of a fixed predictor (Efron ties, cycle strata)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")        # evaluates a fixed predictor; no fit happens here
        model = PHReg(people.exit_age.to_numpy(float), np.asarray(eta)[:, None], status=people.event.to_numpy(int),
                      entry=people.entry_age.to_numpy(float), strata=people.cycle.to_numpy(str), ties="efron")
    return float(model.loglike(np.ones(1)))


def kaplan_meier_risk(years: np.ndarray, event: np.ndarray, horizon: float = 5) -> float:
    """Observed probability of an endpoint death within `horizon` years (1 - Kaplan-Meier survival)."""
    survival = 1.
    for t in np.unique(years[(years <= horizon) & event]):
        survival *= 1 - event[years == t].sum() / (years >= t).sum()
    return 1 - survival


GROUPS = {"all": None, "men": 0, "women": 1}


# -- CALCULATION: pooling the six held-out cycles ----------------------------
# Sums concordant and comparable pairs, log likelihoods and predicted risks
# over the six held-out cycles, for all participants, men and women.
# risk_error = mean predicted risk - observed Kaplan-Meier risk (Table 2's
# "Five-year prediction error" columns, in percentage points after x 100).
def evaluate(people: pd.DataFrame, held_out: dict) -> dict:
    """Pool one model's held-out predictions over the six folds.

    `held_out` maps each cycle to (eta, risk) for that cycle's participants.
    """
    out = {}
    for group, female in GROUPS.items():
        wins = pairs = ll = predicted = 0.
        n = 0
        for cycle, (eta, risk) in held_out.items():
            rows = people.cycle.eq(cycle).to_numpy()
            if female is not None:
                keep = people.female.to_numpy()[rows] == female
            else:
                keep = np.ones(rows.sum(), bool)
            sub = people[rows][keep].reset_index(drop=True)
            w, p = concordance_counts(sub, eta[keep])
            wins, pairs = wins + w, pairs + p
            ll += log_likelihood(sub, eta[keep])
            predicted += risk[keep].sum()
            n += len(sub)
        mask = np.ones(len(people), bool) if female is None else people.female.eq(female).to_numpy()
        observed = kaplan_meier_risk(people.followup_years.to_numpy()[mask], people.event.to_numpy()[mask].astype(bool))
        out[group] = {"C": wins / pairs, "ll": ll, "predicted_risk": predicted / n,
                      "observed_risk": observed, "risk_error": predicted / n - observed}
    return out
