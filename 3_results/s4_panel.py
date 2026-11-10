# =============================================================================
# READ ME FIRST: s4_panel.py  (Results 4)
# -----------------------------------------------------------------------------
# RESULTS 4     "Ten routine markers retained prediction and outperformed
#               Phenotypic Age."
#
# WHAT IT DOES  (1) In each training sample, ranks the markers by elastic-net
#               Cox, refits the top 3-15 as a reduced q score, and tests it on
#               held-out data against the full q score (eight adequacy
#               criteria). (2) In the matched NHANES 1999-2010 cohort, compares
#               the full and ten-marker q scores with Phenotypic Age.
#
# METHODS       Supplementary Methods, "Selection and refitting"; "Held-out
#               prediction"; "Survival models".
#
# MAKES         Figure 4 source tables (a-d); Table 3; Supp. Tables S4.2-S4.4,
#               S4.7; every Results 4 number, except the 32-marker q comparison
#               ("difference 0.0003"; Supp. Tables S4.5-S4.6), which is not in
#               this release.
#
# NOT IN THIS RELEASE  The inner search that chose each selection's penalty and
#               L1 fraction (about 177 CPU-hours). This script reads the choices
#               from inputs/panel_selection_tuning.csv and repeats every step
#               after them.
#
# PSEUDOCODE
#   (1) for each evaluation (252 five-vs-five cycle divisions; 10 one-cycle
#       holdouts; 200 NHANES III half-splits; 10 NHANES III folds):
#     rank markers on TRAINING rows: elastic-net Cox, sex unpenalised, stored
#       penalty; rank by |coefficient|
#     full q  = Cox(female + all markers) on training; direction u = beta/|beta|
#     for size in 3, 5, 7, 10, 12, 15:
#       reduced q = Cox(female + top markers) on training; its direction
#       on HELD-OUT rows: correlation with full q, gain over sex alone,
#         retained gain = reduced gain / full gain, HR in men and in women,
#         sex contrast
#     adequacy: the eight criteria (Supp. Table S4.1) over all evaluations
#   locked ten = the same ranking on all ten cycles (Table 3)
#   (2) for each of the six matched cycles held out:
#     q38, q10 (ten chosen without that cycle), local refit on Phenotypic Age
#       inputs, published Phenotypic Age
#     evaluation Cox(female + age + score) on training; C on the held-out cycle
#   intervals for (2): 500 participant draws, the ten markers held fixed
# =============================================================================
"""Results 4: ten routine markers retained prediction and outperformed Phenotypic Age.

The questions
-------------
(1) How far can the 38-marker panel be reduced while keeping the full score's
behaviour? (2) How does a reduced score compare with Phenotypic Age?

How (1): nested selection
-------------------------
In each training sample, an elastic-net Cox model (female unpenalised; one
penalty and L1 fraction, see below) ranks the markers by the size of their
coefficients. The top 3, 5, 7, 10, 12 or 15 are refitted with an ordinary Cox
model and evaluated on held-out data. Training samples: all 126 divisions of
the ten survey cycles into five and five, each used in both directions (252
evaluations); and, as a sensitivity, one cycle held out at a time (10).
A size is adequate if it meets eight criteria together (Supplementary Table S4.1):
  1-2  the reduced score correlates with the full score in held-out data:
       median >= 0.90, 2.5th percentile >= 0.85;
  3-4  the same, when the selected markers are scored in adverse-excursion
       coordinates: median >= 0.80, 2.5th percentile >= 0.70;
  5    it retains a median >= 0.90 of the full score's held-out gain over sex alone;
  6    its log hazard ratio is positive at the 2.5th percentile in men and in women;
  7    women score lower than men (age- and cycle-adjusted) in >= 95% of evaluations;
  8    it keeps a median >= 0.80 of the full score's sex contrast.
The ten markers of Table 3 come from the same procedure applied once to all
ten cycles.

NHANES III repeats the selection on its own data (34 markers): 100 random
half-splits of the sampling units (one unit of each stratum to each half, seed
20261017), each used in both directions (200 evaluations), and, as a
sensitivity, ten folds of whole sampling units held out one at a time.

The penalty choice. Each selection's penalty and L1 fraction were chosen by
an inner leave-one-cycle-out search over a 7 x 4 grid. That search took about
177 CPU-hours, so by default this script reads the choices from
inputs/panel_selection_tuning.csv and repeats every step after them.

How (2): matched comparison with Phenotypic Age
-----------------------------------------------
In the matched NHANES 1999-2010 cohort, each survey cycle is predicted from
the other five. Scores: the full 38-marker score; a ten-marker score (markers
selected without the held-out cycle); a Cox model refitted locally on
Phenotypic Age's nine blood inputs plus age; and published Phenotypic Age
(Levine et al., 2018). Every evaluation model includes sex and age.
Uncertainty (R064): 500 paired bootstrap draws (numbered from 0) resample the
full 1999-2018 cohort within cycle x sex and keep its matched members; every
weight fit and evaluation is repeated, while each fold's ten selected markers
are held fixed. Intervals are the 2.5th and 97.5th percentiles (numpy's linear
rule), so they do not include the uncertainty of the marker selection.

Scope
-----
Population   the two analysis cohorts; the comparison uses the matched NHANES
             1999-2010 cohort (6,759 participants, 2,112 deaths).
Endpoint     accident-censored mortality.
Estimands    for each panel size, held-out agreement with the full score, the share
             of its held-out gain retained, sex-specific hazard ratios and the sex
             contrast, summarised over evaluations; the smallest size meeting all
             eight criteria; the locked ten. For the comparison: held-out
             concordance of each score and paired differences.
Split        every selection, refit and threshold uses training rows only; metrics
             are computed on the held-out rows.
Uncertainty  the bars of Figure 4a,b are percentiles across evaluations (how much the
             answer varies with the training sample), not confidence intervals. The
             comparison's intervals come from its own bootstrap (above).
Assumptions  the stored penalty choices stand for the inner search (inputs/); the
             transformations and references are fixed.
Permitted    whether ten markers keep the full score's held-out behaviour under the
interpretation  stated criteria, and how a ten-marker score ranks deaths beside
             Phenotypic Age in the same participants. Not a claim that these ten are
             the only adequate set (Table 3 shows how often each was chosen).

What it makes (in outputs/s4/)
------------------------------
  evaluations.csv          every held-out evaluation (selection, size, coordinates) and its metrics
  adequacy.csv             the eight criteria by panel size (full run)
  locked_ten_<era>.csv     the locked ten markers and how often each was selected at size ten (full run)
  table3.csv               Table 3 (full run)
  phenoage.csv             held-out concordance and paired differences, Fig. 4c,d, with intervals
  benchmark_bootstrap.csv  the 500 benchmark draws (with --bootstrap) and its manifest
  numbers.csv              every Results 4 number quoted in the text

    python -B s4_panel.py                      # all 472 evaluations (~30 min on 12 cores); reuses stored draws
    python -B s4_panel.py --bootstrap          # also the 500 benchmark draws
    python -B s4_panel.py --evaluations 2 --only F09 --check-draws 2
                                               # a check: two evaluations of each primary design, the
                                               # NHANES III fold F09, and benchmark draws 0-1
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
from lifelines import CoxPHFitter         # noqa: E402
from statsmodels.duration.hazard_regression import PHReg   # noqa: E402

import data                               # noqa: E402
import models                             # noqa: E402
import prediction as P                    # noqa: E402
import runs                               # noqa: E402
import settings as S                      # noqa: E402

OUT = S.output_folder("s4")
SIZES = (3, 5, 7, 10, 12, 15)
# Validity limits from the reviewed selection specifications (modern selection_spec.csv and
# NHANES_III/config/panel_selection_spec.csv; identical in both eras). A ratio is only meaningful
# when its denominator is clearly positive: if the full 38-marker score predicts no better than
# sex alone in a held-out sample (full gain <= 1e-10), "the share of that gain kept by fewer
# markers" is undefined, and the reviewed analysis records it as missing rather than as a number.
# The same applies to the sex contrast, and a direction with norm <= 1e-12 has no direction at all.
NEAR_ZERO_GAIN = 1e-10
NEAR_ZERO_NORM = 1e-12


# A ratio such as "retained gain" is left missing when the full score has no
# held-out gain to share (see the comment above).
def _ratio(numerator: float, denominator: float) -> float:
    """numerator / denominator, or missing when the denominator is not above NEAR_ZERO_GAIN.

    Medians over evaluations then skip the missing ones, as the reviewed summaries did
    (one NHANES III held-out-fold evaluation at size ten has a negative full gain).
    """
    return numerator / denominator if denominator > NEAR_ZERO_GAIN else np.nan


def _unit(beta: np.ndarray) -> np.ndarray:
    """A fitted marker direction scaled to length one; stops if the fit has no direction."""
    norm = np.linalg.norm(beta)
    if not np.isfinite(norm) or norm <= NEAR_ZERO_NORM:
        raise ValueError(f"fitted marker direction has norm {norm}; the reviewed analysis treats this as a failed fit")
    return beta / norm
TUNING = S.PIPELINE / "inputs" / "panel_selection_tuning.csv"
# Published Phenotypic Age (Levine et al., 2018): weights for albumin, creatinine, glucose, log CRP,
# lymphocyte %, mean cell volume, red cell distribution width, alkaline phosphatase, white cells, age.
BENCHMARK_SEED, BENCHMARK_DRAWS = 20261013, 500      # R064's stream: SeedSequence([seed, draw])
PHENOAGE_WEIGHTS = np.array([-.03359355, .009506491, .1953192, .09536762, -.01199984, .02676401, .3306156,
                             .001868778, .05542406, .08035356])


# ===========================================================================
# Fitting (each era keeps its reviewed conventions)
# ===========================================================================
# NHANES 1999-2018: statsmodels Cox stratified by survey cycle, held-out likelihood with
# inclusive entry. NHANES III: lifelines, unstratified, held-out likelihood with strict
# entry (someone entering at exactly a death age is not yet at risk), as its selection did.
def frame_for(people: pd.DataFrame, values: pd.DataFrame, markers: list) -> pd.DataFrame:
    cols = ["entry_age", "exit_age", "event", "female"] + (["cycle"] if "cycle" in people else [])
    f = people[cols].reset_index(drop=True).copy()
    f[markers] = values[markers].to_numpy(float)
    return f


def era_of(frame: pd.DataFrame) -> str:
    return "modern" if "cycle" in frame else "iii"


def cox(frame: pd.DataFrame, terms: list) -> np.ndarray:
    fit = models.cox(era_of(frame), frame, terms)
    return np.array([fit[t] for t in terms])


# -- CALCULATION: held-out log likelihood at training coefficients -----------
# Methods, "Held-out prediction": "evaluated by the evaluation-set Cox log
#   partial likelihood at training coefficients".
def heldout_loglike(frame: pd.DataFrame, terms: list, beta: np.ndarray) -> float:
    """Cox partial log-likelihood of fixed coefficients on new participants."""
    if era_of(frame) == "modern":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")      # building the model only; no fit happens here
            model = PHReg(frame.exit_age.to_numpy(float), frame[terms].to_numpy(float), status=frame.event.to_numpy(int),
                          entry=frame.entry_age.to_numpy(float), strata=frame.cycle.to_numpy(str), ties="efron")
        return float(model.loglike(beta))
    entry = np.nextafter(frame.entry_age.to_numpy(float), np.inf)     # strict entry
    risk = models.RiskSets(entry, frame.exit_age.to_numpy(float), frame.event.to_numpy(int), np.zeros(len(frame), str))
    return risk.evaluate(frame[terms].to_numpy(float), beta, order=0)


# -- CALCULATION: ranking the markers ----------------------------------------
# Supplementary Methods, "Selection and refitting": "Within each training
#   sample, an elastic-net attained-age Cox model [59, 60] ranked markers by
#   absolute penalised coefficient, with sex unpenalised."
# How: lifelines CoxPHFitter with penalizer 0 for female and the stored
#   penalty for every marker, the stored L1 fraction, attained age with entry,
#   cycle strata in NHANES 1999-2018; markers sorted by |coefficient| (ties by
#   name). Only training rows are passed in.
def elastic_net_ranking(frame: pd.DataFrame, markers: list, penalty: float, l1_fraction: float) -> list:
    """Markers ranked by the size of their elastic-net Cox coefficient (ties broken by name).

    lifelines, female unpenalised (stratified by cycle for 1999-2018), started from
    the unpenalised fit, at most 100 steps to precision 1e-6 (the reviewed settings).
    """
    seed = cox(frame, ["female"] + markers)
    strata = {"strata": ["cycle"]} if era_of(frame) == "modern" else {}
    columns = ["entry_age", "exit_age", "event", "female"] + (["cycle"] if strata else []) + markers
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = CoxPHFitter(penalizer=np.r_[0.0, np.repeat(penalty, len(markers))], l1_ratio=l1_fraction).fit(
            frame[columns], duration_col="exit_age", event_col="event", entry_col="entry_age", initial_point=seed,
            show_progress=False, fit_options={"max_steps": 100, "precision": 1e-6}, **strata)
    beta = fit.params_.loc[markers].to_numpy(float)
    models.note_fit("elastic_net", caught, era=era_of(frame), rows=len(frame), penalty=penalty,
                    l1_fraction=l1_fraction)
    models.require_finite(beta, "elastic-net ranking")
    return [markers[j] for j in sorted(range(len(markers)), key=lambda j: (-abs(beta[j]), markers[j]))]


# -- CALCULATION: sex contrast of a score (adequacy criteria 7 and 8) --------
# Women minus men in the score, adjusted for an age spline and survey period
# (the same model as models.sex_contrast).
def _ols_female(people, score):
    """Women minus men in a (training-standardised) score, adjusted for age and survey period: OLS."""
    from patsy import dmatrix
    age, period = ("RIDAGEYR", "cycle") if "cycle" in people else ("age", "phase")
    X = dmatrix(f"female + bs({age}, df=4, degree=3, include_intercept=False) + C({period})",
                people[["female", age, period]], return_type="dataframe")
    return float((np.linalg.pinv(X.to_numpy(float)) @ np.asarray(score, float))[list(X.columns).index("female")])


# ===========================================================================
# One held-out evaluation
# ===========================================================================
_D = {}


# -- DESIGN: NHANES III training and evaluation sampling units ---------------
# Results 4: "100 sampling-unit half-splits in NHANES III, also evaluated in
#   both directions (200 evaluations)". One unit of each stratum goes to each
#   half (seed 20261017); the ten grouped folds are the one-fold sensitivity.
def iii_splits() -> dict:
    """NHANES III training/test sampling units: 100 seeded half-splits (one unit of each
    stratum to each half; both directions) and the ten grouped folds."""
    from sklearn.model_selection import GroupKFold
    people = _D["iii"]["people"]
    strata = sorted(people.SDPSTRA6.astype(int).unique())
    rng, seen, splits = np.random.default_rng(20261017), set(), {}
    while len(seen) < 100:
        bits = tuple(int(b) for b in rng.integers(0, 2, size=len(strata)))
        if bits in seen or tuple(1 - b for b in bits) in seen:
            continue
        seen.add(bits)
        k = len(seen)
        a = {f"{s}_{1 if b else 2}" for s, b in zip(strata, bits)}
        b_ = {f"{s}_{2 if b else 1}" for s, b in zip(strata, bits)}
        splits[f"P{k:03d}_AtoB"], splits[f"P{k:03d}_BtoA"] = (a, b_), (b_, a)
    groups = _D["iii"]["group"]
    for f, (_, test) in enumerate(GroupKFold(10).split(people, groups=groups), start=1):
        held = set(groups[test])
        splits[f"F{f:02d}"] = (set(groups) - held, held)
    return splits


def _start() -> None:
    c = data.load_modern()
    _D["modern"] = {"people": c.people, "coords": {"HD": c.z, "FI": c.x}, "markers": c.markers}
    c3 = data.load_iii()
    _D["iii"] = {"people": c3.people, "coords": {"HD": c3.z, "FI": c3.x}, "markers": c3.markers,
                 "group": (c3.people.SDPSTRA6.astype(int).astype(str) + "_" + c3.people.SDPPSU6.astype(int).astype(str)).to_numpy()}
    _D["iii_splits"] = iii_splits()


# -- CALCULATION: one held-out evaluation of every panel size ----------------
# Results 4: "Within each training sample, we ranked the markers with an
#   elastic-net Cox model, refitted the top 3, 5, 7, 10, 12 or 15, and tested
#   them on held-out data." Supplementary Methods: "The top 3, 5, 7, 10, 12 or
#   15 markers were then refitted by unpenalised Cox regression with sex, and
#   the training coefficients and score scale were applied unchanged to
#   evaluation data."
# How: ranking, the full fit and every reduced fit use training rows only;
#   the held-out score is standardised with the training mean and SD
#   (c_tr.mean(), c_tr.std()); the retained gain is (reduced held-out gain over
#   sex alone) / (full held-out gain over sex alone), Figure 4a,b.
#   "HD" = z coordinates; "FI" = adverse-excursion coordinates (criteria 3-4).
def evaluate(job: dict) -> list:
    """Rank on the training sample, then score every panel size in both coordinate systems."""
    era = "iii" if job["analysis"].startswith("iii") else "modern"
    d = _D[era]
    people, coords, markers = d["people"], d["coords"], d["markers"]
    if era == "modern":
        train = people.cycle.isin(job["training_cycles"].split("|")).to_numpy()
        test = people.cycle.isin(job["test_cycles"].split("|")).to_numpy()
    else:
        a, b = _D["iii_splits"][job["outer_id"]]
        train, test = np.isin(d["group"], list(a)), np.isin(d["group"], list(b))
    ranking = elastic_net_ranking(frame_for(people[train], coords["HD"][train], markers), markers,
                                  job["penalty"], job["l1_fraction"])
    sex_tr, sex_te = frame_for(people[train], coords["HD"][train], []), frame_for(people[test], coords["HD"][test], [])
    sex_only = heldout_loglike(sex_te, ["female"], cox(sex_tr, ["female"]))
    rows = []
    for rep in ("HD", "FI"):
        tr, te = frame_for(people[train], coords[rep][train], markers), frame_for(people[test], coords[rep][test], markers)
        full = cox(tr, ["female"] + markers)
        full_ll = heldout_loglike(te, ["female"] + markers, full)
        u_full = _unit(full[1:])
        q_tr, q_te = tr[markers].to_numpy() @ u_full, te[markers].to_numpy() @ u_full
        full_z = (q_te - q_tr.mean()) / q_tr.std(ddof=1)
        full_contrast = _ols_female(people[test], full_z)
        for size in SIZES:
            chosen = ranking[:size]
            cand = cox(tr, ["female"] + chosen)
            u = _unit(cand[1:])
            c_tr, c_te = tr[chosen].to_numpy() @ u, te[chosen].to_numpy() @ u
            z = (c_te - c_tr.mean()) / c_tr.std(ddof=1)
            cand_ll = heldout_loglike(te, ["female"] + chosen, cand)
            row = {"analysis": job["analysis"], "outer_id": job["outer_id"], "coordinates": rep, "panel_size": size,
                   "selected": "|".join(chosen), "correlation": float(np.corrcoef(z, full_z)[0, 1]),
                   "full_gain": full_ll - sex_only, "candidate_gain": cand_ll - sex_only,
                   "retained_gain": _ratio(cand_ll - sex_only, full_ll - sex_only)}
            for label, female in (("men", 0), ("women", 1)):
                rows_sex = te.female.eq(female).to_numpy()
                f = te[rows_sex].assign(score=z[rows_sex])
                row[f"log_hr_{label}"] = cox(f, ["score"])[0]
            contrast = _ols_female(people[test], z)
            row.update(sex_contrast=contrast, sex_contrast_retained=_ratio(abs(contrast), abs(full_contrast)))
            rows.append(row)
    return rows


# -- CALCULATION: the eight adequacy criteria (Supp. Table S4.1) -------------
# Results 4: "We counted a panel size as adequate only if it met all eight
#   prespecified criteria: the reduced q score had to reproduce the full q
#   score, retain a median of at least 90% of its held-out gain, predict
#   mortality in both sexes and preserve the full q score's sex contrast".
#   "Ten was the smallest tested size that met all eight criteria in both eras".
# How: each criterion is a median or 2.5th percentile across evaluations,
#   compared with its threshold (listed below, and in the module docstring).
def adequacy(ev: pd.DataFrame) -> pd.DataFrame:
    """The eight criteria by panel size (HD selection; 'cross' = scored in adverse-excursion coordinates)."""
    out = []
    for size in SIZES:
        w = ev[ev.coordinates.eq("HD") & ev.panel_size.eq(size)]
        x = ev[ev.coordinates.eq("FI") & ev.panel_size.eq(size)]
        values = {"within_correlation_median": (w.correlation.median(), .90),
                  "within_correlation_p2_5": (w.correlation.quantile(.025), .85),
                  "cross_correlation_median": (x.correlation.median(), .80),
                  "cross_correlation_p2_5": (x.correlation.quantile(.025), .70),
                  "retained_gain_median": (w.retained_gain.median(), .90),
                  "men_log_hr_p2_5": (w.log_hr_men.quantile(.025), 0.0),
                  "women_log_hr_p2_5": (w.log_hr_women.quantile(.025), 0.0),
                  "sex_contrast_negative_fraction": (w.sex_contrast.lt(0).mean(), .95),
                  "sex_contrast_retained_median": (w.sex_contrast_retained.median(), .80)}
        for name, (value, limit) in values.items():
            passed = value > limit if "log_hr" in name else value >= limit
            out.append({"panel_size": size, "criterion": name, "value": value, "threshold": limit, "passed": passed})
    return pd.DataFrame(out)


# ===========================================================================
# (2) The matched comparison with Phenotypic Age
# ===========================================================================
# -- PHENOTYPIC AGE COMPARISON: set-up ---------------------------------------
# Results 4: "a matched NHANES 1999-2010 cohort (6,759 participants; 2,112
#   deaths), predicting each of six survey cycles from the other five".
# How: published Phenotypic Age from its published formula (PHENOAGE_WEIGHTS,
#   Levine et al., 2018). For each held-out cycle, the ten markers are ranked on
#   the other nine cycles of the full 1999-2018 cohort, so the held-out cycle
#   never contributes to its own ten.
def benchmark_setup(tuning: pd.DataFrame) -> dict:
    """The matched cohort, its inputs, published Phenotypic Age and each fold's ten markers.

    The ten markers of each fold are selected once, on the other nine cycles of the full
    1999-2018 cohort (stored penalty choice), and held fixed in every bootstrap draw.
    """
    c = data.load_modern()
    keep = data.matched_rows(c)
    people = c.people[keep].reset_index(drop=True)
    z = c.z[keep].reset_index(drop=True)
    inputs = data.phenoage_inputs(people.SEQN).assign(age=people.RIDAGEYR.to_numpy()).to_numpy(float)
    xb = -19.90667 + inputs @ PHENOAGE_WEIGHTS
    published = 141.50225 + (np.log(.0055305 * 1.51714 / .007692696) + xb) / .090165
    ten = {}
    for r in tuning[tuning.analysis.eq("phenoage_fold")].itertuples():
        sel = c.people.cycle.isin(r.training_cycles.split("|")).to_numpy()
        ten[r.test_cycles] = elastic_net_ranking(frame_for(c.people[sel], c.z[sel], c.markers), c.markers,
                                                 r.penalty, r.l1_fraction)[:10]
    matched_position = np.full(len(c.people), -1)
    matched_position[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    return {"people": people, "z": z, "inputs": inputs, "published": published, "markers": c.markers,
            "ten": ten, "full": c.people[["cycle", "female"]].reset_index(drop=True), "position": matched_position}


# -- UNCERTAINTY: one bootstrap draw of the comparison -----------------------
# Figure 4 legend: "bars in c and d are 95% intervals from 500
#   participant-bootstrap draws within cycle and sex. They exclude uncertainty
#   from marker selection and preprocessing." The full cohort is resampled
#   within cycle x sex and its matched members kept; the ten markers stay fixed.
def benchmark_rows(b: dict, draw: int) -> np.ndarray:
    """Matched-cohort rows of one paired bootstrap draw (R064): the full cohort resampled within
    cycle x sex (cycles in order, men first), then only its matched members kept."""
    rng = np.random.default_rng(np.random.SeedSequence([BENCHMARK_SEED, draw]))
    full = b["full"]
    blocks = [np.flatnonzero(full.cycle.eq(c).to_numpy() & full.female.eq(f).to_numpy())
              for c in sorted(full.cycle.unique()) for f in (0, 1)]
    rows = np.concatenate([rng.choice(block, len(block), replace=True) for block in blocks])
    rows = b["position"][rows]
    return rows[rows >= 0]


# -- CALCULATION: held-out concordance of each score -------------------------
# Results 4: "Held-out concordance was 0.694 for the full q score, 0.685 for the
#   ten-marker q score, 0.662 for the local refit and 0.649 for published
#   Phenotypic Age. The ten-marker q score exceeded the local refit by 0.024 ...
#   and published Phenotypic Age by 0.037". "every evaluation model included sex
#   and age".
# How, per held-out cycle: q38 and q10 directions from Cox(female + markers) on
#   training; local refit = Cox(female + Phenotypic Age inputs standardised on
#   training); evaluation Cox(female + age + score / training SD) on training;
#   concordance on the held-out cycle (prediction.concordance_counts). Pooled
#   over the six cycles.
#   Female HRs (Supp. Table S4.7; Results 4: "the female hazard ratio was 0.92
#   and 1.02 ... it remained 0.65 and 0.68") are descriptive: one model of
#   female + age + each participant's held-out score in the whole matched cohort.
def phenoage_metrics(b: dict, rows: np.ndarray | None = None) -> dict:
    """Held-out concordance of each score and their differences; descriptive female HRs."""
    people, z, inputs, published = b["people"], b["z"], b["inputs"], b["published"]
    if rows is not None:
        people, z = people.iloc[rows].reset_index(drop=True), z.iloc[rows].reset_index(drop=True)
        inputs, published = inputs[rows], published[rows]
    markers = b["markers"]
    scores = {k: np.full(len(people), np.nan) for k in ("q38", "q10", "clinical_refit")}
    held = {}
    for cycle, ten in b["ten"].items():
        train = people.cycle.ne(cycle).to_numpy()
        risk = P.training_risk_sets(people[train])
        fold = {}
        for name, cols in (("q38", markers), ("q10", ten)):
            X = z[cols].to_numpy(float)
            beta = P.fit(risk, np.column_stack([people.female, X])[train])
            fold[name] = X @ (beta[1:] / np.linalg.norm(beta[1:]))
        mu, sd = inputs[train].mean(axis=0), inputs[train].std(axis=0, ddof=1)
        Zc = (inputs - mu) / sd
        beta = P.fit(risk, np.column_stack([people.female, Zc])[train])
        fold["clinical_refit"] = Zc @ beta[1:]
        fold["phenoage"] = published
        held[cycle] = (train, fold)
    # Held-out evaluation: sex + age + score (score scaled by its training SD).
    names = ("q38", "q10", "clinical_refit", "phenoage")
    specs = {"baseline": (), **{n: (n,) for n in names}, "q10+phenoage": ("q10", "phenoage"),
             "q10+clinical_refit": ("q10", "clinical_refit")}
    counts = {k: np.zeros(2) for k in specs}
    for cycle, (train, fold) in held.items():
        test = ~train
        scaled = {n: fold[n] / fold[n][train].std(ddof=1) for n in names}
        risk = P.training_risk_sets(people[train])
        mv = people[test].reset_index(drop=True)
        for name, terms in specs.items():
            X = np.column_stack([people.female, people.RIDAGEYR] + [scaled[t] for t in terms])
            beta = P.fit(risk, X[train])
            counts[name] += P.concordance_counts(mv, X[test] @ beta)
        for n in scores:                                   # each participant's own held-out score
            scores[n][test] = fold[n][test]
    C = {k: v[0] / v[1] for k, v in counts.items()}
    out = {f"C_{k}": v for k, v in C.items()}
    for a, b_ in (("q10", "clinical_refit"), ("q10", "phenoage"), ("q38", "q10"), ("clinical_refit", "phenoage")):
        out[f"C_{a}_minus_{b_}"] = C[a] - C[b_]
    out["C_phenoage_added_to_q10"] = C["q10+phenoage"] - C["q10"]
    # Descriptive: female HR with each participant's held-out score in a pooled model (sex + age + score).
    risk = P.training_risk_sets(people)
    for n in names:
        s_ = (scores[n] if n != "phenoage" else published)
        s_ = s_ / s_.std(ddof=1)
        beta = P.fit(risk, np.column_stack([people.female, people.RIDAGEYR, s_]))
        out[f"female_hr_{n}"] = float(np.exp(beta[0]))
    return out


_B = {}


def _start_benchmark(ten: dict) -> None:
    """Worker set-up: the matched data, with the point run's ten markers per fold (no reselection)."""
    c = data.load_modern()
    keep = data.matched_rows(c)
    people = c.people[keep].reset_index(drop=True)
    inputs = data.phenoage_inputs(people.SEQN).assign(age=people.RIDAGEYR.to_numpy()).to_numpy(float)
    xb = -19.90667 + inputs @ PHENOAGE_WEIGHTS
    position = np.full(len(c.people), -1)
    position[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    _B.update(people=people, z=c.z[keep].reset_index(drop=True), inputs=inputs, markers=c.markers, ten=ten,
              published=141.50225 + (np.log(.0055305 * 1.51714 / .007692696) + xb) / .090165,
              full=c.people[["cycle", "female"]].reset_index(drop=True), position=position)


def benchmark_draw(draw: int) -> dict:
    return {"draw": draw, **phenoage_metrics(_B, benchmark_rows(_B, draw))}


BENCHMARK_INTERVALS = ("C_q38", "C_q10", "C_clinical_refit", "C_phenoage", "C_q10_minus_clinical_refit",
                       "C_q10_minus_phenoage", "C_q38_minus_q10", "C_clinical_refit_minus_phenoage",
                       "C_phenoage_added_to_q10")


def benchmark_intervals(boot: pd.DataFrame) -> dict:
    """95% percentile interval (numpy linear rule) of each benchmark concordance quantity."""
    return {k: models.percentile_interval(boot[k]) for k in BENCHMARK_INTERVALS}


# Figure 4's source tables, in the layout the accepted renderer reads
# (final_manuscript/figures_tables/section_4/figure4_*_source_data.csv). The drawing itself is done
# by that renderer; the pipeline supplies the numbers.
REDUCTION_DESIGNS = (("a", "LOCO", "loco"), ("a", "five_vs_five", "five_vs_five"),
                     ("b", "half_split", "iii_half_split"), ("b", "one_PSU_fold", "iii_tenfold"))
FOREST = (("c", "q38", "Full-panel q (38)"), ("c", "q10", "Ten-marker q"), ("c", "phenoage", "Published PhenoAge"),
          ("c", "clinical_refit", "Locally refitted PhenoAge inputs"),
          ("d", "q10_minus_phenoage", "Ten-marker q − published PhenoAge"),
          ("d", "q10_minus_clinical_refit", "Ten-marker q − locally refitted inputs"),
          ("d", "clinical_refit_minus_phenoage", "Locally refitted inputs − published PhenoAge"),
          ("d", "q38_minus_q10", "Full-panel q − ten-marker q"))


# Figure 4a,b source table: median retained gain and its 2.5th-97.5th
# percentiles across evaluations (legend: "bars span the 2.5th-97.5th
# percentiles across evaluations and are not confidence intervals").
def figure4_reduction(ev: pd.DataFrame) -> pd.DataFrame:
    """Figure 4a,b: median retained gain and its 2.5th-97.5th percentiles across evaluations
    (HD coordinates), by design and panel size; undefined ratios are left out and counted."""
    rows = []
    for panel, design, analysis in REDUCTION_DESIGNS:
        for size in SIZES:
            g = ev[ev.analysis.eq(analysis) & ev.coordinates.eq("HD") & ev.panel_size.eq(size)].retained_gain
            finite = g[np.isfinite(g.to_numpy(float))]
            rows.append({"panel": panel, "design": design, "panel_size": size, "estimate": finite.median(),
                         "lower": finite.quantile(.025), "upper": finite.quantile(.975), "evaluations": len(g),
                         "valid_ratios": len(finite), "invalid_ratios": len(g) - len(finite)})
    return pd.DataFrame(rows)


def figure4_forest(comparison: dict, bounds: dict) -> pd.DataFrame:
    """Figure 4c,d: held-out concordance of each score and the paired differences, with intervals."""
    return pd.DataFrame([{"panel": panel, "key": key, "label": label, "estimate": comparison[f"C_{key}"],
                          "lower": bounds.get(f"C_{key}", (np.nan, np.nan))[0],
                          "upper": bounds.get(f"C_{key}", (np.nan, np.nan))[1]} for panel, key, label in FOREST])


# Results 4 selection numbers: smallest adequate size, median retained gain
# ("0.913 in NHANES 1999-2018 and 0.934 in NHANES III"; one held out: "0.876 and
# 0.878"), distinct ten-marker sets ("181 distinct sets ... across 200
# evaluations"), and how often each locked marker was selected (Table 3).
def selection_numbers(ev: pd.DataFrame, locks: dict) -> tuple[list, pd.DataFrame, dict]:
    """The Results 4 selection summaries from all evaluations and the two locked panels.

    Returns (numbers, adequacy table, per-era locked-ten tables). Retained-gain medians skip
    evaluations whose ratio is undefined (see _ratio), as the reviewed summaries did.
    """
    numbers, tables, locked = [], [], {}
    for era, primary, secondary in (("modern", "five_vs_five", "loco"), ("iii", "iii_half_split", "iii_tenfold")):
        main_ev = ev[ev.analysis.eq(primary)]
        ad = adequacy(main_ev).assign(era=era)
        tables.append(ad)
        passing = [s_ for s_ in SIZES if ad[ad.panel_size.eq(s_)].passed.all()]
        numbers.append((f"{era}_smallest_adequate_size", min(passing) if passing else np.nan))
        for size in SIZES:
            hd = main_ev[main_ev.coordinates.eq("HD") & main_ev.panel_size.eq(size)]
            numbers.append((f"{era}_retained_gain_median_{size}", hd.retained_gain.median()))
        sec = ev[ev.analysis.eq(secondary) & ev.coordinates.eq("HD") & ev.panel_size.eq(10)]
        numbers.append((f"{era}_heldout_one_retained_gain_median_10", sec.retained_gain.median()))
        at10 = main_ev[main_ev.coordinates.eq("HD") & main_ev.panel_size.eq(10)].selected.str.split("|")
        numbers.append((f"{era}_distinct_sets_10", at10.map(lambda m: tuple(sorted(m))).nunique()))
        numbers.append((f"{era}_evaluations", len(at10)))
        ten = locks[era]
        locked[era] = pd.DataFrame({"rank": range(1, 11), "marker": ten,
                                    "times_selected_at_10": [sum(m in s_ for s_ in at10) for m in ten],
                                    "evaluations": len(at10)})
        if era == "modern":
            for m in ten:
                numbers.append((f"selected_at_10_percent.{m}", 100 * sum(m in s_ for s_ in at10) / len(at10)))
            numbers.append(("markers_selected_in_90pct",
                            sum(100 * sum(m in s_ for s_ in at10) / len(at10) >= 90 for m in ten)))
    return numbers, pd.concat(tables), locked


# ===========================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--evaluations", type=int, default=0, metavar="N",
                        help="run only the first N five-versus-five evaluations (a check)")
    parser.add_argument("--only", metavar="ID[,ID...]",
                        help="run only these evaluations, named by outer_id in inputs/panel_selection_tuning.csv "
                             "(a check; e.g. F09, the NHANES III held-out fold whose full gain is negative)")
    parser.add_argument("--check-draws", type=int, default=0, metavar="N",
                        help="also run benchmark bootstrap draws 0..N-1, compared with the archive")
    parser.add_argument("--bootstrap", action="store_true", help="run all 500 benchmark bootstrap draws")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    run = runs.start("s4")
    tuning = pd.read_csv(TUNING)
    jobs = tuning[tuning.analysis.isin(["five_vs_five", "loco", "iii_half_split", "iii_tenfold"])]
    all_evaluations = int(len(jobs))
    partial = bool(args.evaluations or args.only)
    if partial:                       # a check: the first N of each primary design, plus any named ones
        wanted = args.only.split(",") if args.only else []
        unknown = sorted(set(wanted) - set(jobs.outer_id))
        if unknown:
            parser.error(f"no evaluation named {', '.join(unknown)}")
        chosen = pd.concat([jobs[jobs.analysis.eq("five_vs_five")].head(args.evaluations),
                            jobs[jobs.analysis.eq("iii_half_split")].head(args.evaluations),
                            jobs[jobs.outer_id.isin(wanted)]])
        jobs = chosen[~chosen.index.duplicated()]
    with ProcessPoolExecutor(min(args.workers, len(jobs)), initializer=_start) as pool:
        ev = pd.DataFrame([r for out in pool.map(evaluate, jobs.to_dict("records")) for r in out])
    ev.to_csv(OUT / "evaluations.csv", index=False)
    # Evaluation-by-evaluation comparison with the archive is written with the benchmark draws below.
    # A full run compares all 5,664 rows (472 evaluations x 2 coordinate systems x 6 sizes).
    numbers = [("partial_run", 1)] if partial else []
    if not partial:
        _start()
        locks = {}
        for era, lock_name in (("modern", "lock"), ("iii", "iii_lock")):
            d, lock = _D[era], tuning[tuning.analysis.eq(lock_name)].iloc[0]
            locks[era] = elastic_net_ranking(frame_for(d["people"], d["coords"]["HD"], d["markers"]), d["markers"],
                                             lock.penalty, lock.l1_fraction)[:10]
        selection, adequacy_table, locked = selection_numbers(ev, locks)
        numbers += selection
        adequacy_table.to_csv(OUT / "adequacy.csv", index=False)
        figure4_reduction(ev).to_csv(OUT / "figure4_reduction_source_data.csv", index=False)
        for era, table in locked.items():
            table.to_csv(OUT / f"locked_ten_{era}.csv", index=False)
        t = locked["modern"]
        pd.DataFrame({"Rank": t["rank"], "Marker": t.marker,
                      "Selected in primary evaluations (%)": [f"{k} ({100 * k / n:.1f})" for k, n in
                                                              zip(t.times_selected_at_10, t.evaluations)]}) \
            .to_csv(OUT / "table3.csv", index=False, encoding="utf-8-sig")

    print("Phenotypic Age comparison")
    bench = benchmark_setup(tuning)
    comparison = phenoage_metrics(bench)

    def bench_draws(numbers_):
        with ProcessPoolExecutor(min(args.workers, len(numbers_)), initializer=_start_benchmark,
                                 initargs=(bench["ten"],)) as pool:
            return pd.DataFrame(list(pool.map(benchmark_draw, numbers_)))

    stored, fresh = OUT / "benchmark_bootstrap.csv", None
    if args.check_draws:
        fresh = bench_draws(list(range(args.check_draws)))
    if args.bootstrap:
        boot = bench_draws(list(range(BENCHMARK_DRAWS)))
        run.save_draws(boot, stored, ["draw"], BENCHMARK_SEED)
    else:
        boot = run.stored_draws(stored, ["draw"], BENCHMARK_SEED,
                                runs.draw_keys(None, range(BENCHMARK_DRAWS), era_column=False), fresh)
    bounds = benchmark_intervals(boot) if boot is not None else {}

    # One check file for both kinds of check (check_against_paper compares it with expected/s4_draws.csv):
    # evaluations keyed "<outer_id>|<coordinates>|<size>|<analysis>", benchmark draws keyed "benchmark|<draw>".
    check = ev.assign(draw=ev.outer_id + "|" + ev.coordinates + "|" + ev.panel_size.astype(str) + "|" + ev.analysis)[
        ["draw", "correlation", "retained_gain", "log_hr_men", "log_hr_women", "sex_contrast", "sex_contrast_retained"]]
    if fresh is not None:
        check = pd.concat([check, fresh[["draw", *BENCHMARK_INTERVALS]].assign(draw=lambda f: "benchmark|" + f.draw.astype(str))])
    run.write_check_draws(check, ["draw"])

    pd.DataFrame([{"key": k, "value": v, "low": bounds.get(k, (np.nan, np.nan))[0],
                   "high": bounds.get(k, (np.nan, np.nan))[1]} for k, v in comparison.items()]) \
        .to_csv(OUT / "phenoage.csv", index=False)
    figure4_forest(comparison, bounds).to_csv(OUT / "figure4_forest_source_data.csv", index=False)
    numbers = [(k, v, np.nan, np.nan) for k, v in numbers]
    numbers += [("matched_participants", len(bench["people"]), np.nan, np.nan),
                ("matched_deaths", int(bench["people"].event.sum()), np.nan, np.nan)]
    numbers += [(k, v, *bounds.get(k, (np.nan, np.nan))) for k, v in comparison.items()]
    pd.DataFrame(numbers, columns=["key", "value", "low", "high"]).assign(era="modern") \
        .to_csv(OUT / "numbers.csv", index=False)
    # full_draws: what an interval or a selection summary uses in a full run (checked against
    # expected/interval_verification.csv by check_against_paper).
    run.finish(evaluations=int(len(jobs)),
               full_draws={"interval": BENCHMARK_DRAWS, "summary": all_evaluations})
    print(f"done: outputs in {OUT}")

if __name__ == "__main__":
    main()
