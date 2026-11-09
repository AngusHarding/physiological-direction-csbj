"""Results 6: a within-sex screen selected distant, low-risk women and men.

The question
------------
Can pairing physiological distance with the mortality score select a cohort
that is (1) substantially distant from health, (2) at comparatively low
near-term risk of death, and (3) balanced between women and men?

How
---
Each era builds its own 32-marker scores, in NHANES 1999-2014 (9,362
participants) and NHANES III (4,440). Participants are split into five folds
of whole sampling units; for each fold, everything is learned from the other
four: the young-healthy reference, the mortality direction q, the caps of the
bounded distance, and the thresholds. Within each sex, the thresholds are the
60th percentile of distance and the 60th percentile of q among the training
participants of that sex. Groups:

  G1  distance below its cut-off, q at or below its cut-off   (lower distance, lower score)
  G2  distance at or above, q at or below                     (the candidates)
  G3  distance at or above, q above                           (higher distance, higher score)
  G4  distance below, q above                                 (lower distance, higher score)
  low score   q at or below its cut-off (G1 + G2), the screen without distance

Two distances: residual distance R (the part of the Euclidean length of the
coordinates not along q; the specified primary screen) and bounded distance
D_B (per marker, the raw distance from the same-sex reference median towards
the adverse side, capped at the training 99th percentile and scaled to 0-1;
the mean over markers).

Outcome: five-year all-cause mortality (Kaplan-Meier) within each group and
sex, its average over the sexes ("balanced"), and the risk expected for the
group's age mix (six age bands of the eligible participants of that sex).

Uncertainty: a Rao-Wu survey bootstrap of whole sampling units, 1,000 draws
per era, numbered from 0. Every learned stage (reference, direction, caps,
thresholds, and in NHANES 1999-2014 the 24-input rule) is refitted in every
draw; 95% intervals are the 2.5th and 97.5th percentiles of the draws (numpy's
linear rule), as in the reviewed runs (E7c/E7g stage 5, R115, R119, which all
used these same draws). The intervals therefore describe resampling of the
survey's sampling units with the whole procedure repeated; they do not cover
the choice of markers, transformations or the 60th-percentile level.

Scope
-----
Population   aged 52-79 with the 32 common markers: NHANES 1999-2014 (9,362; the
             cycles in which accidental deaths are identifiable) and NHANES III (4,440).
Endpoints    direction training: accident-censored mortality; group evaluation:
             five-year all-cause mortality (Kaplan-Meier, month-zero deaths count).
Estimands    five-year risk by group and sex, its average over the sexes, the risk
             expected from the group's age mix, selection fractions, and paired
             differences between groups, screens and the 24-input rule.
Split        every learned stage uses the other four folds only.
Assumptions  transformations fixed; point estimates do not use NHANES sample weights.
             Within-sex thresholds give each sex a similar selected share in
             training; they do not guarantee equal numbers of women and men in a trial.
Permitted    whether pairing distance with a low score selects distant, lower-risk
interpretation  participants of both sexes, compared with the score alone. Not a
             statement about trial costs or about benefit from an intervention.

What it makes (in outputs/s6/)
------------------------------
  groups.csv            every participant's out-of-fold scores and group under both screens
  metrics.csv           risks, age-expected risks, selection fractions and paired differences
  table4.csv / .xlsx    Table 4
  figure6.svg / .png    Figure 6 (interval bars once the bootstrap exists)
  screen_bootstrap.csv  the 2 x 1,000 survey-bootstrap draws (with --bootstrap) and its manifest
  intervals.csv         the 95% interval of every drawn quantity
  numbers.csv           every Results 6 number quoted in the text

    python -B s6_screen.py                  # point estimates (a few minutes); reuses stored draws
    python -B s6_screen.py --check-draws 2  # also draws 0-1 per era, compared with the archive
    python -B s6_screen.py --bootstrap      # all 2 x 1,000 draws (hours)
"""
from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import data                               # noqa: E402
import models                             # noqa: E402
import runs                               # noqa: E402
import settings as S                      # noqa: E402

OUT = S.output_folder("s6")
ERAS = ("modern", "iii")
RECIPE = {"modern": "modern_recipe", "iii": "nhanes3_recipe"}
SCREEN_CYCLES = S.CYCLES[:8]              # 1999-2014: accidental deaths identifiable
FOLD_SEED, BOOTSTRAP_SEED, DRAWS = 2026101402, 2026101401, 1000
LEVEL, CAP = 0.60, 0.99
AGE_CUTS = np.array([52, 55, 60, 65, 70, 75, 80])
GROUPS = ("E", "low_score", "high_distance", "G1", "G2", "G3", "G4")
# The exploratory 24-input rule: markers available from a comprehensive metabolic panel,
# a complete blood count with differential and platelet indices, and the examination.
NOT_ON_CMP_CBC = ("hdl_cholesterol", "hba1c", "ACR", "phosphorus", "uric_acid", "cholesterol", "triglycerides", "ldh")


# ===========================================================================
# Data
# ===========================================================================
def assign_folds(older: pd.DataFrame) -> dict:
    """Whole sampling units into five folds, balancing the number of men and of women.

    Units are visited largest first (after a seeded shuffle that breaks ties)
    and each goes to the fold where it least increases the squared relative
    imbalance. Uses only counts by sex, never outcomes.
    """
    counts = pd.crosstab(older.psu, older.sex).reindex(columns=[1, 2], fill_value=0)
    rng = np.random.default_rng(FOLD_SEED)
    counts = counts.loc[counts.index.to_numpy()[rng.permutation(len(counts))]]
    order = counts.sum(axis=1).sort_values(ascending=False, kind="stable").index
    target = counts.sum(axis=0).to_numpy(float) / 5
    loads, fold_of = np.zeros((5, 2)), {}
    for unit in order:
        n = counts.loc[unit].to_numpy(float)
        fold = int(np.argmin((((loads + n) / target) ** 2 - (loads / target) ** 2).sum(axis=1)))
        loads[fold] += n
        fold_of[unit] = fold
    return fold_of


def load_era(era: str) -> dict:
    frame = data.load_common32(era)
    markers = data.common32_markers()
    keys = markers.key.tolist()
    eligible = frame.age.between(52, 79) & frame.complete
    if era == "modern":
        frame["cycle"] = frame.stratum.str.split("|").str[0]
        mort = pd.read_csv(S.PREP / "01_imported/mortality.csv", usecols=["SEQN", "mortstat", "ucod_leading", "permth_exm"])
        eligible &= frame.cycle.isin(SCREEN_CYCLES) & frame.SEQN.isin(mort.SEQN[mort.permth_exm.notna()])
    else:
        mort = data.mortality_iii()[["SEQN", "mortstat", "ucod_leading", "permth_exm"]]
    older = frame[eligible].merge(mort, on="SEQN", how="left", validate="one_to_one").reset_index(drop=True)
    assert len(older) == {"modern": 9_362, "iii": 4_440}[era]
    older["fold"] = older.psu.map(assign_folds(older)).astype(int)
    months = older.permth_exm.to_numpy(float)
    position = pd.Series(np.arange(len(frame)), index=frame.SEQN).loc[older.SEQN].to_numpy()
    kinds = markers[RECIPE[era]].tolist()
    values = frame[keys].to_numpy(float)
    return {
        "era": era, "frame": frame, "older": older, "position": position,
        "T": np.column_stack([data.transform(values[:, j], k) for j, k in enumerate(kinds)]),
        "raw": values, "side": markers.adverse_side.to_numpy(),
        # Training the direction: principal endpoint, attained age; a death at month 0 counts at half a month.
        "entry": older.age.to_numpy(float),
        "exit": older.age.to_numpy(float) + np.where(months == 0, 1 / 24, months / 12),
        "event": data.principal_endpoint(older.mortstat, older.ucod_leading).to_numpy(),
        "strata": older.cycle.to_numpy(str) if era == "modern" else np.zeros(len(older), dtype=str),
        # Evaluating the groups: all-cause death within 60 months.
        "months": months, "died": older.mortstat.to_numpy(int), "sex": older.sex.to_numpy(int),
        "band": np.searchsorted(AGE_CUTS, older.age.to_numpy(), side="right") - 1}


# ===========================================================================
# One screen: learned fold by fold
# ===========================================================================
def weighted_mean_sd(X, w):
    n = w.sum()
    mu = (w[:, None] * X).sum(axis=0) / n
    return mu, np.sqrt((w[:, None] * (X - mu) ** 2).sum(axis=0) / (n - 1))


def crossfit(d: dict, frame_weights: np.ndarray, markers: list | None = None) -> dict:
    """Out-of-fold q, R and D_B for every eligible participant, and the two screens' groups.

    `markers` (column numbers) restricts the panel, as in the 24-input rule; every
    learned stage is refitted on that panel. Also returns the R screen under a
    pooled policy: one pair of thresholds for both sexes together.
    """
    if markers is not None:
        d = {**d, "T": d["T"][:, markers], "raw": d["raw"][:, markers], "side": d["side"][markers]}
    older, frame = d["older"], d["frame"]
    W = frame_weights[d["position"]]
    n = len(older)
    out = {k: np.full(n, np.nan) for k in ("q", "R", "DB", "q_cut", "R_cut", "DB_cut")}   # cut = own fold and sex
    high = {"R": np.zeros(n, bool), "DB": np.zeros(n, bool)}
    low = np.zeros(n, bool)
    pooled_high, pooled_low = np.zeros(n, bool), np.zeros(n, bool)
    for fold in range(5):
        test = older.fold.eq(fold).to_numpy()
        train = ~test
        held = set(older.psu[test])
        ref = (frame.is_reference & frame.complete & ~frame.psu.isin(held)).to_numpy() & (frame_weights > 0)
        sex_ref = frame.sex.to_numpy()[ref]
        Z = np.empty((n, d["T"].shape[1]))
        medians = np.empty((d["raw"].shape[1], 3))
        for s in (1, 2):
            rows = d["sex"] == s
            mu, sd = weighted_mean_sd(d["T"][ref][sex_ref == s], frame_weights[ref][sex_ref == s])
            Z[rows] = (d["T"][d["position"][rows]] - mu) / sd
            for j in range(d["raw"].shape[1]):          # bounded distance: same-sex reference medians
                medians[j, s] = models.weighted_quantile(d["raw"][ref][sex_ref == s, j],
                                                         frame_weights[ref][sex_ref == s], 0.5)
        # Mortality direction on training rows (female + coordinates), stratified by cycle.
        keep = train & (W > 0)
        risk = models.RiskSets(d["entry"][keep], d["exit"][keep], d["event"][keep], d["strata"][keep], W[keep])
        X = np.column_stack([older.sex.eq(2).to_numpy(float), Z])
        beta, ok = models.ridge_cox(risk, X[keep], np.zeros(X.shape[1]))
        if not ok:
            raise ArithmeticError("direction fit did not converge")
        u = beta[1:] / np.linalg.norm(beta[1:])
        q = Z @ u
        R = np.sqrt(np.maximum((Z * Z).sum(axis=1) - q * q, 0))
        # Bounded distance: adverse raw distance, capped at the training 99th percentile, scaled 0-1.
        raw = d["raw"][d["position"]]
        e = np.empty_like(raw)
        for s in (1, 2):
            rows = d["sex"] == s
            gap = raw[rows] - medians[:, s]
            e[rows] = np.where(d["side"] == "high_only", np.maximum(gap, 0),
                               np.where(d["side"] == "low_only", np.maximum(-gap, 0), np.abs(gap)))
        DB = np.empty(n)
        for s in (1, 2):
            rows, tr = d["sex"] == s, train & (d["sex"] == s)
            caps = np.array([models.weighted_quantile(e[tr, j], W[tr], CAP) for j in range(e.shape[1])])
            DB[rows] = np.minimum(e[rows] / caps, 1.0).mean(axis=1)
        for name, value in (("q", q), ("R", R), ("DB", DB)):
            out[name][test] = value[test]
        for s in (1, 2):                                   # within-sex training thresholds
            tr, ev = train & (d["sex"] == s), test & (d["sex"] == s)
            out["q_cut"][ev] = models.weighted_quantile(q[tr], W[tr], LEVEL)
            low[ev] = q[ev] <= out["q_cut"][ev]
            for name, value in (("R", R), ("DB", DB)):
                out[f"{name}_cut"][ev] = models.weighted_quantile(value[tr], W[tr], LEVEL)
                high[name][ev] = value[ev] >= out[f"{name}_cut"][ev]
        pooled_low[test] = q[test] <= models.weighted_quantile(q[train], W[train], LEVEL)     # both sexes together
        pooled_high[test] = R[test] >= models.weighted_quantile(R[train], W[train], LEVEL)
    groups = {}
    for rule in ("R", "DB"):
        h = high[rule]
        groups[rule] = {"E": np.ones(n, bool), "low_score": low, "high_distance": h,
                        "G1": ~h & low, "G2": h & low, "G3": h & ~low, "G4": ~h & ~low}
    return {"scores": out, "groups": groups, "pooled_candidates": pooled_high & pooled_low}


def evaluate(d: dict, groups: dict, frame_weights: np.ndarray) -> dict:
    """Five-year risk, age-expected risk and selection fraction by group and sex; sex averages."""
    W = frame_weights[d["position"]]
    t, died, sex, band = d["months"], d["died"], d["sex"], d["band"]
    expected_by_band = {s: np.array([models.kaplan_meier_risk(t[(sex == s) & (band == b)], died[(sex == s) & (band == b)],
                                                               W[(sex == s) & (band == b)]) for b in range(6)])
                        for s in (1, 2)}
    m = {}
    for rule, members in groups.items():
        for group, inside in members.items():
            for s in (1, 2):
                ix = inside & (sex == s)
                key = f"{rule}.{group}.sex{s}"
                m[f"{key}.risk"] = models.kaplan_meier_risk(t[ix], died[ix], W[ix])
                share = np.bincount(band[ix], weights=W[ix], minlength=6) / W[ix].sum()
                m[f"{key}.age_expected"] = float(share @ expected_by_band[s])
                m[f"{key}.fraction"] = W[ix].sum() / W[sex == s].sum()
            for q in ("risk", "age_expected"):
                m[f"{rule}.{group}.balanced.{q}"] = (m[f"{rule}.{group}.sex1.{q}"] + m[f"{rule}.{group}.sex2.{q}"]) / 2
            m[f"{rule}.{group}.balanced.age_residual"] = m[f"{rule}.{group}.balanced.risk"] - m[f"{rule}.{group}.balanced.age_expected"]
        for other in ("G3", "E", "low_score", "G1"):
            for s in ("sex1", "sex2", "balanced"):
                m[f"{rule}.G2_minus_{other}.{s}.risk"] = m[f"{rule}.G2.{s}.risk"] - m[f"{rule}.{other}.{s}.risk"]
    for s in ("sex1", "sex2", "balanced"):
        m[f"DB_minus_R.G2.{s}.risk"] = m[f"DB.G2.{s}.risk"] - m[f"R.G2.{s}.risk"]
    return m


def screen_metrics(d: dict, frame_weights: np.ndarray) -> tuple[dict, dict]:
    """Every drawn Results 6 quantity for one set of weights, and the 32-marker screen itself.

    The same function gives the point estimates (all weights one) and each bootstrap
    draw, so the two cannot drift apart. In NHANES 1999-2014 it also refits the
    exploratory 24-input rule with the same weights, so the 24-minus-32 differences
    are paired within each draw, as in R119.
    """
    screen = crossfit(d, frame_weights)
    m = evaluate(d, screen["groups"], frame_weights)
    if d["era"] == "modern":
        keys = data.common32_markers().key.tolist()
        subset = [j for j, k in enumerate(keys) if k not in NOT_ON_CMP_CBC]
        short = crossfit(d, frame_weights, subset)
        s24 = evaluate(d, short["groups"], frame_weights)
        m.update({f"cmpcbc24.{k}": v for k, v in s24.items() if k.startswith(("R.", "DB."))})
        for rule in ("R", "DB"):
            for side in ("sex1", "sex2", "balanced"):
                m[f"cmpcbc24_minus_32.{rule}.G2.{side}.risk"] = s24[f"{rule}.G2.{side}.risk"] - m[f"{rule}.G2.{side}.risk"]
        screen["cmpcbc24"] = short
    return m, screen


interval_table = models.interval_table          # 95% percentile intervals of the draws, per era


# ===========================================================================
# Table 4: clinical and lifestyle profiles of the bounded-distance groups (NHANES 1999-2014)
# ===========================================================================
CONDITIONS = {"hypertension": "Hypertension", "diabetes": "Diabetes", "coronary_heart_disease": "Coronary heart disease",
              "angina": "Angina", "heart_attack": "Heart attack", "congestive_heart_failure": "Congestive heart failure",
              "stroke": "Stroke", "kidney_fail": "Weak/failing kidneys", "liver_condition": "Liver condition",
              "arthritis": "Arthritis", "osteoporosis": "Osteoporosis", "gout": "Gout",
              "thyroid_disease": "Thyroid disease", "emphysema": "Emphysema", "chronic_bronchitis": "Chronic bronchitis",
              "copd": "COPD", "asthma": "Asthma", "cancer": "Cancer"}
LABEL_SUFFIX = {"osteoporosis": "‡", "gout": " [2007-14]", "copd": " [2013-14]"}
STANDARD_BANDS = ["52-54", "55-59", "60-64", "65-69", "70-74", "75-79"]


def recreational_activity(raw: pd.DataFrame) -> pd.Series:
    """Recreational physical activity (GPAQ, 2007-2014): None, Some or Higher (>= 1,200 MET-min/week).

    Vigorous (8 METs) and moderate (4 METs) recreation: a "no" answer counts
    as zero; a "yes" counts days x minutes x METs when days are 1-7 and minutes
    10-1439. A participant is classed only when the category is certain.
    """
    parts = []
    for gate, days, minutes, mets in (("PAQ650", "PAQ655", "PAD660", 8), ("PAQ665", "PAQ670", "PAD675", 4)):
        value = np.full(len(raw), np.nan)
        yes = (raw[gate].eq(1) & raw[days].between(1, 7) & raw[days].eq(raw[days].round())
               & raw[minutes].ge(10) & raw[minutes].lt(1440)).to_numpy()
        value[raw[gate].eq(2).to_numpy()] = 0.0
        value[yes] = (mets * raw.loc[yes, days] * raw.loc[yes, minutes]).to_numpy()
        parts.append(value)
    parts = np.column_stack(parts)
    complete, total = np.isfinite(parts).all(axis=1), np.nansum(parts, axis=1)
    out = np.full(len(raw), None, dtype=object)
    out[complete & (total == 0)] = "None"
    out[complete & (total > 0) & (total < 1200)] = "Some"
    out[total >= 1200] = "Higher"
    return pd.Series(out, index=raw.index)


def profile_data(d: dict, groups: dict, scores: dict) -> pd.DataFrame:
    """One row per eligible participant: group, age band and every profile characteristic."""
    older = d["older"]
    p = older[["SEQN", "sex", "age", "cycle"]].copy()
    p["bin"] = np.select([groups["G1"], groups["G2"]], ["G1", "G2"], "G3+G4")
    p["age_band"] = pd.cut(p.age, [51.5, 54.5, 59.5, 64.5, 69.5, 74.5, 79.5], labels=STANDARD_BANDS).astype(str)
    p["q"], p["DB"] = scores["q"], scores["DB"]
    p["died_5y"], p["months"] = d["died"], d["months"]
    flags = pd.read_csv(S.PREP / "02_harmonised/self_reported_disease_flags.csv", usecols=["SEQN"] + list(CONDITIONS))
    p = p.merge(flags, on="SEQN", how="left", validate="one_to_one")
    smoking = pd.read_csv(S.PREP / "07_master/analysis_master_raw.csv", usecols=["SEQN", "Smoking_status"])
    p = p.merge(smoking, on="SEQN", how="left", validate="one_to_one")
    p["current_smoking"] = p.Smoking_status.map({"Current": 1.0, "Former": 0.0, "Never": 0.0})
    # Diet quality: HEI-2020, 2005-2014, complete two-day recalls only; thirds of the eligible participants' scores.
    hei = pd.read_csv(S.PREP / "02_harmonised/hei_2020_total.csv")
    hei = hei.merge(pd.read_csv(S.PREP / "02_harmonised/hei_2020_component_scores.csv",
                                usecols=["SEQN", "WTDR2D", "all_components_scored"]), on="SEQN", validate="one_to_one")
    hei = hei.merge(pd.read_csv(S.PREP / "01_imported/hei_package_fped_inputs.csv",
                                usecols=["SEQN", "WTDR2D", "DR1TKCAL", "DR2TKCAL"]).rename(columns={"WTDR2D": "raw_weight"}),
                    on="SEQN", validate="one_to_one")
    valid = (hei.HEI_2020_total.between(0, 100) & hei.missing_component_count.eq(0) & hei.all_components_scored.eq(True)
             & hei.WTDR2D.gt(0) & hei.raw_weight.gt(0) & hei.DR1TKCAL.ge(0) & hei.DR2TKCAL.ge(0))
    p = p.merge(hei.loc[valid, ["SEQN", "HEI_2020_total"]], on="SEQN", how="left", validate="one_to_one")
    year = p.cycle.str[:4].astype(int)
    p["HEI"] = p.HEI_2020_total.where(year.between(2005, 2014))
    cut1, cut2 = np.quantile(p.HEI.dropna().to_numpy(), [1 / 3, 2 / 3])
    p["lower_HEI"] = np.where(p.HEI.notna(), (p.HEI < cut1).astype(float), np.nan)
    p["higher_HEI"] = np.where(p.HEI.notna(), (p.HEI >= cut2).astype(float), np.nan)
    activity = pd.read_csv(S.PREP / "01_imported/physical_activity_raw.csv",
                           usecols=["SEQN", "PAQ650", "PAQ655", "PAD660", "PAQ665", "PAQ670", "PAD675"])
    activity["PA"] = recreational_activity(activity)
    p = p.merge(activity[["SEQN", "PA"]], on="SEQN", how="left", validate="one_to_one")
    p.loc[~year.between(2007, 2014), "PA"] = None
    p["no_MVPA"] = np.where(p.PA.notna(), p.PA.eq("None").astype(float), np.nan)
    p["higher_MVPA"] = np.where(p.PA.notna(), p.PA.eq("Higher").astype(float), np.nan)
    return p


def table4(p: pd.DataFrame) -> pd.DataFrame:
    """Age-standardised percentages within sex: each age band's share of all eligible participants
    of that sex is the standard; percentages use those who answered.

    Markers: § fewer than 10 cases, or fewer than 20 responses in some age band;
    † the value differs from the same-sex G1 value by more than 2 SD, where the SD
    is the G1 within-sex SD pooled across men and women (the colour scale of the
    published table is clipped there).
    """
    standard = {s: p[p.sex.eq(s)].age_band.value_counts(normalize=True).reindex(STANDARD_BANDS).to_numpy()
                for s in (1, 2)}

    def standardised(cell, column, s):
        rates, smallest = [], np.inf
        for band in STANDARD_BANDS:
            x = cell.loc[cell.age_band.eq(band), column].dropna()
            rates.append(x.mean() if len(x) else np.nan)
            smallest = min(smallest, len(x))
        sparse = smallest < 20 or cell[column].eq(1).sum() < 10
        return 100 * float(np.dot(standard[s], rates)), sparse

    def g1_sd(column, scale=1.0):
        parts = [p.loc[p.sex.eq(s) & p.bin.eq("G1"), column].dropna() for s in (1, 2)]
        pooled = sum((len(x) - 1) * x.var(ddof=1) for x in parts) / sum(len(x) - 1 for x in parts)
        return scale * np.sqrt(pooled)

    rows = [("Mean q", "q", "mean", 3), ("Mean bounded index (0–1)", "DB", "mean", 3)]
    rows += [(f"{label}, %{LABEL_SUFFIX.get(key, '')}", key, "percent", 1) for key, label in CONDITIONS.items()]
    rows += [("Lower HEI third, % [2005-14]", "lower_HEI", "percent", 1),
             ("Higher HEI third, % [2005-14]", "higher_HEI", "percent", 1),
             ("No recreational MVPA, % [2007-14]", "no_MVPA", "percent", 1),
             ("Higher recreational MVPA, % [2007-14]", "higher_MVPA", "percent", 1),
             ("Current smoking, %", "current_smoking", "percent", 1),
             ("Five-year mortality, %", None, "mortality", 2)]
    columns = {}
    for s, sex in ((1, "Men"), (2, "Women")):
        for b in ("G1", "G2", "G3+G4"):
            columns[f"{sex} {b} (n={int((p.sex.eq(s) & p.bin.eq(b)).sum())})"] = []
    for label, column, kind, digits in rows:
        for s, sex in ((1, "Men"), (2, "Women")):
            reference = None
            for b in ("G1", "G2", "G3+G4"):
                cell = p[p.sex.eq(s) & p.bin.eq(b)]
                sparse = False
                if kind == "mean":
                    value = cell[column].mean()
                elif kind == "percent":
                    value, sparse = standardised(cell, column, s)
                else:
                    value = 100 * models.kaplan_meier_risk(cell.months.to_numpy(), cell.died_5y.to_numpy())
                text = f"{value:.{digits}f}"
                if b == "G1":
                    reference = value
                elif kind != "mortality":
                    sd = g1_sd(column, 1.0 if kind == "mean" else 100.0)
                    if abs(value - reference) / sd > 2:
                        text += "†"
                if sparse:
                    text += "§"
                columns[f"{sex} {b} (n={len(cell)})"].append(text)
    return pd.DataFrame({"Measure": [r[0] for r in rows], **columns})


def figure6(groups: pd.DataFrame, metrics: pd.DataFrame, intervals: pd.DataFrame | None = None) -> None:
    """Figure 6: a, b every participant relative to their own training cut-offs (bounded-distance
    screen); c five-year mortality by group and sex; d three selection rules, sexes weighted equally."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import displays as V

    colour = {"G1": V.MEN, "G2": V.GREEN, "G3": "#B5452F", "G4": "#A866B0"}
    fig = plt.figure(figsize=(185 * V.MM, 200 * V.MM), facecolor="white")
    for panel, era, x0 in (("a", "modern", .10), ("b", "iii", .61)):
        ax = fig.add_axes([x0, .63, .34, .27])
        g = groups[groups.era.eq(era)]
        for name, c in colour.items():
            cell = g[g.group_DB.eq(name)]
            ax.scatter(cell.q - cell.q_cut, cell.DB - cell.DB_cut, s=1.5, color=c, alpha=.5, lw=0, rasterized=True)
        ax.axhline(0, color=V.MUTED, lw=.8, ls=(0, (3, 2)))
        ax.axvline(0, color=V.MUTED, lw=.8, ls=(0, (3, 2)))
        ax.set(xlim=(-4.5, 8.5), ylim=(-.17, .32), xlabel="q minus training 60th-percentile cut-off",
               ylabel="Bounded distance minus\ntraining 60th-percentile cut-off")
        for name, (x, y, ha, va) in {"G2": (.02, .97, "left", "top"), "G3": (.98, .97, "right", "top"),
                                     "G1": (.02, .03, "left", "bottom"), "G4": (.98, .03, "right", "bottom")}.items():
            ax.text(x, y, name, transform=ax.transAxes, ha=ha, va=va, color=colour[name], fontweight="bold")
        m = metrics[metrics.era.eq(era)].set_index("quantity").estimate
        ax.set_title(f"{'NHANES 1999–2014' if era == 'modern' else 'NHANES III'}", loc="left", pad=16, fontweight="bold")
        ax.text(0, 1.02, f"{int(m['participants']):,} participants · {int(m['deaths_within_5y']):,} five-year deaths",
                transform=ax.transAxes, fontsize=6.7, color=V.MUTED, va="bottom")
        fig.text(x0 - .075, .93, panel, fontweight="bold", fontsize=10)
    fig.legend(handles=[Line2D([0], [0], marker="o", lw=0, color=c, label=label) for c, label in
                        ((colour["G1"], "G1: lower distance / lower q"), (colour["G2"], "G2: candidates"),
                         (colour["G3"], "G3: higher distance / higher q"), (colour["G4"], "G4: lower distance / higher q"))],
               loc="center", bbox_to_anchor=(.55, .555), ncol=2, frameon=False, fontsize=7)

    def value(era, key):
        m = metrics[metrics.era.eq(era)].set_index("quantity").estimate
        low = high = np.nan
        if intervals is not None and ((intervals.era == era) & (intervals.key == key)).any():
            r = intervals[(intervals.era == era) & (intervals.key == key)].iloc[0]
            low, high = r.low, r.high
        return 100 * m[key], 100 * low, 100 * high

    c = fig.add_axes([.17, .10, .30, .36])
    ticks = []
    for block, era in enumerate(("modern", "iii")):
        for i, group in enumerate(("G1", "G2", "G3", "G4")):
            y = 8.5 - block * 4.6 - i
            ticks.append((y, group))
            for s, col, offset in ((1, V.MEN, .12), (2, V.WOMEN, -.12)):
                v, lo, hi = value(era, f"DB.{group}.sex{s}.risk")
                c.hlines(y + offset, lo, hi, color=col, lw=1.4)
                c.plot(v, y + offset, "o", color=col, ms=5)
    c.text(34, 8.95, "NHANES 1999–2014", ha="right", fontsize=7, fontweight="bold")
    c.text(34, 4.35, "NHANES III", ha="right", fontsize=7, fontweight="bold")
    c.axhline(4.75, color=V.GRID, lw=1)
    c.set_yticks([y for y, _ in ticks])
    c.set_yticklabels([t for _, t in ticks])
    c.set(xlim=(0, 35), ylim=(0.4, 9.3), xlabel="Five-year all-cause mortality, %")
    c.legend(handles=[Line2D([0], [0], marker="o", lw=0, color=V.MEN, label="Men"),
                      Line2D([0], [0], marker="o", lw=0, color=V.WOMEN, label="Women")],
             loc="lower center", bbox_to_anchor=(.5, 1.02), ncol=2, frameon=False)
    d = fig.add_axes([.70, .10, .27, .36])
    ticks = []
    for block, era in enumerate(("modern", "iii")):
        for i, (key, label, col) in enumerate((("DB.low_score.balanced.risk", "Low q only", V.MUTED),
                                               ("R.G2.balanced.risk", "R + low q\n(specified primary)", V.PURPLE),
                                               ("DB.G2.balanced.risk", "Bounded distance\n+ low q", V.GREEN))):
            y = 6.4 - block * 3.6 - i
            ticks.append((y, label))
            v, lo, hi = value(era, key)
            d.hlines(y, lo, hi, color=col, lw=1.4)
            d.plot(v, y, "o", color=col, ms=6)
    d.text(11.8, 6.85, "NHANES 1999–2014", ha="right", fontsize=7, fontweight="bold")
    d.text(11.8, 3.25, "NHANES III", ha="right", fontsize=7, fontweight="bold")
    d.axhline(3.6, color=V.GRID, lw=1)
    d.set_yticks([y for y, _ in ticks])
    d.set_yticklabels([t for _, t in ticks])
    d.set(xlim=(0, 12), ylim=(.4, 7.2), xlabel="Five-year all-cause mortality, %")
    d.set_title("Sexes weighted equally", fontsize=7, color=V.MUTED)
    for ax in (c, d):
        ax.grid(axis="x", color=V.GRID, lw=.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0, pad=6)
    for key, x, y, title in (("c", .025, .52, "Mortality across the four groups"),
                             ("d", .52, .52, "The distance criterion and risk")):
        fig.text(x, y, key, fontweight="bold", fontsize=10, va="center")
        fig.text(x + .033, y, title, fontsize=8.6, fontweight="bold", va="center")
    note = "" if intervals is not None else " Interval bars appear once the survey bootstrap has been run."
    fig.text(.5, .02, "a–b: each participant relative to their own sex- and fold-specific training cut-offs." + note,
             ha="center", fontsize=6.5, color=V.MUTED)
    V.save_figure(fig, OUT, "figure6")


def replicate_weights(frame: pd.DataFrame, era_index: int, draw: int) -> np.ndarray:
    """One survey-bootstrap draw (Rao-Wu): n - 1 of each stratum's n units drawn with replacement.

    Each draw has its own random stream (seed, era, draw number), so any draw
    can be regenerated alone. Draws are numbered from 0, as in the archive.
    """
    rng = np.random.default_rng(np.random.SeedSequence([BOOTSTRAP_SEED, era_index, draw]))
    weight = {}
    for _, units in frame.groupby("stratum", sort=True).psu.unique().items():
        units = sorted(units)
        n = len(units)
        counts = np.bincount(rng.integers(0, n, size=n - 1), minlength=n)
        weight.update({u: c * n / (n - 1) for u, c in zip(units, counts)})
    return frame.psu.map(weight).to_numpy(float)


_LOADED = {}


def one_draw(job: tuple) -> dict:
    """Every learned stage refitted with one draw's weights, then every quantity evaluated."""
    era, draw = job
    if era not in _LOADED:                       # each worker process loads an era once
        _LOADED[era] = load_era(era)
    d = _LOADED[era]
    m, _ = screen_metrics(d, replicate_weights(d["frame"], ERAS.index(era), draw))
    return {"era": era, "draw": draw, **m}


def race(d: dict) -> pd.Series:
    """Race and ethnicity of the eligible participants, in the four groups both eras share.

    NHANES III folds other Hispanic participants into "Other", so in NHANES
    1999-2018 other Hispanic and other or multiracial participants are joined too.
    """
    if d["era"] == "modern":
        r = pd.read_csv(S.PREP / "01_imported/socioeconomic.csv", usecols=["SEQN", "RIDRETH1"]).set_index("SEQN").RIDRETH1
        labels = {1: "Mexican American", 2: "Other", 3: "Non-Hispanic White", 4: "Non-Hispanic Black", 5: "Other"}
    else:
        r = data.read_nhanes3("adult", ["SEQN", "DMARETHN"]).set_index("SEQN").DMARETHN
        labels = {1: "Non-Hispanic White", 2: "Non-Hispanic Black", 3: "Mexican American", 4: "Other"}
    return r.loc[d["older"].SEQN].astype(int).map(labels).reset_index(drop=True)


def composition(d: dict, screen: dict) -> dict:
    """Who is selected: the share of women under each policy; race shares; selection by age and race."""
    sex, m = d["sex"], {}
    for rule in ("R", "DB"):
        g2 = screen["groups"][rule]["G2"]
        m[f"{rule}.within_sex.women_share_of_candidates"] = (g2 & (sex == 2)).sum() / g2.sum()
    pooled = screen["pooled_candidates"]
    m["R.pooled.women_share_of_candidates"] = (pooled & (sex == 2)).sum() / pooled.sum()
    groups = race(d)
    age_band = pd.cut(d["older"].age, [51.5, 59.5, 69.5, 79.5], labels=["52-59", "60-69", "70-79"]).astype(str)
    for rule in ("R", "DB"):
        g2 = screen["groups"][rule]["G2"]
        for label in groups.unique():
            inside = (groups == label).to_numpy()
            m[f"{rule}.race_share_candidates.{label}"] = (g2 & inside).sum() / g2.sum()
            m[f"{rule}.race_share_eligible.{label}"] = inside.mean()
            m[f"{rule}.selected_fraction_by_race.{label}"] = g2[inside].mean()
        for label in ("52-59", "60-69", "70-79"):
            m[f"{rule}.selected_fraction_by_age.{label}"] = g2[(age_band == label).to_numpy()].mean()
        m[f"{rule}.selected_fraction"] = g2.mean()
        m[f"{rule}.candidates"] = int(g2.sum())
    return m


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check-draws", type=int, default=0, metavar="N",
                        help="run only draws 0..N-1 per era and compare them with the archive (minutes)")
    parser.add_argument("--bootstrap", action="store_true", help="run all 2 x 1,000 survey-bootstrap draws (hours)")
    parser.add_argument("--workers", type=int, default=8, help="processes for the draws")
    args = parser.parse_args()
    run = runs.start("s6")

    rows, metrics = [], []
    for era in ERAS:
        d = load_era(era)
        ones = np.ones(len(d["frame"]))
        m, screen = screen_metrics(d, ones)
        m.update(composition(d, screen))
        if era == "modern":                       # counts for the exploratory 24-input rule
            short = screen["cmpcbc24"]
            for rule in ("R", "DB"):
                g2 = short["groups"][rule]["G2"]
                m[f"cmpcbc24.{rule}.candidates"] = int(g2.sum())
                m[f"cmpcbc24.{rule}.selected_fraction"] = g2.mean()
                m[f"cmpcbc24.{rule}.women_share_of_candidates"] = (g2 & (d["sex"] == 2)).sum() / g2.sum()
                m[f"cmpcbc24.{rule}.candidates_also_selected_by_32"] = int((g2 & screen["groups"][rule]["G2"]).sum())
        if era == "modern":
            import displays as V
            V.write_table(table4(profile_data(d, screen["groups"]["DB"], screen["scores"])), OUT, "table4")
        metrics.append(pd.DataFrame({"era": era, "quantity": list(m), "estimate": list(m.values())}))
        g = d["older"][["SEQN", "sex", "age", "fold"]].assign(era=era, **screen["scores"])
        for rule, members in screen["groups"].items():
            for group in ("G1", "G2", "G3", "G4"):
                g.loc[members[group], f"group_{rule}"] = group
        rows.append(g)
        five_year = int(((d["months"] <= 60) & (d["died"] == 1)).sum())
        metrics.append(pd.DataFrame({"era": era, "quantity": ["participants", "deaths_within_5y"],
                                     "estimate": [len(d["older"]), five_year]}))
        print(f"  {era}: done", flush=True)
    pd.concat(rows).to_csv(OUT / "groups.csv", index=False)

    def survey_draws(numbers):
        jobs = [(era, b) for era in ERAS for b in numbers]
        with ProcessPoolExecutor(min(len(jobs), args.workers)) as pool:
            return pd.DataFrame(list(pool.map(one_draw, jobs, chunksize=4)))

    stored, seed = OUT / "screen_bootstrap.csv", {"bootstrap": BOOTSTRAP_SEED, "folds": FOLD_SEED}
    fresh = None
    if args.check_draws:
        fresh = survey_draws(range(args.check_draws))
        run.write_check_draws(fresh, ["era", "draw"])
        print(f"{args.check_draws} draws per era written to check_draws.csv")
    if args.bootstrap:
        boot = survey_draws(range(DRAWS))
        run.save_draws(boot, stored, ["era", "draw"], seed)
    else:
        boot = run.stored_draws(stored, ["era", "draw"], seed, runs.draw_keys(ERAS, range(DRAWS)), fresh)
    intervals = interval_table(boot) if boot is not None else None
    if intervals is not None:
        intervals.to_csv(OUT / "intervals.csv", index=False)

    metrics = pd.concat(metrics)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    figure6(pd.concat(rows), metrics, intervals)
    # numbers.csv: risks, fractions and shares in percent, as the paper prints them.
    numbers = metrics.rename(columns={"quantity": "key", "estimate": "value"})
    if intervals is not None:
        numbers = numbers.merge(intervals[["era", "key", "low", "high"]], on=["era", "key"], how="left")
    else:
        numbers = numbers.assign(low=np.nan, high=np.nan)
    percent = numbers.key.str.contains("risk|age_expected|age_residual|fraction|share")
    numbers.loc[percent, ["value", "low", "high"]] *= 100
    numbers.to_csv(OUT / "numbers.csv", index=False)
    run.finish(full_draws={"interval": DRAWS})       # per era
    print(f"done: outputs in {OUT}")


if __name__ == "__main__":
    main()
