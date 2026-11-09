"""Results 1: physiological distance did not track women's lower mortality.

The question
------------
If women outlive men because their physiology has drifted less far from
health, then (1) women should sit closer than men to a young-healthy
reference, and (2) allowing for that distance should bring women's mortality
towards men's. This script tests both in the two analysis cohorts.

Scope
-----
Population   the two analysis cohorts (data.py: aged 52-79, complete marker panel,
             linked follow-up); Part B uses broader adult samples (ages 18-79).
Endpoint     accident-censored mortality (data.py).
Estimands    (a) women minus men in each distance, in cohort SDs, adjusted for an
             age spline and survey period (ordinary least squares); (b) the female
             hazard ratio with sex alone and with each distance added; (c) the hazard
             ratio per SD of each distance. Cox models on attained age, each era with
             its reviewed program (models.py).
Split        none: every estimate describes the whole cohort; nothing learned from
             outcomes is evaluated on new participants here.
Uncertainty  500 participant-bootstrap draws (within cycle x sex, or sex x phase in
             NHANES III); distances are not re-derived in a draw, so the intervals
             hold the reference and the transformations fixed.
Assumptions  proportional hazards on attained age; the prepared transformations and
             young-healthy references are taken as given.
Permitted    associations: whether women are closer to the reference on a given
interpretation  metric, and whether allowing for a distance moves the female hazard
             ratio. Not causal, and not a statement about any metric not computed.

What it makes (in outputs/s1/)
------------------------------
  table1.csv / .xlsx       Table 1: characteristics of women and men
  figure1.svg / .png       Figure 1, panels a-d
  ten_year_curves.csv      Figure 1a,b: ten-year probability of death by age and sex
  distance_estimates.csv   Figure 1c,d: sex differences in distance and the
                           female-to-male hazard ratios, with 95% intervals
  distance_bootstrap.csv   the 500 bootstrap draws behind those intervals
                           (made only with --bootstrap; minutes, not hours)
  numbers.csv              every Results 1 number quoted in the text

Run from this folder:
    python -B s1_sex_gap.py                  # estimates, table and figure (under a minute); reuses stored draws
    python -B s1_sex_gap.py --check-draws 2  # also draws 1-2, compared with the stored and archived draws
    python -B s1_sex_gap.py --bootstrap      # redo all 500 draws per era (minutes)

The four parts
--------------
  A  Distance and mortality (Fig. 1c,d): for each cohort, measure each
     participant's adverse distance (DA) and Mahalanobis distance (DM); compare
     women with men after allowing for age and survey period; then fit Cox
     models with sex alone and with sex plus each distance.
  B  Ten-year mortality by age (Fig. 1a,b), in broader adult samples.
  C  Table 1.
  D  Figure 1.
"""
from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")       # one thread per process; parallelism comes from workers

import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
import statsmodels.api as sm              # noqa: E402
import statsmodels.formula.api as smf     # noqa: E402

import data                               # noqa: E402
import displays as V                      # noqa: E402
import measures as M                      # noqa: E402
import models                             # noqa: E402
import runs                               # noqa: E402
import settings as S                      # noqa: E402

OUT = S.output_folder("s1")
DISTANCES = {"DA": "Adverse distance", "DM": "Mahalanobis distance", "DE": "Euclidean distance"}


# ===========================================================================
# Part A. Distance and mortality
# ===========================================================================
def distances(c: data.Cohort) -> pd.DataFrame:
    """Each participant's three distances, in their natural units."""
    return pd.DataFrame({"DA": M.adverse_distance(c.x),
                         "DM": M.mahalanobis_distance(c.z, c.people.female, c.ref_z),
                         "DE": M.euclidean_distance(c.z)})


def estimates(era: str, people: pd.DataFrame, dist: pd.DataFrame) -> dict:
    """All Part A quantities for one sample (the cohort itself or one bootstrap draw).

    For each distance:
      contrast    women minus men, in cohort SDs, adjusted for age and period
      log_hr      female log hazard ratio with that distance in the model
      per_sd      log hazard ratio per SD of the distance
    plus sex_only, the female log hazard ratio with no distance.
    Distances are standardised within the sample, so "per SD" means per SD of
    that sample.
    """
    people, dist = people.reset_index(drop=True), dist.reset_index(drop=True)
    frame = people.copy()
    for name in DISTANCES:
        frame[name] = M.standardise(dist[name].to_numpy(float))
    out = {"sex_only": models.cox(era, frame, ["female"])["female"]}
    for name in DISTANCES:
        fit = models.cox(era, frame, ["female", name])
        out[f"{name}|contrast"] = models.sex_contrast(era, people, dist[name].to_numpy(float))
        out[f"{name}|log_hr"] = fit["female"]
        out[f"{name}|per_sd"] = fit[name]
    return out


# The bootstrap runs in several worker processes; each loads the cohort once.
_COHORT = {}


def _start_worker(era: str) -> None:
    c = data.load(era)
    _COHORT.update(era=era, people=c.people, dist=distances(c),
                   strata=models.bootstrap_strata(era, c.people))


def _one_draw(draw: int) -> dict:
    era = _COHORT["era"]
    rows = models.bootstrap_rows(era, _COHORT["strata"], draw)
    return {"era": era, "draw": draw,
            **estimates(era, _COHORT["people"].iloc[rows], _COHORT["dist"].iloc[rows])}


def run_bootstrap(workers: int, draw_numbers) -> pd.DataFrame:
    """Participant-bootstrap draws for both cohorts (draw numbers 1..500 for the full set).

    Each draw resamples participants with replacement within groups (see
    models.bootstrap_strata) and recomputes every Part A quantity. Distances
    are not re-derived in a draw: the reference is held fixed, so the
    intervals describe sampling of participants given that reference.
    """
    results, draw_numbers = [], list(draw_numbers)
    for era in S.ERAS:
        with ProcessPoolExecutor(min(workers, len(draw_numbers)), initializer=_start_worker,
                                 initargs=(era,)) as pool:
            results += list(pool.map(_one_draw, draw_numbers, chunksize=4))
        print(f"  {era}: {len(draw_numbers)} draws done", flush=True)
    return pd.DataFrame(results)


def summarise(point: dict, draws: pd.DataFrame | None) -> pd.DataFrame:
    """One row per quantity: estimate, 95% interval, and hazard ratios where relevant."""
    rows = []
    for era, values in point.items():
        for key, value in values.items():
            distance, quantity = key.split("|") if "|" in key else ("none", key)
            row = {"era": era, "distance": distance, "quantity": quantity, "estimate": value,
                   "low": np.nan, "high": np.nan}
            if draws is not None:
                row["low"], row["high"] = models.percentile_interval(draws.loc[draws.era.eq(era), key])
            if quantity in ("sex_only", "log_hr", "per_sd"):     # report as hazard ratios too
                row.update(hr=np.exp(value), hr_low=np.exp(row["low"]), hr_high=np.exp(row["high"]))
            rows.append(row)
    return pd.DataFrame(rows)


# ===========================================================================
# Part B. Ten-year mortality by age, in broader adult samples
# ===========================================================================
def broad_modern() -> pd.DataFrame:
    """NHANES 1999-2008 adults aged 18-79 with complete Table 1 covariates and follow-up.

    1999-2008 only: these are the cycles in which every survivor has ten years
    of follow-up. No biomarker panel is required.
    """
    p = S.PREP
    d = pd.read_csv(p / "04_cohorts/eligibility_flags.csv", usecols=["SEQN", "cycle", "RIDAGEYR",
                                                                     "RIAGENDR", "eligible_adult"])
    d = d[d.eligible_adult.astype(bool)]
    for file, cols in (("01_imported/socioeconomic.csv", ["RIDRETH1", "DMDEDUC2", "INDFMPIR"]),
                       ("01_imported/mortality.csv", ["mortstat", "ucod_leading", "permth_exm"]),
                       ("02_harmonised/smoking_status.csv", ["Smoking_status"]),
                       ("02_harmonised/alcohol_drinks_per_week.csv", ["drinks_per_week"]),
                       ("02_harmonised/physical_activity_inclusion_preserving.csv", ["physically_active"])):
        d = d.merge(pd.read_csv(p / file, usecols=["SEQN"] + cols), on="SEQN", how="left",
                    validate="one_to_one")
    covariates = ["RIDRETH1", "DMDEDUC2", "INDFMPIR", "mortstat", "permth_exm", "Smoking_status",
                  "drinks_per_week", "physically_active"]
    d = d[d[covariates].notna().all(axis=1) & d.permth_exm.gt(0)]
    d = d[d.cycle.isin(S.CYCLES[:5])]
    return pd.DataFrame({"female": d.RIAGENDR.eq(2).astype(int), "age": d.RIDAGEYR.astype(float),
                         "months": d.permth_exm, "dead": d.mortstat,
                         "event": data.principal_endpoint(d.mortstat, d.ucod_leading)})


def broad_iii(adult: pd.DataFrame) -> pd.DataFrame:
    """NHANES III adults aged 18-79 with linked follow-up."""
    d = adult[adult.HSAGEIR.between(18, 79)].merge(data.mortality_iii(), on="SEQN", how="inner")
    d = d[d.permth_exm.gt(0)]
    return pd.DataFrame({"female": d.HSSEX.eq(2).astype(int), "age": d.HSAGEIR.astype(float),
                         "months": d.permth_exm, "dead": d.mortstat, "event": d.event})


def ten_year_curves(era: str, d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Probability of an endpoint death within ten years, by age, separately for each sex.

    Survivors followed for less than ten years are left out (their ten-year
    status is unknown). A logistic regression on a smooth age curve (natural
    cubic spline, 4 degrees of freedom) gives the probability and its 95%
    model-based band at each year of age.
    """
    d = d[~(d.dead.eq(0) & d.months.lt(120))].copy()
    d["died_10y"] = (d.event.eq(1) & d.months.le(120)).astype(int)
    ages = pd.DataFrame({"age": np.arange(d.age.min(), d.age.max() + 1)})
    curves, counts = [], {"era": era, "n": len(d), "deaths_10y": int(d.died_10y.sum())}
    for female, sex in ((0, "Men"), (1, "Women")):
        s = d[d.female.eq(female)]
        fit = smf.glm("died_10y ~ cr(age, df=4, constraints='center')", data=s,
                      family=sm.families.Binomial()).fit()
        p = fit.get_prediction(ages).summary_frame()
        curves.append(pd.DataFrame({"era": era, "sex": sex, "age": ages.age, "probability": p["mean"],
                                    "low": p.mean_ci_lower, "high": p.mean_ci_upper}))
    return pd.concat(curves), counts


# ===========================================================================
# Part C. Table 1
# ===========================================================================
def characteristics_modern(people: pd.DataFrame) -> pd.DataFrame:
    p = S.PREP
    d = people[["SEQN", "female", "RIDAGEYR", "followup_years", "event"]].rename(columns={"RIDAGEYR": "age"})
    for file, cols in (("01_imported/socioeconomic.csv", ["RIDRETH1", "DMDEDUC2", "INDFMPIR"]),
                       ("02_harmonised/smoking_status.csv", ["Smoking_status"]),
                       ("02_harmonised/alcohol_drinks_per_week.csv", ["drinks_per_week"]),
                       ("02_harmonised/physical_activity.csv", ["physically_active"]),
                       ("02_harmonised/anthropometry.csv", ["BMXWT", "BMXHT"]),
                       ("02_harmonised/self_reported_disease_flags.csv",
                        ["historical_13_complete", "historical_13_count"])):
        d = d.merge(pd.read_csv(p / file, usecols=["SEQN"] + cols), on="SEQN", how="left",
                    validate="one_to_one")
    complete = d.historical_13_complete.fillna(False).astype(bool)
    return pd.DataFrame({
        "female": d.female, "age": d.age,
        "white": d.RIDRETH1.eq(3), "black": d.RIDRETH1.eq(4), "mexican": d.RIDRETH1.eq(1),
        "college": d.DMDEDUC2.eq(5), "poor": d.INDFMPIR.lt(1.0),
        "never_smoked": d.Smoking_status.eq("Never"), "current_smoker": d.Smoking_status.eq("Current"),
        "drinks": d.drinks_per_week,
        "inactive": d.physically_active.isin([0, False]),
        "bmi": d.BMXWT / (d.BMXHT / 100) ** 2,
        "conditions": d.historical_13_count.where(complete),        # of 13 conditions
        "followup": d.followup_years, "death": d.event})


# NHANES III questionnaire items (adult file) and examination items (exam file).
ACTIVITY_ITEMS = ["HAT2", "HAT4", "HAT6", "HAT8", "HAT10", "HAT12", "HAT14", "HAT16", "HAT18"]
CONDITION_ITEMS = ["HAC1A", "HAC1C", "HAC1D", "HAC1F", "HAC1G", "HAC1J", "HAC1K", "HAC1L",
                   "HAC1M", "HAC1N", "HAC1O"]                               # 11 conditions
ALCOHOL_ITEMS = ["MAPE1", "MAPE2", "MAPE3S", "MAPE4"]


def characteristics_iii(people: pd.DataFrame, adult: pd.DataFrame) -> pd.DataFrame:
    exam = data.read_nhanes3("exam", ["SEQN", "BMPBMI"] + ALCOHOL_ITEMS)
    d = people[["SEQN", "female", "age", "followup_years", "event"]] \
        .merge(adult, on="SEQN", how="left", validate="one_to_one") \
        .merge(exam, on="SEQN", how="left", validate="one_to_one")
    # Alcohol: drinking days in the past year / 52 x drinks per drinking day (capped at 15),
    # the same rule as NHANES 1999-2018; zero for non-drinkers.
    none = d.MAPE1.eq(2) | d.MAPE2.eq(2) | d.MAPE3S.eq(0)
    drinks = (d.MAPE3S.where(d.MAPE3S.between(1, 365)) / 52
              * d.MAPE4.where(d.MAPE4.gt(0)).clip(upper=15)).mask(none, 0.0)
    # Inactive: walked no miles and answered no to every other leisure activity in the past month.
    # HAT1S values of 9,000 or more are missing codes.
    walked = d.HAT1S.where(d.HAT1S.lt(9000))
    inactive = walked.eq(0) & d[ACTIVITY_ITEMS].eq(2).all(axis=1)
    answered = d[CONDITION_ITEMS].isin([1, 2]).all(axis=1)
    return pd.DataFrame({
        "female": d.female, "age": d.age,
        "white": d.DMARETHN.eq(1), "black": d.DMARETHN.eq(2), "mexican": d.DMARETHN.eq(3),
        "college": d.HFA8R.ge(16), "poor": d.DMPPIR.lt(1.0),
        "never_smoked": d.HAR1.eq(2), "current_smoker": d.HAR1.eq(1) & d.HAR3.eq(1),
        "drinks": drinks, "inactive": inactive, "bmi": d.BMPBMI,
        "conditions": d[CONDITION_ITEMS].eq(1).sum(axis=1).where(answered),
        "followup": d.followup_years, "death": d.event})


def table1(groups: dict) -> pd.DataFrame:
    """Table 1 cells. Percentages are of everyone of that sex (missing answers included
    in the denominator); means, SDs, medians and IQRs use those with a value."""
    def pct(s, col): return f"{100 * s[col].sum() / len(s):.1f}"
    def mean_sd(s, col, k=1): return f"{s[col].mean():.{k}f} ({s[col].std(ddof=1):.{k}f})"
    def median_iqr(s, col):
        q1, q2, q3 = s[col].dropna().quantile([.25, .5, .75])
        return f"{q2:.1f} ({q1:.1f}–{q3:.1f})"
    rows = [("Participants, n", lambda s: f"{len(s):,}"),
            ("Age, years, mean (SD)", lambda s: mean_sd(s, "age")),
            ("Non-Hispanic White, %", lambda s: pct(s, "white")),
            ("Non-Hispanic Black, %", lambda s: pct(s, "black")),
            ("Mexican American, %", lambda s: pct(s, "mexican")),
            ("College graduate, %", lambda s: pct(s, "college")),
            ("Income below poverty line, %", lambda s: pct(s, "poor")),
            ("Never smoked, %", lambda s: pct(s, "never_smoked")),
            ("Current smoker, %", lambda s: pct(s, "current_smoker")),
            ("Alcohol, drinks per week, median (IQR)", lambda s: median_iqr(s, "drinks")),
            ("Physically inactive, % [a]", lambda s: pct(s, "inactive")),
            ("Body mass index, kg/m², mean (SD)", lambda s: mean_sd(s, "bmi")),
            ("Doctor-diagnosed conditions, mean (SD) [b]", lambda s: mean_sd(s, "conditions", 2)),
            ("Follow-up, years, median (IQR)", lambda s: median_iqr(s, "followup")),
            ("Deaths, n (%) [c]", lambda s: f"{int(s.death.sum()):,} ({100 * s.death.mean():.1f})"),
            ("Deaths per 1,000 person-years", lambda s: f"{1000 * s.death.sum() / s.followup.sum():.1f}")]
    table = {"Characteristic": [label for label, _ in rows]}
    for column, s in groups.items():
        table[column] = [cell(s) for _, cell in rows]
    return pd.DataFrame(table)


# ===========================================================================
# Part D. Figure 1
# ===========================================================================
def figure1(curves: pd.DataFrame, counts: list[dict], est: pd.DataFrame) -> None:
    """Draw this run's results with the author-approved Figure 1 layout.

    Convert names only: DA/DM are the historical FI_lin/D display keys.
    Values and interval endpoints come from the arguments, never the saved
    publication inputs. Keep small input tables beside the plotted data so
    every point can be traced to this run.
    """
    from pathlib import Path

    renderer = V.publication_renderer()
    # The curve table's historical key uses an ASCII hyphen; its plotted
    # title uses an en dash. Match the data key here so no curve is omitted.
    era_names = {"modern": "NHANES 1999-2018", "iii": "NHANES III"}
    curve_rows = curves.rename(columns={"low": "ci_low", "high": "ci_high"}).to_dict("records")
    for row in curve_rows:
        row["era"] = era_names[row["era"]]
    samples = [{"era": era_names[c["era"]], "sex": "All", "n": c["n"],
                "deaths_within_10y": c["deaths_10y"]} for c in counts]
    summary = []
    measures = {"none": "none", "DA": "FI_lin", "DM": "D"}
    quantities = {"sex_only": "b0_sex_only", "contrast": "contrast", "log_hr": "b1"}
    for row in est.to_dict("records"):
        if row["distance"] not in measures or row["quantity"] not in quantities:
            continue
        summary.append({"era": row["era"], "measure": measures[row["distance"]],
                        "quantity": quantities[row["quantity"]], "observed": row["estimate"],
                        "ci_low": row["low"], "ci_high": row["high"],
                        "hr": row["hr"], "hr_low": row["hr_low"], "hr_high": row["hr_high"]})
    renderer.write_csv(OUT / "figure1_curve_input.csv", curve_rows)
    renderer.write_csv(OUT / "figure1_summary_input.csv", summary)
    repo = S.LOCAL_REPO                   # provenance paths are relative to this checkout
    renderer.render_figure1(curve_rows, samples, summary, OUT,
        source_paths={"curves": (OUT / "figure1_curve_input.csv").relative_to(repo).as_posix(),
                      "summary": (OUT / "figure1_summary_input.csv").relative_to(repo).as_posix()})


# ===========================================================================
# Numbers quoted in the Results 1 text
# ===========================================================================
def quoted_numbers(est: pd.DataFrame, groups: dict, counts: list[dict], cohorts: dict) -> pd.DataFrame:
    rows = []

    def add(key, era, value, low=np.nan, high=np.nan):
        rows.append({"key": key, "era": era, "value": value, "low": low, "high": high})

    for era, c in cohorts.items():
        add("participants", era, len(c.people))
        add("deaths", era, int(c.people.event.sum()))
        add("markers", era, len(c.markers))
    for column, s in groups.items():
        era = "modern" if column.startswith("1999") else "iii"
        sex = column.split(": ")[1].lower()
        add(f"death_rate_per_1000py_{sex}", era, 1000 * s.death.sum() / s.followup.sum())
        add(f"conditions_mean_{sex}", era, s.conditions.mean())
    for k in counts:
        add("broad_sample_n", k["era"], k["n"])
        add("broad_sample_deaths_10y", k["era"], k["deaths_10y"])
    for r in est.itertuples():
        if r.quantity == "contrast":
            add(f"{r.distance}_sex_difference_sd", r.era, r.estimate, r.low, r.high)
        elif r.quantity == "sex_only":
            add("female_hr_sex_only", r.era, r.hr, r.hr_low, r.hr_high)
        elif r.quantity == "log_hr":
            add(f"female_hr_adjusted_{r.distance}", r.era, r.hr, r.hr_low, r.hr_high)
        elif r.quantity == "per_sd":
            add(f"hr_per_sd_{r.distance}", r.era, r.hr, r.hr_low, r.hr_high)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--bootstrap", action="store_true", help="redo the 500 bootstrap draws")
    parser.add_argument("--check-draws", type=int, default=0, metavar="N",
                        help="run only draws 1..N and compare them with the stored and archived draws (seconds)")
    parser.add_argument("--workers", type=int, default=6, help="processes for the bootstrap")
    args = parser.parse_args()
    run = runs.start("s1")

    print("A. distance and mortality")
    cohorts = {era: data.load(era) for era in S.ERAS}
    point = {era: estimates(era, c.people, distances(c)) for era, c in cohorts.items()}
    stored = OUT / "distance_bootstrap.csv"
    seed = {era: S.ERAS[era]["seed"] for era in S.ERAS}
    fresh = None
    if args.check_draws:
        fresh = run_bootstrap(args.workers, range(1, args.check_draws + 1))
        run.write_check_draws(fresh, ["era", "draw"])
    if args.bootstrap:
        draws = run_bootstrap(args.workers, range(1, S.BOOTSTRAP_DRAWS + 1))
        run.save_draws(draws, stored, ["era", "draw"], seed)
    else:
        draws = run.stored_draws(stored, ["era", "draw"], seed,
                                 runs.draw_keys(S.ERAS, range(1, S.BOOTSTRAP_DRAWS + 1)), fresh)
    if draws is None:
        print("  no bootstrap draws yet: intervals left blank (run with --bootstrap)")
    est = summarise(point, draws)
    est.to_csv(OUT / "distance_estimates.csv", index=False)

    print("B. ten-year mortality by age")
    adult = data.read_nhanes3("adult", ["SEQN", "HSAGEIR", "HSSEX", "DMARETHN", "HFA8R", "DMPPIR", "HAR1",
                                        "HAR3", "HAT1S"] + ACTIVITY_ITEMS + CONDITION_ITEMS)
    curves, counts = zip(*(ten_year_curves(era, d) for era, d in
                           (("modern", broad_modern()), ("iii", broad_iii(adult)))))
    curves = pd.concat(curves)
    curves.to_csv(OUT / "ten_year_curves.csv", index=False)

    print("C. Table 1")
    chars = {"modern": characteristics_modern(cohorts["modern"].people),
             "iii": characteristics_iii(cohorts["iii"].people, adult)}
    groups = {f"{label}: {sex}": chars[era][chars[era].female.eq(f)]
              for era, label in (("modern", "1999–2018"), ("iii", "NHANES III"))
              for f, sex in ((1, "Women"), (0, "Men"))}
    V.write_table(table1(groups), OUT, "table1")

    print("D. Figure 1")
    if draws is not None:
        figure1(curves, list(counts), est)
    quoted_numbers(est, groups, list(counts), cohorts).to_csv(OUT / "numbers.csv", index=False)
    run.finish()
    print(f"done: outputs in {OUT}")


if __name__ == "__main__":
    main()
