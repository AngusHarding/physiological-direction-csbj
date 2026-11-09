"""Results 2: sex-specific references let the score register women's lower mortality.

The question
------------
The mortality score q is a weighted sum of the markers, with weights from a
Cox model of death (models.mortality_score). Each marker enters relative to
the young-healthy reference. How does q treat women and men, and how much of
the answer is the choice of reference rather than physiology?

Scope
-----
Population   the two analysis cohorts; Table 2 uses the matched NHANES 1999-2010
             cohort (6,759 participants with every Phenotypic Age input).
Endpoint     accident-censored mortality.
Estimands    Part A: under each reference (SS, SC, PC), the women-minus-men
             difference in q and the female hazard ratio with q in the model, and
             paired differences between references. Part C: young-healthy
             women-minus-men marker means in pooled SDs, with a 0.20 SD equivalence
             margin. Part D (Table 2): held-out concordance and five-year risk error
             with and without an explicit sex term.
Split        q is cross-fitted (each participant scored by a direction learned without
             their cycle, or without their NHANES III fold); Table 2 holds out one
             cycle at a time and re-estimates the reference without it.
Uncertainty  Part B: 500 participant draws per cohort with the reference held fixed
             and every direction refitted. Table 2: 500 draws resampling participants
             and reference members within cycle x sex.
Assumptions  the prepared transformations are fixed; proportional hazards.
Permitted    how the choice of reference changes what q says about women and men, and
interpretation  whether a sex term adds held-out prediction. The hazard ratios are
             descriptive associations, not causal effects.

What it makes (in outputs/s2/)
------------------------------
  figure2.svg / .png        Figure 2: a,b distance vs score; c,d the score's sex
                            contrast and the female hazard ratio under each reference
  reference_estimates.csv   Figure 2c,d numbers, with 95% intervals, and the paired
                            differences between references
  reference_bootstrap.csv   the 500 draws per cohort behind those intervals (--bootstrap)
  reference_equivalence.csv do healthy young women and men differ, marker by marker?
  table2.csv / .xlsx        Table 2: held-out prediction with and without a sex term
  prediction_bootstrap.csv  the 500 draws behind Table 2's intervals (--bootstrap)
  numbers.csv               every Results 2 number quoted in the text

Run from this folder:
    python -B s2_reference.py                  # estimates, tables, figure (a few minutes); reuses stored draws
    python -B s2_reference.py --check-draws 2  # also draws 1-2 of both bootstraps, compared with the archive
    python -B s2_reference.py --bootstrap      # redo both bootstraps (about an hour on 8 cores)

The four parts
--------------
  A  The score under three references (Fig. 2): SS, the primary, measures each
     sex from its own young-healthy mean and SD; SC keeps each sex's own mean
     but uses one SD per marker pooled across the sexes; PC pools both. For
     each, q is cross-fitted, then women and men are compared and the female
     hazard ratio is re-estimated with q in the model.
  B  The bootstrap behind Part A's intervals: participants are resampled and
     every direction is refitted; the reference itself is held fixed.
  C  Reference equivalence: are healthy young women's and men's marker means
     within 0.20 SD of each other (two one-sided tests)?
  D  Table 2: in the matched NHANES 1999-2010 cohort, each survey cycle is
     predicted from the other five, with the reference re-estimated without
     it, with and without an explicit sex term.
"""
from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
from scipy import stats                   # noqa: E402

import data                               # noqa: E402
import displays as V                      # noqa: E402
import measures as M                      # noqa: E402
import models                             # noqa: E402
import runs                               # noqa: E402
import prediction as P                    # noqa: E402
import settings as S                      # noqa: E402

OUT = S.output_folder("s2")
REF_LABEL = {"SS": "Sex-specific centre\nand SD (SS)", "SC": "Sex-specific centre,\ncommon SD (SC)",
             "PC": "Pooled centre\nand SD (PC)"}
PREDICTION_SEED = 20261023            # the Table 2 bootstrap's own seed


# ===========================================================================
# Part A. The score under three references
# ===========================================================================
def coordinates(c: data.Cohort) -> dict:
    """Marker coordinates under SS, SC and PC. SS is the stored primary coordinate set."""
    values, ref = data.transformed_values(c)
    sex = np.where(c.people.female.eq(1), 2, 1)
    X = {kind: M.reference_coordinates(values, sex, ref, kind) for kind in ("SC", "PC")}
    X["SS"] = c.z.to_numpy(float)
    return X


def scores(era: str, people: pd.DataFrame, X: dict) -> dict:
    """q under each reference.

    PC needs no fit of its own: PC coordinates are SC coordinates shifted by a
    constant for each sex, and the model contains a sex term, which absorbs
    that shift. So PC and SC have the same marker weights, and q under PC is
    the PC coordinates projected on the SC direction.
    """
    q_ss, _ = models.mortality_score(era, people, X["SS"])
    q_sc, units = models.mortality_score(era, people, X["SC"])
    q_pc, _ = models.mortality_score(era, people, X["PC"], units)
    return {"SS": q_ss, "SC": q_sc, "PC": q_pc}


def reference_estimates(era: str, people: pd.DataFrame, q: dict) -> dict:
    """For each reference: women minus men in q, and the female log hazard ratio with q.

    change = (female log HR with q) - (female log HR with sex only).
    """
    frame = people.reset_index(drop=True).copy()
    sex_only = models.cox(era, frame, ["female"])["female"]
    out = {"sex_only": sex_only}
    for kind, values in q.items():
        frame["q"] = M.standardise(values)
        log_hr = models.cox(era, frame, ["female", "q"])["female"]
        out[f"{kind}|contrast"] = models.sex_contrast(era, people.reset_index(drop=True), values)
        out[f"{kind}|log_hr"] = log_hr
        out[f"{kind}|change"] = log_hr - sex_only
    # Paired differences between references (Supplementary Table S2.2)
    for a, b in (("SC", "SS"), ("PC", "SS"), ("PC", "SC")):
        for quantity in ("contrast", "change"):
            out[f"{a}-{b}|{quantity}"] = out[f"{a}|{quantity}"] - out[f"{b}|{quantity}"]
    return out


# ===========================================================================
# Part B. Bootstrap for Part A
# ===========================================================================
_W = {}


def _start_reference_worker(era: str) -> None:
    c = data.load(era)
    _W.update(era=era, people=c.people, X=coordinates(c), strata=models.bootstrap_strata(era, c.people))


def _reference_draw(draw: int) -> dict:
    era = _W["era"]
    rows = models.bootstrap_rows(era, _W["strata"], draw)
    people = _W["people"].iloc[rows].reset_index(drop=True)
    X = {k: v[rows] for k, v in _W["X"].items()}
    return {"era": era, "draw": draw, **reference_estimates(era, people, scores(era, people, X))}


# ===========================================================================
# Part C. Reference equivalence
# ===========================================================================
def reference_equivalence(margin: float = 0.20, alpha: float = 0.05) -> pd.DataFrame:
    """For each marker, healthy young women minus men, in pooled within-sex SDs, and
    whether two one-sided t-tests (Welch) place the difference within +/- margin."""
    c = data.load_modern()
    _, ref = data.transformed_values(c)
    rows = []
    for marker in c.markers:
        men, women = ref.loc[ref.sex.eq(1), marker].to_numpy(), ref.loc[ref.sex.eq(2), marker].to_numpy()
        nm, nw = len(men), len(women)
        vm, vw = men.var(ddof=1), women.var(ddof=1)
        pooled_sd = np.sqrt(((nm - 1) * vm + (nw - 1) * vw) / (nm + nw - 2))
        difference = women.mean() - men.mean()
        se = np.sqrt(vm / nm + vw / nw)
        df = (vm / nm + vw / nw) ** 2 / ((vm / nm) ** 2 / (nm - 1) + (vw / nw) ** 2 / (nw - 1))
        p = max(stats.t.sf((difference + margin * pooled_sd) / se, df),
                stats.t.cdf((difference - margin * pooled_sd) / se, df))
        rows.append({"marker": marker, "women_minus_men_pooled_sd": difference / pooled_sd,
                     "tost_p": p, "equivalent": p < alpha})
    return pd.DataFrame(rows)


# ===========================================================================
# Part D. Table 2: held-out prediction with and without a sex term
# ===========================================================================
def matched_data() -> dict:
    """The matched cohort, its transformed markers, and the reference members of 1999-2010."""
    c = data.load_modern()
    keep = data.matched_rows(c)
    values, ref = data.transformed_values(c)
    spec = pd.read_csv(S.PREP_CONFIG / "fi_scoring_spec.csv")
    assert spec.marker.tolist() == c.markers
    ref = ref[ref.cycle.isin(data.MATCHED_CYCLES)].reset_index(drop=True)
    return {"people": c.people[keep].reset_index(drop=True), "values": values[keep].reset_index(drop=True),
            "ref": ref, "orientation": spec.orientation.to_numpy(), "side": spec.adverse_direction.to_numpy()}


def resample_matched(d: dict, draw: int) -> dict:
    """Bootstrap draw: older participants and reference members both resampled within
    cycle x sex blocks, from one random stream."""
    rng = np.random.default_rng(np.random.SeedSequence([PREDICTION_SEED, draw]))

    def rows(frame, female):
        blocks = [np.flatnonzero(frame.cycle.eq(c).to_numpy() & (female == s))
                  for c in data.MATCHED_CYCLES for s in (0, 1)]
        return np.concatenate([rng.choice(b, len(b), replace=True) for b in blocks])
    i = rows(d["people"], d["people"].female.to_numpy())
    j = rows(d["ref"], (d["ref"].sex.to_numpy() == 2).astype(int))
    return {**d, "people": d["people"].iloc[i].reset_index(drop=True),
            "values": d["values"].iloc[i].reset_index(drop=True), "ref": d["ref"].iloc[j].reset_index(drop=True)}


def prediction_metrics(d: dict) -> dict:
    """All Table 2 models, each predicted cycle by cycle; returns pooled held-out metrics."""
    people, values = d["people"], d["values"]
    sex = np.where(people.female.eq(1), 2, 1)
    held = {}
    for cycle in data.MATCHED_CYCLES:
        train = people.cycle.ne(cycle).to_numpy()
        test = ~train
        ref = d["ref"][d["ref"].cycle.ne(cycle)]            # reference learned without this cycle
        risk_sets = P.training_risk_sets(people[train])
        for kind in M.REFERENCES:
            X = M.reference_coordinates(values, sex, ref, kind)
            da = M.robust_adverse_distance(values, sex, ref, kind, d["orientation"], d["side"])
            for with_sex in (0, 1):
                q = P.learn_weights(people, X, train, risk_sets, bool(with_sex))
                for family, score in (("q", q), ("DA", da)):
                    eta = P.evaluation_model(people, score, train, risk_sets, bool(with_sex))
                    risk = P.five_year_risk(people[train].reset_index(drop=True), eta[train],
                                            people[test].reset_index(drop=True), eta[test])
                    held.setdefault(f"{family}|{kind}|{with_sex}", {})[cycle] = (eta[test], risk)
    return {model: P.evaluate(people, h) for model, h in held.items()}


def flatten(metrics: dict) -> dict:
    """The Table 2 quantities from one set of metrics."""
    out = {}
    for family in ("q", "DA"):
        for kind in M.REFERENCES:
            m0, m1 = metrics[f"{family}|{kind}|0"], metrics[f"{family}|{kind}|1"]
            key = f"{family}|{kind}"
            out[f"{key}|C_without_sex"] = m0["all"]["C"]
            out[f"{key}|C_with_sex"] = m1["all"]["C"]
            out[f"{key}|C_change"] = m1["all"]["C"] - m0["all"]["C"]
            for group in ("all", "women", "men"):
                out[f"{key}|risk_error_{group}"] = 100 * m0[group]["risk_error"]     # percentage points
                out[f"{key}|risk_error_with_sex_{group}"] = 100 * m1[group]["risk_error"]
    out["observed_risk_women"] = 100 * metrics["q|SS|0"]["women"]["observed_risk"]
    out["observed_risk_men"] = 100 * metrics["q|SS|0"]["men"]["observed_risk"]
    return out


_D = {}


def _start_prediction_worker() -> None:
    _D["data"] = matched_data()


def _prediction_draw(draw: int) -> dict:
    return {"draw": draw, **flatten(prediction_metrics(resample_matched(_D["data"], draw)))}


def table2(point: dict, draws: pd.DataFrame | None) -> pd.DataFrame:
    def interval(key, digits):
        text = f"{point[key]:.{digits}f}"
        if draws is not None:
            low, high = models.percentile_interval(draws[key])
            text += f" ({low:.{digits}f} to {high:.{digits}f})"
        return text.replace("-", "−")                      # a true minus sign, as in print
    rows = []
    for family, name in (("q", "Mortality score, q"), ("DA", "Adverse distance, DA")):
        for kind in M.REFERENCES:
            key = f"{family}|{kind}"
            rows.append({"Score": name, "Reference": "SS (primary)" if (kind == "SS" and family == "q") else kind,
                         "C, no sex term": f"{point[key + '|C_without_sex']:.4f}",
                         "C, with sex term": f"{point[key + '|C_with_sex']:.4f}",
                         "Change in C (95% interval)": interval(key + "|C_change", 4),
                         "Five-year prediction error without sex term: women, pp (95% interval)":
                             interval(key + "|risk_error_women", 2),
                         "Men, pp (95% interval)": interval(key + "|risk_error_men", 2)})
    return pd.DataFrame(rows)


# ===========================================================================
# Figure 2
# ===========================================================================
def figure2(scatter: dict, est: pd.DataFrame) -> None:
    """Draw this run's scores and intervals using the accepted Figure 2 layout.

    The drawing routine is shared with the publication rebuild. It receives
    these computed coordinates and estimates; saved publication estimates
    are not substituted. Store the renamed display inputs for traceability.
    """
    import math

    renderer = V.publication_renderer()
    participants = []
    for era in ("modern", "iii"):
        frame = scatter[era]
        # Each exported point carries the participant's own NHANES identifier (SEQN), so the
        # display file can be traced back to the cohort; never a row number (R120-F5).
        if "SEQN" not in frame or frame.SEQN.isna().any() or frame.SEQN.duplicated().any():
            raise ValueError(f"Figure 2 {era}: every point needs one unique participant identifier (SEQN)")
        for row in frame.itertuples(index=False):
            participants.append({"cohort": era, "SEQN": int(row.SEQN),
                                 "sex": int(row.female) + 1,
                                 "continuous_departure_index": float(row.DA), "q": float(row.q)})
    ledger, baseline = {}, []
    for row in est.to_dict("records"):
        if row["key"] == "sex_only":
            recipe, metric = "SS", "b0"
            baseline.append({"era": row["era"], "quantity": "b0_sex_only",
                             "hr": math.exp(float(row["estimate"]))})
        else:
            recipe, name = row["key"].split("|")
            if name not in {"contrast", "log_hr"}:
                continue
            metric = "contrast_q" if name == "contrast" else "b1_q"
        ledger[(row["era"], recipe, metric)] = {
            "cohort": row["era"], "recipe": recipe, "metric": metric,
            **{key: row[key] for key in ["estimate", "low", "high"]}}
    renderer.write_csv(OUT / "figure2_reference_input.csv", list(ledger.values()))
    renderer.write_csv(OUT / "figure2_participant_coordinates.csv", participants)
    repo = S.LOCAL_REPO                   # provenance paths are relative to this checkout
    renderer.render_figure2(participants, ledger, baseline, OUT,
        source_paths={"ledger": (OUT / "figure2_reference_input.csv").relative_to(repo).as_posix()})


# ===========================================================================
def summarise(point: dict, draws: pd.DataFrame | None, era_column: bool = True) -> pd.DataFrame:
    rows = []
    for era, values in point.items():
        for key, value in values.items():
            low = high = np.nan
            if draws is not None:
                sub = draws[draws.era.eq(era)] if era_column else draws
                low, high = models.percentile_interval(sub[key])
            rows.append({"era": era, "key": key, "estimate": value, "low": low, "high": high})
    return pd.DataFrame(rows)


def quoted_numbers(ref_est: pd.DataFrame, equivalence: pd.DataFrame, pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in ref_est.itertuples():
        kind, _, quantity = r.key.partition("|")
        if r.key == "sex_only":
            rows.append(("female_hr_sex_only", r.era, np.exp(r.estimate), np.exp(r.low), np.exp(r.high)))
        elif quantity == "contrast":
            rows.append((f"q_sex_difference_sd_{kind}", r.era, r.estimate, r.low, r.high))
        elif quantity == "log_hr":
            rows.append((f"female_hr_adjusted_q_{kind}", r.era, np.exp(r.estimate), np.exp(r.low), np.exp(r.high)))
        elif quantity == "change":
            rows.append((f"female_log_hr_change_{kind}", r.era, r.estimate, r.low, r.high))
    hb = equivalence.set_index("marker").loc["LBXHGB", "women_minus_men_pooled_sd"]
    rows += [("equivalent_markers", "modern", int(equivalence.equivalent.sum()), np.nan, np.nan),
             ("markers_differing_over_0.5sd", "modern",
              int(equivalence.women_minus_men_pooled_sd.abs().gt(.5).sum()), np.nan, np.nan),
             ("haemoglobin_difference_sd", "modern", abs(hb), np.nan, np.nan)]
    for r in pred.itertuples():
        rows.append((f"prediction_{r.key.replace('|', '_')}", "matched", r.estimate, r.low, r.high))
    return pd.DataFrame(rows, columns=["key", "era", "value", "low", "high"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--bootstrap", action="store_true", help="redo both bootstraps")
    parser.add_argument("--check-draws", type=int, default=0, metavar="N",
                        help="run only draws 1..N of both bootstraps and compare them with the stored draws")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    run = runs.start("s2")

    print("A. the score under three references")
    point, scatter = {}, {}
    for era in S.ERAS:
        c = data.load(era)
        q = scores(era, c.people, coordinates(c))
        point[era] = reference_estimates(era, c.people, q)
        scatter[era] = pd.DataFrame({"SEQN": c.people.SEQN.to_numpy(), "DA": M.adverse_distance(c.x),
                                     "q": q["SS"], "female": c.people.female, "event": c.people.event})

    def reference_draws(numbers):
        results = []
        for era in S.ERAS:
            with ProcessPoolExecutor(min(args.workers, len(numbers)), initializer=_start_reference_worker,
                                     initargs=(era,)) as pool:
                results += list(pool.map(_reference_draw, numbers, chunksize=2))
            print(f"  {era}: {len(numbers)} draws done", flush=True)
        return pd.DataFrame(results)

    def prediction_draws(numbers):
        with ProcessPoolExecutor(min(args.workers, len(numbers)), initializer=_start_prediction_worker) as pool:
            out = pd.DataFrame(list(pool.map(_prediction_draw, numbers, chunksize=2)))
        print(f"  Table 2: {len(numbers)} draws done", flush=True)
        return out

    full, check = list(range(1, S.BOOTSTRAP_DRAWS + 1)), list(range(1, args.check_draws + 1))
    ref_seed = {era: S.ERAS[era]["seed"] for era in S.ERAS}
    ref_fresh = pred_fresh = None
    if check:
        print("B. check draws for both bootstraps")
        ref_fresh, pred_fresh = reference_draws(check), prediction_draws(check)
        both = pd.concat([ref_fresh.assign(bootstrap="reference"),
                          pred_fresh.assign(bootstrap="prediction", era="matched")])
        run.write_check_draws(both, ["bootstrap", "era", "draw"])
    stored = OUT / "reference_bootstrap.csv"
    if args.bootstrap:
        print("B. bootstrap for the three references")
        run.save_draws(reference_draws(full), stored, ["era", "draw"], ref_seed)
    ref_draws = pd.read_csv(stored) if args.bootstrap else \
        run.stored_draws(stored, ["era", "draw"], ref_seed, runs.draw_keys(S.ERAS, full), ref_fresh)
    ref_est = summarise(point, ref_draws)
    ref_est.to_csv(OUT / "reference_estimates.csv", index=False)

    print("C. reference equivalence")
    equivalence = reference_equivalence()
    equivalence.to_csv(OUT / "reference_equivalence.csv", index=False)

    print("D. Table 2: held-out prediction")
    pred_point = flatten(prediction_metrics(matched_data()))
    stored = OUT / "prediction_bootstrap.csv"
    if args.bootstrap:
        run.save_draws(prediction_draws(full), stored, ["draw"], PREDICTION_SEED)
    pred_draws = pd.read_csv(stored) if args.bootstrap else \
        run.stored_draws(stored, ["draw"], PREDICTION_SEED, runs.draw_keys(None, full, era_column=False), pred_fresh)
    V.write_table(table2(pred_point, pred_draws), OUT, "table2")
    pred_est = summarise({"matched": pred_point}, pred_draws, era_column=False)

    if ref_draws is not None:
        figure2(scatter, ref_est)
    quoted_numbers(ref_est, equivalence, pred_est).to_csv(OUT / "numbers.csv", index=False)
    run.finish()
    print(f"done: outputs in {OUT}")


if __name__ == "__main__":
    main()
