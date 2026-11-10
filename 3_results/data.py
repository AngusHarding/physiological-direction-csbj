# =============================================================================
# READ ME FIRST: data.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Loads the analysis cohorts from the prepared tables: who is in
#               each cohort, their endpoint, follow-up and marker coordinates.
#               No model is fitted here.
#
# METHODS       "Participants and endpoints": "Mortality cohorts" and
#               "Endpoints"; "Statistical analysis": "Survival models"
#               (attained age, delayed entry).
#
# USED BY       Every Results script. load_modern / load_iii: Results 1-4 and 7.
#               transformed_values: Results 2. matched_rows: Results 2 (Table 2)
#               and 4 (Phenotypic Age comparison). load_common32 and transform:
#               Results 5 and 6.
#
# PSEUDOCODE
#   NHANES 1999-2018 cohort:
#     start from the 29,053 adults with a complete 38-marker panel
#     attach the linked mortality file
#     keep age >= 52 AND follow-up months > 0   (no mortality link -> dropped)
#     event = died AND cause recorded AND cause != accident
#     entry = examination age; exit = entry + months / 12
#     check: 11,497 people, 2,494 events
#   NHANES III cohort: the same, with interview age, linkage-eligible
#     participants only, and the complete 34-marker panel;
#     check: 4,474 people, 3,401 events
#   Matched 1999-2010 cohort: cycles 1999-2010 with all Phenotypic Age inputs;
#     check: 6,759 people, 2,112 events
# =============================================================================
"""Load the two analysis cohorts, one row per participant.

What a cohort is
----------------
Both cohorts are adults examined at age 52 or older who have a complete
biomarker panel and linked mortality follow-up:

  NHANES 1999-2018  11,497 participants, 2,494 deaths, 38 markers
  NHANES III        4,474 participants, 3,401 deaths, 34 markers

Each loader returns a `Cohort` with:

  people   one row per participant: SEQN (the NHANES participant number),
           female (1 = woman), age in years (NHANES 1999-2018: RIDAGEYR, age at
           screening; NHANES III: HSAGEIR, age at interview; the variables the
           reviewed analysis used), the period (survey cycle or phase),
           follow-up, and the endpoint;
  z        each marker transformed towards a normal shape and expressed in SDs
           of that participant's own sex's young-healthy reference
           ("z-score coordinates");
  x        each marker's distance beyond the young-healthy median in the
           direction that is bad for health, zero inside the healthy side
           ("adverse-excursion coordinates");
  ref_z    the z coordinates of the young-healthy reference members, by sex
           (1 = man, 2 = woman), used for the Mahalanobis distance.

The rows of people, z and x line up: row i is the same participant in all three.

The endpoint
------------
The principal endpoint is death from any cause except an identified accident:
a death counts if NHANES linked a death with a recorded cause and the cause is
not "accident" (NCHS leading-cause code 4). Accidental deaths, and the few
deaths with no recorded cause, are treated as censored. In 2015-2018 NCHS does
not release the cause categories that identify accidents, so there every death
with a cause code counts (Methods).

Time scale
----------
Every Cox model uses attained age: a participant enters the risk set at the
age above and leaves at that age plus follow-up (months from examination).

Assumptions: the prepared tables (transformations, references, coordinates)
are taken as given; the upper age limit of 79 is inherited from them.

These loaders start from tables made by the data-preparation stage (they do
not rebuild the transformations or the reference); see README.md.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import settings as S


@dataclass
class Cohort:
    era: str
    people: pd.DataFrame
    z: pd.DataFrame
    x: pd.DataFrame
    ref_z: dict
    markers: list


# -- CALCULATION: the principal endpoint -------------------------------------
# Methods, "Endpoints": "The principal endpoint counted deaths with a recorded
#   underlying cause other than an identified accident. Identified accidental
#   deaths and deaths without a recorded cause were censored at death.
#   Accidents cannot be identified in the public-use files for 2015-2018, where
#   all deaths with a recorded cause were counted"
# How: event = 1 if died (mortstat 1) AND a cause is recorded AND the cause is
#   not code 4 (accident; settings.ACCIDENT_CODE); otherwise 0. Follow-up still
#   ends at death, so a 0 for an accidental or no-cause death is a censoring at
#   death. The 2015-2018 files never carry code 4, so all deaths with a cause
#   count there.
def principal_endpoint(mortstat: pd.Series, cause: pd.Series) -> pd.Series:
    """1 for a linked death with a recorded cause other than accident, else 0."""
    died = mortstat.eq(1)
    return (died & cause.notna() & cause.ne(S.ACCIDENT_CODE)).astype(int)


# -- CALCULATION: attained-age entry and exit --------------------------------
# Methods, "Survival models": "Cox models with attained age as the time scale,
#   delayed entry at examination age (interview age in NHANES III)".
# How: entry = age at examination (interview age in NHANES III);
#   exit = entry + months of follow-up / 12.
def add_follow_up(people: pd.DataFrame, age: str) -> None:
    """Add follow-up in years and the attained-age entry and exit times (in place)."""
    people["followup_years"] = people.permth_exm / 12.0         # months to years
    people["entry_age"] = people[age].astype(float)
    people["exit_age"] = people.entry_age + people.followup_years


# ---------------------------------------------------------------------------
# NHANES 1999-2018
# ---------------------------------------------------------------------------
def mortality_1999_2018() -> pd.DataFrame:
    """The linked mortality file (all ages), with the survey cycle written as 1999-2000."""
    m = pd.read_csv(S.PREP / "01_imported/mortality.csv",
                    usecols=["SEQN", "cycle", "mortstat", "ucod_leading", "permth_exm"])
    m["SEQN"] = m.SEQN.astype("int64")
    m["cycle"] = m.cycle.astype(str).str.replace("_", "-", regex=False)
    return m


# -- COHORT: NHANES 1999-2018, ages 52-79 ------------------------------------
# Methods, "Mortality cohorts": "11,497 of 14,312 eligible participants in
#   NHANES 1999-2018 (2,494 principal-endpoint deaths)".
# How: complete-panel adults (29,053; the upper age limit of 79 is applied in
#   data preparation), restricted to age >= 52 and follow-up > 0. The
#   follow-up rule also removes the 12 participants with no mortality link
#   (11,509 -> 11,497; see the methods checklist, item 3).
def load_modern() -> Cohort:
    # The two coordinate tables cover all 29,053 adults with a complete 38-marker
    # panel, in the same participant order; the healthy reference is flagged.
    z_all = pd.read_csv(S.PREP / "09_standardised/final_38_panel_transformed_standardised_hd.csv")
    x_all = pd.read_csv(S.PREP / "11_fi/fi_signed_graded_deficits.csv")
    markers = list(z_all.columns[4:])                     # after SEQN, RIAGENDR, RIDAGEYR, is_reference
    assert z_all.SEQN.equals(x_all.SEQN) and list(x_all.columns[4:]) == markers

    ref_z = {sex: z_all.loc[z_all.is_reference.eq(1) & z_all.RIAGENDR.eq(sex), markers].to_numpy(float)
             for sex in (1, 2)}

    people = z_all[["SEQN", "RIAGENDR", "RIDAGEYR"]].merge(
        mortality_1999_2018(), on="SEQN", how="left", validate="one_to_one")
    keep = people.RIDAGEYR.ge(S.AGE_MIN) & people.permth_exm.gt(0)
    people = people[keep].reset_index(drop=True)
    people["female"] = people.RIAGENDR.eq(2).astype(int)
    people["event"] = principal_endpoint(people.mortstat, people.ucod_leading)
    add_follow_up(people, "RIDAGEYR")

    z = z_all.loc[keep.to_numpy(), markers].reset_index(drop=True)
    x = x_all.loc[keep.to_numpy(), markers].reset_index(drop=True)
    assert (len(people), people.event.sum()) == (11_497, 2_494), "NHANES 1999-2018 cohort changed"
    return Cohort("modern", people, z, x, ref_z, markers)


# ---------------------------------------------------------------------------
# NHANES III
# ---------------------------------------------------------------------------
def read_nhanes3(file: str, variables: list[str]) -> pd.DataFrame:
    """Selected NHANES III variables, from the "adult" questionnaire or the "exam" (examination).

    2_prepare_nhanes_iii/05_extract_results_variables.py copies them from the raw
    files, with NHANES III's 8- and 9-fills (88888 in a five-character field)
    already set to missing.
    """
    return pd.read_csv(S.III_PREP / f"results_{file}_variables.csv", usecols=variables)[variables]  # public_code


# NHANES III linked mortality: only participants eligible for linkage
# (eligstat 1), with the same endpoint rule as above.
def mortality_iii() -> pd.DataFrame:
    """The NCHS 2019 linked mortality file for NHANES III (participants eligible for linkage)."""
    m = pd.read_csv(S.III_PREP / "mortality_linked.csv")  # public_code: copied by 2_prepare_nhanes_iii/05
    m = m[m.eligstat.eq(1)].drop(columns="eligstat")
    m["event"] = principal_endpoint(m.mortstat, m.ucod_leading)
    return m


# -- COHORT: NHANES III, ages 52-79 ------------------------------------------
# Methods, "Mortality cohorts": "4,474 in NHANES III, drawn from 13,195 adults
#   with complete 34-marker panels (3,401 deaths)". "In NHANES III, interview
#   age served as the entry age".
# How: linked participants with follow-up > 0 and interview age >= 52, kept
#   only if they have the complete 34-marker panel.
def load_iii() -> Cohort:
    out = S.III_PREP  # public_code
    markers = pd.read_csv(out / "nhanes3_panel_markers.csv").marker.tolist()
    panel = pd.read_csv(out / "nhanes3_panel_raw.csv",
                        usecols=["SEQN", "HSAGEIR", "HSSEX", "SDPPHASE", "SDPSTRA6", "SDPPSU6"])
    z_all = pd.read_csv(out / "coordinates_hd.csv")
    x_all = pd.read_csv(out / "coordinates_fi.csv")
    x_cols = [f"signed_{m}" for m in markers]

    people = panel.merge(mortality_iii(), on="SEQN", how="inner")
    people = people[people.permth_exm.gt(0) & people.HSAGEIR.ge(S.AGE_MIN)].copy()
    people["female"] = people.HSSEX.eq(2).astype(int)
    people["age"] = people.HSAGEIR.astype(float)
    people["phase"] = people.SDPPHASE
    add_follow_up(people, "age")
    # Keep only participants with both coordinate sets (the complete 34-marker panel).
    people = people.merge(z_all[["SEQN"] + markers], on="SEQN", how="inner")
    people = people.merge(x_all[["SEQN"] + x_cols], on="SEQN", how="inner").reset_index(drop=True)
    z = people[markers].copy()
    x = people[x_cols].set_axis(markers, axis=1)
    people = people.drop(columns=markers + x_cols)

    members = pd.read_csv(out / "reference_members.csv")
    members = members[members.is_reference.eq(1)]
    ref_z = {sex: z_all.merge(members.loc[members.HSSEX.eq(sex), ["SEQN"]], on="SEQN")[markers]
             .to_numpy(float) for sex in (1, 2)}
    assert (len(people), people.event.sum()) == (4_474, 3_401), "NHANES III cohort changed"
    return Cohort("iii", people, z, x, ref_z, markers)


def load(era: str) -> Cohort:
    return load_modern() if era == "modern" else load_iii()


# ---------------------------------------------------------------------------
# Transformed marker values, before any reference is applied
# ---------------------------------------------------------------------------
def _apply_transforms(raw: pd.DataFrame, rule: pd.Series, markers: list, reciprocal_sign: int) -> pd.DataFrame:
    functions = {"identity": lambda v: v, "sqrt": np.sqrt, "cbrt": np.cbrt, "log": np.log,
                 "reciprocal": lambda v: reciprocal_sign / v}
    return pd.DataFrame({m: functions[rule[m]](raw[m].to_numpy(float)) for m in markers}, index=raw.index)


# Transformed marker values before any reference is applied. Results 2 uses
# them to rebuild coordinates under the SS, SC and PC references
# (measures.reference_coordinates) and, for Table 2, to re-estimate the
# reference without the held-out cycle.
def transformed_values(c: Cohort) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Each marker after its shape transformation (log, square root, ...) but before
    it is centred and scaled by any reference.

    Returns (cohort, reference): the cohort's rows in cohort order, and the
    young-healthy reference members with their sex (1 = man, 2 = woman) and,
    for 1999-2018, their survey cycle. A reference is then just a choice of
    centre and scale applied to these values (measures.reference_coordinates).
    """
    if c.era == "modern":
        t = pd.read_csv(S.PREP / "06_transform/final_38_panel_transformed_unstandardised_v2.csv").set_index("SEQN")
        z = pd.read_csv(S.PREP / "09_standardised/final_38_panel_transformed_standardised_hd.csv",
                        usecols=["SEQN", "is_reference"])
        members = z.loc[z.is_reference.eq(1), "SEQN"].to_numpy()
        cycles = pd.read_csv(S.PREP / "04_cohorts/eligibility_flags.csv", usecols=["SEQN", "cycle"]) \
            .set_index("SEQN").cycle
        ref = t.loc[members, c.markers].assign(sex=t.loc[members, "RIAGENDR"].astype(int).to_numpy(),
                                               cycle=cycles.loc[members].to_numpy())
    else:
        # NHANES III stores the recipe per marker; its reciprocal is taken as -1/x so
        # that, like every other transform, a larger raw value gives a larger result.
        raw = pd.read_csv(S.III_PREP / "nhanes3_panel_raw.csv").set_index("SEQN")  # public_code
        rule = pd.read_csv(S.III_PREP / "reference_parameters.csv").drop_duplicates("marker") \
            .set_index("marker")["transform"]
        t = _apply_transforms(raw, rule, c.markers, reciprocal_sign=-1)
        members = pd.read_csv(S.III_PREP / "reference_members.csv").SEQN.to_numpy()  # public_code
        ref = t.loc[members].assign(sex=raw.loc[members, "HSSEX"].astype(int).to_numpy())
    return t.loc[c.people.SEQN, c.markers].reset_index(drop=True), ref.reset_index(drop=True)


# ---------------------------------------------------------------------------
# The matched NHANES 1999-2010 cohort (prediction tests and the PhenoAge comparison)
# ---------------------------------------------------------------------------
MATCHED_CYCLES = S.CYCLES[:6]                 # 1999-2000 to 2009-2010: C-reactive protein measured


# -- INPUTS: published Phenotypic Age (Levine et al., 2018) ------------------
# The nine blood markers of Phenotypic Age, converted to the units of the
# published formula. Results 4: "published Phenotypic Age" in the matched
# NHANES 1999-2010 comparison. The formula itself is in s4_panel.py.
def phenoage_inputs(seqn: pd.Series) -> pd.DataFrame:
    """The nine blood markers of published Phenotypic Age, in its units (Levine et al., 2018).

    Albumin g/L, creatinine umol/L, glucose mmol/L, log C-reactive protein (mg/dL),
    lymphocyte %, mean cell volume fL, red cell distribution width %, alkaline
    phosphatase U/L, white cells 1000/uL. Missing where any input is missing.
    """
    p = S.PREP / "02_harmonised"
    chem = pd.read_csv(p / "standard_biochemistry.csv",
                       usecols=["SEQN", "LBXSAL", "LBXSCR", "LBXSGL", "LBXSAPSI"]).set_index("SEQN").reindex(seqn)
    cbc = pd.read_csv(p / "whole_blood_count.csv",
                      usecols=["SEQN", "LBXWBCSI", "LBXLYPCT", "LBXMCVSI", "LBXRDW"]).set_index("SEQN").reindex(seqn)
    crp = pd.read_csv(S.EXTRACTS / "candidate_crp_raw.csv").set_index("SEQN").reindex(seqn)
    # 2015-2018 reports high-sensitivity CRP in mg/L; earlier cycles in mg/dL.
    crp_mgdl = np.where(crp.assay_era.eq("LBXHSCRP_2015_2018"), crp.crp_raw / 10, crp.crp_raw)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = pd.DataFrame({"albumin": chem.LBXSAL * 10, "creatinine": chem.LBXSCR * 88.4,
                            "glucose": chem.LBXSGL * .0555, "log_crp": np.log(crp_mgdl),
                            "lymphocyte_pct": cbc.LBXLYPCT, "mcv": cbc.LBXMCVSI, "rdw": cbc.LBXRDW,
                            "alk_phos": chem.LBXSAPSI, "wbc": cbc.LBXWBCSI})
    return out.replace([np.inf, -np.inf], np.nan).reset_index(drop=True)


# -- COHORT: matched NHANES 1999-2010 ----------------------------------------
# Results 2: "a matched NHANES 1999-2010 cohort (6,759 participants; 2,112
#   deaths)", used for Table 2 and for the Phenotypic Age comparison in Results 4.
# How: the 1999-2018 cohort, restricted to 1999-2010 (when C-reactive protein
#   was measured) and to people with every Phenotypic Age input.
def matched_rows(c: Cohort) -> np.ndarray:
    """Rows of the 1999-2018 cohort in the matched cohort: 1999-2010, all PhenoAge inputs present."""
    assert c.era == "modern"
    keep = c.people.cycle.isin(MATCHED_CYCLES).to_numpy() & phenoage_inputs(c.people.SEQN).notna().all(axis=1).to_numpy()
    assert (keep.sum(), c.people.event[keep].sum()) == (6_759, 2_112), "matched cohort changed"
    return keep


# ---------------------------------------------------------------------------
# The 32 markers measured comparably in both eras (Results 5 and 6)
# ---------------------------------------------------------------------------
def common32_markers() -> pd.DataFrame:
    """inputs/common32_markers.csv: each marker's column in each era and its transformation choices."""
    return pd.read_csv(S.PIPELINE / "inputs" / "common32_markers.csv")


# -- COHORT INPUT: the 32 markers measured in both eras ----------------------
# Methods, "Distances": the bounded distance DB is computed "over the 32
#   markers measured in both eras". Results 5 (transport between eras) and
#   Results 6 (the screen) use this table, which holds all adults aged 18-79,
#   so that each script can rebuild its own reference and choose its own
#   older cohort.
def load_common32(era: str) -> pd.DataFrame:
    """Every adult of an era's prepared panel (ages 18-79) with the 32 common markers.

    Columns: SEQN, sex (1 men, 2 women), age, is_reference, stratum and psu
    (the survey's sampling strata and units, for the survey bootstrap), phase
    (NHANES III), the 32 markers under their keys, and `complete` (all 32 present).

    NHANES 1999-2018: granulocytes = white cells x (neutrophil + eosinophil +
    basophil %) / 100. A stratum number is unique only within a cycle, so the
    stratum is cycle|stratum.
    NHANES III: blood pressure rebuilt as the mean of the examination-centre
    readings 2 and 3 (zero readings set missing; values outside 50-300 systolic
    or 20-200 diastolic set missing), as continuous NHANES averages its readings.
    """
    table = common32_markers()
    if era == "modern":
        panel = pd.read_csv(S.PREP / "05_panel/final_38_biomarker_panel_raw.csv")
        cbc = pd.read_csv(S.PREP / "02_harmonised/whole_blood_count.csv",
                          usecols=["SEQN", "LBXWBCSI", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT"])
        panel = panel.merge(cbc, on="SEQN", how="left", validate="one_to_one")
        panel["granulocytes"] = panel.LBXWBCSI * (panel.LBXNEPCT + panel.LBXEOPCT + panel.LBXBAPCT) / 100.0
        design = pd.read_csv(S.PREP / "01_imported/survey_design_1999_2018.csv",
                             usecols=["SEQN", "cycle", "SDMVSTRA", "SDMVPSU"])
        panel = panel.merge(design, on="SEQN", how="left", validate="one_to_one")
        stratum = panel.cycle.astype(str) + "|" + panel.SDMVSTRA.astype(int).astype(str)
        members = pd.read_csv(S.PREP / "08_reference/healthy_reference_members.csv")
        assert members.SEQN.equals(panel.SEQN)
        out = pd.DataFrame({"SEQN": panel.SEQN.astype(int), "sex": panel.RIAGENDR.astype(int),
                            "age": panel.RIDAGEYR.astype(float), "is_reference": members.is_reference.astype(bool),
                            "stratum": stratum, "psu": stratum + "|" + panel.SDMVPSU.astype(int).astype(str)})
        for r in table.itertuples():
            out[r.key] = panel[r.modern_column].to_numpy(float)
    else:
        panel = pd.read_csv(S.III_PREP / "nhanes3_panel_raw.csv")  # public_code
        exam = read_nhanes3("exam", ["SEQN", "PEP6H1", "PEP6I1", "PEP6H3", "PEP6I3"])
        exam = panel[["SEQN"]].merge(exam, on="SEQN", how="left", validate="one_to_one")
        rebuilt = {}
        for name, pair, (low, high) in (("systolic", ["PEP6H1", "PEP6I1"], (50.0, 300.0)),
                                        ("diastolic", ["PEP6H3", "PEP6I3"], (20.0, 200.0))):
            mean = exam[pair].mask(exam[pair].eq(0)).mean(axis=1)
            rebuilt[name] = mean.mask(mean.notna() & ~mean.between(low, high))
        members = pd.read_csv(S.III_PREP / "reference_members.csv")  # public_code
        out = pd.DataFrame({"SEQN": panel.SEQN.astype(int), "sex": panel.HSSEX.astype(int),
                            "age": panel.HSAGEIR.astype(float),
                            "is_reference": panel.SEQN.isin(members.loc[members.is_reference.eq(1), "SEQN"]),
                            "stratum": panel.SDPSTRA6.astype(int).astype(str),
                            "psu": panel.SDPSTRA6.astype(int).astype(str) + "|" + panel.SDPPSU6.astype(int).astype(str),
                            "phase": panel.SDPPHASE.astype(int)})
        for r in table.itertuples():
            out[r.key] = (rebuilt[r.key] if r.key in rebuilt else panel[r.nhanes3_column]).to_numpy(float)
    out["complete"] = out[table.key.tolist()].notna().all(axis=1)
    return out


# -- CALCULATION: the shape transformations ----------------------------------
# Methods, "Panel and transformations": "one transformation per marker
#   (identity, square root, cube root, logarithm or reciprocal)". The choice
#   for each marker is made in data preparation; here it is applied. Values
#   outside a transformation's domain become missing, never shifted or clipped.
def transform(values, kind: str) -> np.ndarray:
    """A shape transformation; values outside its domain become missing (never shifted or clipped)."""
    x = np.asarray(values, dtype=float)
    out = np.full(x.shape, np.nan)
    ok = np.isfinite(x) & {"identity": True, "cbrt": True, "sqrt": x >= 0, "log": x > 0, "reciprocal": x > 0}[kind]
    f = {"identity": lambda v: v, "cbrt": np.cbrt, "sqrt": np.sqrt, "log": np.log, "reciprocal": lambda v: 1.0 / v}[kind]
    out[ok] = f(x[ok])
    return out
