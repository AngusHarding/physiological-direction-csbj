# NHANES III external replication: pre-analysis specification

**Status: SIGNED OFF 9 August 2026 by the manuscript lead, who confirmed having
read sections 7 and 8 before any mortality-linked analysis was run.** Fixed from
this point; changes after results exist require a dated amendment.

Predictions in §7 were recorded **before** any NHANES III result bearing on them
was generated. Work completed before sign-off was parsing, quality control,
transformation selection, young-healthy reference construction and coordinate
building — all outcome-blind, none touching any prediction.

**Execution order, set at sign-off.** Claim 1 (the sex–mortality gap, prediction
P1) runs first. It is the premise of claims 2–5 and the simplest analysis in the
set, which makes it a calibration check: a P1 failure would indicate a fault in
mortality linkage or cohort construction rather than in the science, and is better
found before any geometry is built on top of it.

**Feasibility was established before drafting** (§3). No result bearing on any
prediction has been computed.

---

## Correction note (9 October 2026): the marker count is 34, not 30

**Status: recorded AFTER results, at the direction of the manuscript lead. It
corrects a count. It changes no analysis, no threshold and no prediction, and
the text below it is left exactly as signed off.**

**What the text says.** This specification gives the primary marker set as 30
in five places: section 4 ("30 are directly available and up to 34 including
the examination measures"), Decision 1, Decision 3, prediction P7 ("that
cohort's 30-marker" displacement; random panels "from the same 30") and the
sign-off table.

**What the rule in Decision 1 defines, and what was analysed.** Decision 1
defines the set by a rule: the largest complete-case set that does not incur
phase-driven attrition, with the immune panel as the Coulter three-part
differential and with gamma-glutamyl transferase and globulin excluded. Applied
to the 38 markers of the main analysis, that rule removes two markers and
replaces five immune counts by three, which leaves **34**. Those 34 are the set
that was built, analysed and reported.

**Where the 30 comes from.** Thirty of the 34 are in the laboratory file. The
other four are examination measures: systolic blood pressure, diastolic blood
pressure, pulse rate and waist-to-height ratio. The count of laboratory-file
markers was written where the count of the whole set belonged.

**What agrees with 34.** Decision 3 lists pulse rate, an examination measure,
among the eight locked markers to be transported, so the eight could not have
been drawn from the 30 laboratory markers alone. The marker list used by every
analysis (`outputs/nhanes3_panel_markers.csv`) holds 34 names and is dated 9
August 2026, the day of sign-off, from the outcome-blind panel build. The
transport script states that random panels are "drawn from the same 34". The
manuscript's Methods state that 34 of the 38 markers were available.

**How to read this document.** Wherever it says 30 for the primary marker set,
for the projection that the locked panel is compared with, or for the pool that
random and native panels are drawn from, read 34. Predictions P7 and P7b were
evaluated with a 34-marker target and a 34-marker pool, and are reported that
way.

**What is not claimed.** No note made at the time records why 30 was written.
That it was a miscount is the simplest reading of the evidence above; it is an
inference. A comparison restricted to the 30 laboratory markers has not been
run.

---

## Amendment 1 (9 August 2026): functional physiological measures

**Status: recorded BEFORE the analysis was run, AFTER the primary replication was
complete.** This is an addition of scope, not a change to anything already
specified. Predictions in A1.4 were written before any model was fitted.

### A1.1 What changes

§2 lists system coverage as **out of scope, deliberately**. That decision was made
on the assumption — carried over from the development cohort — that functional
measures in this data would be too sparse to test. **That assumption was wrong.**
NHANES III measured spirometry, bone densitometry, bioimpedance and a balance task
on most of the age 52–79 cohort, and 26-year follow-up gives every one of them
thousands of events.

| Measure | System | n in panel cohort | disease deaths |
|---|---|---:|---:|
| FEV₁, FVC | Pulmonary | 4,187 | 3,163 |
| Femur neck BMD, total bone mineral content | Musculoskeletal | 4,198 | 3,179 |
| Bioimpedance resistance, reactance | Body composition | 4,275 | 3,234 |
| Tandem stand time | Balance / neuromuscular | 3,032 | 2,587 |
| ~~Simple reaction time~~ | ~~CNS~~ | 559 | 280 — **excluded, too few** |

For comparison, the development cohort tested FEV₁ on 3,409 participants with
**530** deaths across three folds, and reported it as the largest effect resting on
the thinnest evidence. This cohort has **six times the events**. Three of the four
systems were never testable in the development cohort at all, and all three are
named in that manuscript's Discussion as domains it cannot see.

### A1.2 Why this is not a specification breach, and where it must be labelled

It tests **no registered prediction**. P1–P7b concern the direction, the signal
location and the panel; none is affected by adding an exploratory analysis of
different measurements.

But it was scoped **after** the primary replication was complete and after the
power was known to be adequate, so it is **exploratory, not confirmatory**, and
must be reported under that label wherever it appears. It does not acquire the
standing of the pre-registered replication by sitting next to it.

### A1.3 Framing

The question is not "what does the routine panel miss". It is **whether a
combination of cheap, reproducible blood markers and cheap, quickly administered
physiological measures carries information that either alone does not** — which is
the practical question for healthspan trial design.

**A deployability distinction must be stated and not blurred.** Bioimpedance,
tandem stand and spirometry are inexpensive, fast, portable and require no
venepuncture. **Bone densitometry is not**: it needs a DXA scanner and is neither
cheap nor field-deployable. If BMD adds information it is scientifically
interesting; it is not part of a deployable screen, and the paper must say so.

### A1.4 Predictions, recorded before any model was fitted

**P8.** FEV₁ will retain a positive association with disease mortality conditional
on the 34-marker projection, and will improve held-out partial log-likelihood.
*Confidence: high — this replicates a development-cohort finding with six times the
events.*

**P9.** At least one of the three previously untestable systems — musculoskeletal,
body composition, balance — will add held-out information beyond the 34-marker
projection. *Confidence: high.*

**P10.** Tandem stand will add information. *Confidence: moderate-high. Balance is
entirely unrepresented in the panel and is an established mortality predictor.*

**P11.** Femur neck BMD will add information. *Confidence: moderate. Musculoskeletal
capacity is unrepresented, but BMD covaries with body size and nutritional status,
which the panel partly captures.*

**P12.** Bioimpedance will add information. *Confidence: low-moderate. Body
composition is the system most nearly represented already, through
waist-to-height ratio.*

**P13.** The functional measures will collectively add more held-out information
than the panel's own weakest blood markers contribute. *Confidence: moderate. This
is the information-density claim, and it is the one most likely to fail.*

### A1.5 Reporting commitment

Every candidate is reported whichever way it falls, with the confidence levels
above quoted as recorded. A measure that adds nothing is as informative as one that
does, because the claim being tested is about which systems the routine panel
already covers.

---

## 1. Objective

Test whether the central finding of the main manuscript replicates in an
independent cohort:

> Equal physiological distance from young-health does not imply equal biological
> consequence. The **direction** of physiological change carries information about
> mortality that **magnitude** alone misses, and it is direction rather than
> magnitude that accounts for the female survival advantage.

NHANES III (1988–1994) provides disjoint participants, a different measurement
era, different assay platforms and independent mortality linkage through 2019.

## 2. What is and is not in scope

**In scope.**

| # | Claim to be tested | Main-paper section |
|---|---|---|
| 1 | Women have lower disease mortality despite comparable or greater morbidity | §S1 |
| 2 | Dysregulation **magnitude** does **not** account for the female advantage | §S2 |
| 3 | Mortality-associated **direction** does account for most of it | §S3 |
| 4 | Most positive mortality-associated signal lies within clinical limits, much of it inside the young-healthy interval | §S1b |
| 5 | The locked reduced panel transports to an independent cohort | §S4 |

**Claims 2 and 3 must be reported together.** The paper's contribution is the
*contrast* between them. Reporting a replicated positive without the accompanying
null would misrepresent the finding.

**Out of scope, deliberately.** Trial-arm sorting (§S6), lifestyle construct
validity (§S7) and system coverage (§S8). These are applications and internal
construct checks; they would substantially increase the work without testing the
central claim. Their omission is not a judgement on their validity.

**Also out of scope.** Re-running biomarker *selection* in NHANES III. See
Decision 2.

## 3. Data and established feasibility

Source files, all public and downloaded 9 August 2026 without application or fee:
`lab.dat`, `exam.dat`, `adult.dat` with their SAS layouts, and
`NHANES_III_MORT_2019_PUBLIC.dat`.

Feasibility already confirmed:

| Quantity | Value |
|---|---|
| Adults 18–79 in the laboratory file | 16,260 |
| Complete-case on 31 available marker variables | 13,972 |
| ... with usable examination follow-up | **13,961** |
| Disease deaths in that cohort | **4,905** |
| Age ≥52 subset | **4,802** (2,464 women, 2,338 men) |
| Disease deaths, age ≥52 | **3,670** |
| Follow-up from examination | median 26.3 y, maximum 31.1 y |

Reference-state and panel dependencies confirmed present: `BMPBMI`, `BMPWAIST`,
`BMPHT`, `PEPMNK1R` (mean systolic), `PEP6DR` (pulse), the `HAC1x`
doctor-diagnosed condition battery, and smoking, alcohol and diabetes items.

**No mortality-linked analysis has been run.** The counts above are coverage
arithmetic only.

## 4. Marker set

Of the 38 markers in the main panel, **30 are directly available** and up to 34
including the examination measures. Three constraints are recorded because they
bound what can be claimed.

**C1 — the immune panel collapses from five to three.** NHANES III's automated
differential is Coulter three-part: lymphocyte, mononuclear and granulocyte
counts, at 93–95% coverage. The five-part manual differential exists but only for
a **27% subsample**, which is not usable for a complete-case panel. Neutrophils —
rank 4 in the locked panel — are therefore available only as granulocytes
(neutrophils + eosinophils + basophils).

**C2 — GGT and globulin are phase-dependent.** Coverage is 48% and 44% in Phase 1
(1988–91) against 95% in Phase 2 (1991–94). Including GGT costs approximately 21%
of the cohort. GGT is rank 5 in the locked panel.

**C3 — ACR and eGFR are derived**, as in the main analysis, from urine albumin and
creatinine, and from serum creatinine with age and sex.

**DECISION 1 — primary marker set.** The primary analysis uses the largest
complete-case set that does not incur phase-driven attrition: the 30 markers
available at ≥93% coverage, with the immune panel represented by the Coulter
three-part differential. GGT and globulin are **excluded from the primary** and
examined in a Phase 2 sensitivity.

*Rationale:* the central claim is about the geometry of a multi-marker
physiological state, not about any individual marker. Trading 21% of the cohort
for two markers weakens the test of the claim actually being made.

## 5. Design: what is re-derived and what is transported

**This is the most important decision in the specification.** Two different
questions require two different designs, and conflating them is the most common
way replications become uninterpretable.

**DECISION 2 — re-derive for the direction claim; transport for the panel claim.**

**Re-derivation (claims 1–4).** Everything is rebuilt inside NHANES III: its own
young-healthy reference under the same criteria, its own transformation selection,
its own standardisation, its own fitted mortality direction. Nothing is imported
from the main analysis except the *procedure* and the adverse-direction
specification.

*Rationale:* the claim is about human physiology, not about a fitted object. If it
holds with an independently estimated direction in a different measurement era,
that is evidence about biology. Re-derivation is also robust to the assay
differences between 1988–94 and 1999–2018, because every coordinate is referenced
to a within-cohort young-healthy state — which is precisely the portability
property the instrument is claimed to have.

**Transport (claim 5).** The locked marker identity from the main analysis is
applied unchanged. No reselection. This tests the *product*, and is what the main
manuscript states is required before applied use.

**DECISION 3 — panel transport is tested at eight markers, not ten.** Because of
C1 and C2, the primary transport test uses the eight locked markers available at
full coverage: alanine aminotransferase, red cell distribution width, mean
corpuscular volume, urine albumin-to-creatinine ratio, platelet count, albumin,
phosphorus and pulse rate. A substituted ten-marker version — granulocytes for
neutrophils, Phase 2 only for GGT — is reported as a sensitivity.

**What "transported" means, precisely.** The **marker identity** is imported and
applied unchanged; **coefficients are refitted within NHANES III**. Nothing is
reselected. Importing the coefficients as well would be a different and much
harsher test — it would require the two eras' assay calibrations to agree in
absolute terms, which the reference-relative design deliberately does not assume.

**Transport is judged against two comparators, never as a bare correlation**
(predictions P7 and P7b): a **random null** of 1,000 eight-marker panels drawn from
the 30 available, and a **native ceiling** of eight markers selected within
NHANES III itself. Both are computed in the same cross-fitting scheme as the locked
panel, so the three numbers are directly comparable.

*Rationale:* an honest eight-marker test is worth more than a strained ten-marker
one. If eight of ten reproduce the direction in an independent cohort, the claim
is supported; if the substituted ten does better or worse, that is informative
about the substitutions rather than about the panel.

## 6. Analysis specification

**Cohort.** Adults aged 18–79 at examination with complete data on the primary
marker set and positive examination follow-up. Mortality-model analyses use the
age ≥52 frame to match the main analysis; the full adult range is used for
descriptive and reference-construction steps.

**Young-healthy reference.** Rebuilt within NHANES III under the criteria of the
main analysis: under 30 years, BMI 18.5–24.9, never smoker, not a heavy drinker,
no included chronic condition or functional impairment, no diabetes or
kidney-failure flag. Estimated **separately by sex**, as in the main analysis.

**Coordinates.** HD as Mahalanobis distance from the same-sex young-healthy
centroid; sensitised FI as the mean adverse excursion, using the **unchanged**
adverse-direction specification (§S5 of the main supplement) for the markers that
carry over.

**Geometry.** Magnitude *m*, alignment *a*, projection *q* and orthogonal
displacement *R*, defined exactly as in the main analysis.

**Survival models.** Cox proportional hazards, **attained age with delayed entry**
at examination age, Efron ties, baseline hazard stratified by survey phase.

**DECISION 4 — cross-fitting unit.** The main analysis used
leave-one-survey-cycle-out across ten cycles. NHANES III has only two phases, so
that scheme is unavailable. Cross-fitting is therefore **grouped k-fold by survey
primary sampling unit**, which preserves the principle that correlated sampling
units are never split across training and test. **Leave-one-phase-out** is
reported as a coarser sensitivity that tests transport across measurement era.

*Rationale:* the purpose of grouped cross-fitting is to prevent a participant's own
outcome from influencing their coordinates and to stop the model exploiting
cluster-specific artefacts. PSU grouping serves both. Two-fold phase splitting
alone would be too coarse to estimate a 30-dimensional direction stably.

**DECISION 5 — follow-up horizon.** Full available follow-up is primary. A
20-year truncated sensitivity is reported, approximating the exposure window of
the main analysis.

*Rationale:* longer follow-up gives more events, and the attained-age time scale
handles it naturally. But the main analysis found the projection association
declines with attained age, so a 26-year median follow-up may attenuate the
per-SD estimate relative to the original. That must be anticipated rather than
discovered, and it is a reason to compare *patterns* rather than point estimates.

**Uncertainty.** Participant bootstrap stratified by sex and phase, with the
direction refitted within each replicate, as in the main analysis.

## 6a. Data quality control and transformation selection

NHANES III is fixed-width text read by column position, so a parsing error is
silent and produces plausible numbers. Quality control is therefore a specified
analysis step, not an informal precaution, and it runs **before** any coordinate is
constructed.

**DECISION 6 — QC is layout-driven and four-layered.**

**Layer 1 — parse from the published layout, not from hard-coded positions.** Column
positions are read at runtime from `lab.sas`, `exam.sas` and `adult.sas`. The layout
file is the single source of truth; the code cannot silently disagree with it.

**Layer 2 — fill codes cleared, then range assertions.** NHANES III fills
unavailable numerics with 8s or 9s across the whole field width. Each field is
cleared of its own-width 8-fill and 9-fill, then asserted to fall within a
prespecified plausible physiological interval recorded per marker in
`config/nhanes3_marker_map.csv`. A marker failing its range assertion **halts the
build**; it is not warned about and carried forward.

*Rationale:* distributional statistics detect contamination and shape problems but
are weak against pure displacement. A field parsed one character out often retains
its shape while sitting in the wrong place, and only a location assertion catches
that.

**Layer 3 — full distributional assessment on raw values.** Per marker and sex: n,
missing fraction, mean, median, range, interquartile range, skewness, excess
kurtosis, D'Agostino K² with its p-value, and a density plot. This is the direct
analogue of §S0.2 of the main supplement, run earlier in the chain.

*Rationale:* the fill-code failure mode produces a large spike far out in the tail.
It is invisible in a mean, unmistakable in a density plot, and would corrupt any
transformation fitted over it.

**Layer 4 — side-by-side against continuous NHANES.** Each marker's NHANES III
distribution is tabulated against the same marker in the 1999–2018 panel, which is
already computed. A large discrepancy is either a parsing error or a genuine
assay-era difference, and the two must be distinguished before proceeding.

*Rationale:* this is the most informative diagnostic available and costs nothing,
because both sides already exist. It also has evidential value beyond QC —
documenting where the two eras genuinely differ is what demonstrates that
reference-relative standardisation is doing real work, which is the argument the
whole replication rests on.

**DECISION 7 — transformations are re-selected within NHANES III.** Following
Decision 2, the transformation **procedure** transfers unchanged — the same
candidate set, the same unsupervised minimise-absolute-skewness rule, selected on
the NHANES III analysis cohort rather than on its young-healthy subset — but the
**choices are re-made**. No transformation is imported from the main analysis.

**Markers may select a different transformation than they did in continuous
NHANES. That is a result of re-derivation, not a defect**, and both transformation
tables are reported side by side. Silently adopting the original choices would
convert a re-derivation into a partial transport and would obscure any genuine
distributional difference between the eras.

Order of operations, fixed:

1. parse → clear fill codes → range assertions
2. distributional assessment on raw values
3. side-by-side comparison against continuous NHANES
4. transformation selection on the NHANES III analysis cohort
5. young-healthy reference, standardisation, coordinates

## 6b. Young-healthy reference: criteria and harmonisation deviations

**DECISION 8 — the reference is rebuilt in NHANES III under the main analysis
criteria, with two unavoidable deviations and one confirmed match.**

Criteria applied: age under 30; BMI 18.5–24.9; never smoker (`HAR1`, the
100-cigarette question, identical in definition to the main analysis); no
diabetes, myocardial infarction or osteoporosis; none of the chronic conditions in
the `HAC1` battery; and no functional impairment. Estimated separately by sex.

Result: **1,053 participants — 467 men, 586 women**, against 1,238 (547 M, 691 F)
in the main analysis. The ratio that governs invertibility of the sex-specific
covariance matrix is comparable: **31.0 participants per marker here against 32.6
there**. Attrition is dominated by BMI (−1,648) and smoking (−650), as in the main
analysis; every clinical criterion together removes 101.

**Deviation 1 — no kidney-failure item exists.** Continuous NHANES asks about weak
or failing kidneys (`KIQ022`); NHANES III asks only about kidney *stones*. The
criterion is dropped. **Empirically immaterial: 2 of 3,488 participants under 30
have eGFR below 60, or 0.06%.** An eGFR-based substitute was rejected because it
would tighten the eGFR reference distribution using a marker that the reference
then defines — a circularity the self-reported criterion does not have.

**Deviation 2 — alcohol quantity is not measured.** The main criterion is
`drinks_per_week <= 7`. NHANES III records alcohol only within a food-frequency
block (`HAN6HS`, `HAN6IS`, `HAN6JS`) as *times per month*, with no quantity item
anywhere in the questionnaire, so frequency and volume cannot be interconverted.
**Adopted: fewer than 30 drinking occasions per month**, the closest analogue to
seven drinks a week under a one-drink-per-occasion reading. The reference is
insensitive to the choice across the whole plausible range — 945 at <15, 995 at
<30, 1,020 at <60, 1,023 unrestricted.

**Confirmed match — asthmatics are retained**, as in the main analysis, whose
reference criteria list thirteen conditions not including asthma. Exclusion was
considered on the grounds that asthma might skew the immune subsystem, and was
tested: it removes 59 participants and moves every immune reference median and
standard deviation by under 1% (lymphocyte median 2.050 unchanged in men,
2.150→2.175 in women; granulocyte 3.825→3.850 and 4.100→4.125). Part of the reason
is mechanical — NHANES III has no separate eosinophil count, so the eosinophilic
signature is diluted into the far more abundant granulocyte total before it can
reach the reference. **The deviation was therefore rejected as buying nothing
measurable at the cost of departing from the main reference definition.**

**DECISION 9 — adverse direction transfers; orientation is recomputed.** The
`adverse_direction` field of the main analysis (§S5 of the main supplement) is a
clinical prior about physiology and transfers unchanged. The `orientation` field
does **not**: it encodes whether the *selected transformation* reverses the raw
value axis, and is −1 for reciprocal transformations only. Since transformations
are re-selected within NHANES III (Decision 7), orientation must be recomputed
from the NHANES III choice. Importing it would silently invert every marker whose
transformation differs between the eras — ten markers do.

## 7. Predictions, recorded before any result

Stated so that they can fail.

**P1.** Women will have lower disease mortality than men in NHANES III, with an
age-adjusted female hazard ratio below 1 and its confidence interval excluding 1.

**P2.** Women will **not** be less physiologically displaced than men. The adjusted
female–male difference in coordinate-vector magnitude *m* will be **≥ 0**.

**P3.** Women will be **lower in alignment** *a* and in projection *q*, both
differences negative with intervals excluding 0.

**P4.** Adjusting for magnitude will move the female hazard ratio **away from 1, or
leave it materially unchanged** — it will not attenuate the advantage.

**P5.** Adjusting for projection *q* will move the female hazard ratio
**substantially toward 1**, removing more than half the sex difference on the
log-hazard scale.

**P6.** More than half of positive mortality-associated signal will fall within
conventional clinical limits.

**P7.** The eight locked markers, **transported by identity with coefficients
refitted within NHANES III**, will reproduce that cohort's 30-marker
mortality-directed displacement **substantially better than eight-marker panels
drawn at random from the same 30** — specifically, above the **95th percentile of a
1,000-draw random null**.

**P7b.** The locked eight will achieve **at least 90% of the reproduction attained
by eight markers selected natively within NHANES III**, that native selection being
the practical ceiling for a panel of that size in this cohort.

*No absolute correlation threshold is predicted, deliberately.* An earlier draft of
this specification predicted "above 0.80", which is not a test. Three problems: a
bare correlation has no null, so it cannot distinguish a panel that carries the
direction from one that merely has eight dimensions; the NHANES III comparator is
**30** markers rather than the 38 of the main analysis, making reproduction a
strictly easier task and any threshold calibrated on the original too lenient; and
the main analysis never evaluated a size-8 panel, so no anchor exists. The
comparator-based form above is falsifiable in a way the threshold was not.

*What the two comparators mean.* The random null answers "would any eight markers
do this?" The native ceiling answers "how much is lost by importing a panel rather
than fitting one here?" A transported panel that sits near the ceiling and far above
the null has genuinely transported. One that sits inside the null has not, whatever
its absolute correlation.

**Predictions deliberately not made.** No prediction is recorded about the *size*
of the projection association per SD, because the longer follow-up and older
attained-age distribution are expected to attenuate it (Decision 5). Nor about
which markers would be selected, since selection is not re-run.

## 8. What would falsify the original finding

Recorded so the conclusion is not negotiated after the fact.

**Clear failure to replicate**, requiring the main manuscript to be qualified:

- P3 fails — women are **not** lower on alignment or projection; **or**
- P5 fails — adjusting for projection does not move the female hazard ratio
  toward 1; **or**
- P4 reverses — adjusting for magnitude *does* substantially attenuate the female
  advantage, which would mean magnitude explains the gap after all.

**Partial replication**, to be reported as such rather than framed as success:

- direction attenuates the gap but by materially less than the ~80% seen in the
  main analysis;
- the pattern holds in one representation but not the other;
- P6 fails while P3–P5 hold;
- **P7 holds but P7b fails** — the locked panel beats the random null but falls
  well short of a natively selected panel. That would mean these markers carry
  real directional information but are not the right eight for this cohort, which
  is a qualified transport rather than a failed one.

**Failure of panel transport specifically**, to be reported as such without
qualifying claims 1–4:

- **P7 fails** — the locked eight fall inside the random null. The panel does not
  transport. This would **not** falsify the direction finding, which is tested by
  P3–P5 on re-derived coordinates and does not depend on the panel. It would
  falsify the claim that this particular marker set is portable, and the main
  manuscript's panel claims would need qualifying accordingly.

*Recorded because the two are separable and will be reported separately.* The
direction claim is about physiology and is re-derived; the panel claim is about a
product and is transported. Conflating a panel-transport failure with a failure of
the central finding would misrepresent both.

**Not falsification:** a smaller per-SD hazard ratio for *q*, a smaller absolute
sex difference, or a lower reproduction correlation for the reduced panel. These
are expected consequences of longer follow-up, a smaller cohort, a reduced marker
set and a different assay era.

## 9. Reporting commitment

The result is reported whichever way it falls, with the predictions above quoted
as recorded. If the finding does not replicate, that is reported in the main
manuscript and the claims are qualified accordingly. This specification is cited
in the write-up so that its date precedes the analysis.

## 10. Sign-off

| Item | Setting | Confirmed |
|---|---|---|
| Scope | claims 1–5; trial sorting, lifestyle and system coverage excluded | |
| Primary marker set | 30 markers, Coulter three-part immune, GGT/globulin excluded | |
| Design | re-derive for direction, transport for panel | |
| Panel transport | eight markers primary, substituted ten as sensitivity | |
| Cross-fitting | grouped k-fold by PSU; leave-one-phase-out as sensitivity | |
| Follow-up | full follow-up primary, 20-year truncation as sensitivity | |
| Quality control | four layers; range failure halts the build | |
| Transformations | re-selected within NHANES III; both tables reported | |
| Predictions | recorded in §7 before any result | |
| Falsification | criteria fixed in §8 | |

---

**Open questions for the manuscript lead before sign-off**

1. **Cross-fitting unit (Decision 4).** PSU-grouped k-fold is my recommendation,
   but leave-one-phase-out is the closer analogue of the original cycle-grouping.
   Phase splitting tests transport across measurement era, which is arguably the
   more meaningful axis — at the cost of only two folds. Confirm the primary.
2. ~~**Age frame.**~~ **RESOLVED 9 August 2026 by the manuscript lead: age ≥52 is
   the primary frame**, matching the main analysis (4,802 participants, 3,670
   disease deaths). If that proves statistically inadequate, widening the frame is
   permitted — but any widening must be **reported as a deviation from this
   specification, with the ≥52 result shown alongside**, not substituted for it.
   Widening is a response to insufficient precision, never to an unwelcome result.
3. **Where this lands.** A supplementary section of the main manuscript, a new
   main-text Results section, or a separate short paper. This affects how much is
   written up but not what is run.
