# Sensitised linear frailty index methods

**Version 2 — amended 5 August 2026.** See "Amendment 1" immediately below.
Version 1 remains readable in version control; its Step 6 results are marked
superseded in place.

## Amendment 1 (5 August 2026): the index is made genuinely linear

**Change.** Both grading constants are removed from Step 6. The dead zone changes
from `a = 1.0` to `a = 0`, and the cap `c = 3.0` is removed. The deficit severity
becomes the adverse excursion itself:

```
version 1:  d = clip((e - 1) / (3 - 1), 0, 1)
version 2:  d = e
```

No other step, parameter, reference population, directionality assignment or
transformation is altered. Step 6 now has **no free parameters**.

**Trigger — a specification defect, not a performance observation.** This
document is titled *"Sensitised **linear** frailty index methods"*. Version 1's
Step 6 assigned exactly zero to any adverse excursion within one robust-scale
unit of the same-sex healthy median, and assigned exactly one to every excursion
beyond three. An index that is flat at both ends is piecewise-thresholded, not
linear. The contradiction is verifiable from this document alone, without
reference to any result.

**Why the dead zone had to go.** It set 58.3% of all participant-marker values to
exactly zero in the age 52–79 analysis cohort, precisely in the band nearest the
healthy reference where mortality signal is present. It also redistributed
deficit burden away from mortality-informative markers: alanine
aminotransferase, the highest-ranked marker in the reduced panel, received only
62% of the burden share a linear index gives it, while eGFR and systolic blood
pressure gained roughly a third more than their linear share.

**Why the cap had to go.** The cap flattened 7.0% of all values overall, but its
effect was concentrated in exactly the markers where the population has drifted
furthest from young-healthy physiology: waist-to-height ratio **68.7%**, systolic
blood pressure 29.5%, eGFR 26.3%, glycated haemoglobin 22.3%, urine
albumin-creatinine ratio 14.4%. For waist-to-height the cap was not outlier
control; it made the marker effectively binary for two thirds of the cohort while
that marker carried the largest single share of deficit burden.

This also contradicted a decision already taken elsewhere in this pipeline.
`TRIMMING_ASSESSMENT.md` rejected deleting extreme values because excluded
participants were **substantially enriched for mortality** (death-rate ratio
2.42), concluding that deletion "would remove informative pathological
observations rather than random noise". The cap flattened those same
mortality-enriched values one step later. Retaining extremes and then capping
them are inconsistent positions; the trimming assessment contains the evidence
that settles which is correct.

No ceiling can repair this. Most 52–79-year-olds are genuinely beyond three
young-healthy scale units on waist-to-height ratio, so any bound — 3, 6,
min–max, or a clinical limit — saturates. Only a scale that never saturates
preserves discrimination.

**Why unboundedness is acceptable.** The `[0,1]` bound bought comparability with
published frailty indices, which this instrument does not have in any case: its
deficits are continuous and reference-relative rather than binary and
clinically-defined. The feared outlier domination does not occur, because
averaging across 38 markers dilutes any single extreme: the uncapped index has
median 1.069, interquartile range 0.899–1.264 and maximum 2.67, so the most
extreme participant sits 2.5 times the median rather than at some degenerate
value.

**Interpretation of the version-2 index.** `FI` is the mean adverse excursion
across the 38 markers, in same-sex young-healthy robust-scale units:
*"on average, this participant sits X young-healthy scale units into the adverse
tail."* Standardisation is complete before aggregation — Step 4 divides by the
sex-specific robust scale — so markers are already on a common scale and the mean
is well defined. Contributions are not equal, however: an unweighted mean
implicitly weights each marker by how far the population has drifted from
young-healthy, and waist-to-height ratio accordingly carries 10.2% of total
burden and explains 23.9% of index variance. That is the classical
deficit-accumulation bargain, in which every deficit counts once. Mortality-
weighted variants are used where explicitly stated and are a different
instrument.

**Why this is not tuning.** The rationale is definitional and independent of any
outcome, and the change cannot rescue a preferred result: the HD branch does not
pass through deficit coding, so the primary reduced-panel finding is computed
entirely from HD coordinates and does not move. Only the FI branch and the
cross-representation comparisons are affected.

**Prespecified sensitivity.** The version-1 bounded grading
`d = clip((e - a)/(c - a), 0, 1)` is retained in code as a parameterised option
and reported as a sensitivity, so the effect of bounding can be measured rather
than argued.

**Consequences.** All FI-derived artefacts must be regenerated from Step 6
onward, and every downstream module that hashes
`outputs/11_fi/fi_signed_graded_deficits.csv` — analysis modules 04, 05 and 06 —
must be rerun. Fingerprints recorded under version 1 are superseded and must not
be reused. Deficit severities and signed deficits are no longer bounded to
`[0,1]` and `[-1,1]` respectively; downstream bound assertions are updated
accordingly.

## Scope

This document records the prospective construction of the sensitised linear
frailty index (FI) from the locked 38-biomarker transformed panel. The FI is a
separate scoring branch from health distance (HD). It must not use the
HD-standardised biomarker panel, HD means, standard deviations, covariance
matrices or precision matrices.

## Step 1 — Lock the FI scoring contract

Before estimating any reference statistics or participant deficits, construct
and validate one canonical FI scoring specification containing exactly one row
for each of the 38 locked biomarkers. The scoring specification must contain:

- `marker`: biomarker identifier in the locked panel order;
- `selected_transform`: transformation already selected on the full cohort;
- `orientation`: `-1` for a transformation that reverses the raw-value axis
  and `+1` otherwise;
- `adverse_direction`: one of `high_only`, `low_only` or `two_sided`;
- `physiological_system`: the locked subsystem assignment; and
- a brief physiological justification or source for the adverse direction.

The transformation and orientation must come from
`config/final_38_transformations_v2.csv`. The adverse direction and subsystem
assignment may be brought forward from `config/biomarker_registry.csv`, but
must be reviewed and copied into a dedicated FI specification. The registry's
historical `frozen_transform` field must not be used because transformation
selection has subsequently been updated.

Adverse direction is defined on the original physiological measurement scale.
It must not be re-estimated because a transformation changed. During later FI
standardisation, multiplication by the locked transformation orientation will
restore transformed coordinates to the raw-value direction before the
adverse-tail rule is applied.

Directionality must not be selected or revised using associations with age,
sex, mortality, disease outcomes or the observed transformed-data
distribution. Any correction made during review must be justified from the
biomarker's physiological interpretation and documented before FI scores are
calculated.

### Step 1 validation gate

Proceed only when:

1. the transformed panel and transformation recipe are explicitly identified
   and fingerprinted;
2. the participant identifier columns and all 38 biomarker columns occur in
   their expected order;
3. the scoring specification has exactly 38 unique markers and matches the
   transformed panel one-to-one;
4. every transformation and orientation agrees with the locked transformation
   recipe;
5. every marker has one valid adverse-direction category and one subsystem;
6. every reciprocal transformation has orientation `-1`, and all other locked
   transformations have orientation `+1`; and
7. no FI reference parameters, deficits or scores have yet been estimated.

The resulting scoring specification becomes an immutable input to all later FI
steps. Any subsequent change requires a documented version change and complete
recalculation of the FI branch.

### Step 1 decision record

The canonical specification is `config/fi_scoring_spec.csv`. Transformation
and orientation were copied from the locked version-2 recipe. Directionality
was reconciled against the previously reviewed FI directionality specification.

Two differences were found between that reviewed specification and the current
general biomarker registry. The registry labelled diastolic blood pressure and
pulse rate `two_sided`; the reviewed FI specification labelled both
`high_only`. The FI specification retains `high_only`: hypertension is the
prespecified diastolic disease axis, while low resting pulse is clinically
ambiguous without fitness and rhythm information. These are physiological,
outcome-blind decisions and are recorded in the specification's
`direction_basis` field.

The executable validation gate is
`src/10_validate/validate_fi_scoring_spec.py`, with regression tests in
`tests/test_fi_scoring_spec.py`. Passing this gate completes Step 1 but does not
calculate reference parameters or participant FI values.

Certified SHA-256 fingerprints for Step 1 are:

- transformed panel: `E2F404E78414F9C33B9DCB8DA5BC905C3AD25F4C47472E406F13564CB9243EDC`;
- transformation recipe: `5D5E4189589304E417C88F99236C24ADE9BE2574BB46183D5962F62C998E6569`;
- FI scoring specification: `2E096673A2AC23932A566FF93DC48C8755AD77CC86F635DF30DDDC13C10304FD`.

## Step 2 — Lock the FI healthy-reference population

Use the previously reconstructed availability-aware healthy reference without
reselection or FI-specific exclusions. Join
`outputs/08_reference/healthy_reference_members.csv` to the locked transformed
panel by `SEQN`, preserving the transformed panel's participant order.

The FI reference is the same 1,238-person young-healthy reference used to
construct the HD branch: 547 men and 691 women. Membership was defined before
FI calculation using age under 30, BMI 18.5-24.9, absence of affirmative
chronic disease, never smoking, no more than seven alcoholic drinks per week,
and absence of affirmative functional impairment. Survey non-administration
was handled using the locked availability-aware rules documented in
`HEALTHY_REFERENCE.md`.

This step creates an FI-specific copy of membership, not a new reference
definition. No participant may be added or removed based on biomarker values,
future FI values, age trends, sex differences, mortality or any other outcome.

### Step 2 validation gate

Proceed only when:

1. the source reference fingerprint is
   `9CFDA42A07F3DB275547A7014EFF7C11EBA8BD2F59DF93332B30CB218FEC15B2`;
2. all 29,053 transformed participants match the reference artifact one-to-one;
3. participant sex and age agree exactly between artifacts;
4. participant order is unchanged after the join;
5. the FI reference contains exactly 547 men and 691 women, all younger than
   30 years;
6. all 38 transformed biomarkers are finite and complete in both the reference
   and full FI cohort; and
7. no medians, robust scales, deficits or FI scores have yet been estimated.

The executable is `src/11_fi/lock_fi_reference.py`. It writes
`outputs/11_fi/fi_reference_membership.csv` and
`outputs/11_fi/fi_reference_validation.csv`. Regression tests are in
`tests/test_fi_reference.py`.

Step 2 passed with 14,916 men and 14,137 women in the full cohort. The reference
contains 547 men and 691 women, spans ages 18-29 in both sexes, and has complete
finite data for all 38 biomarkers. Certified SHA-256 fingerprints are:

- FI reference membership: `6D2D62F7C390A1CFB1A88F6C7FC6746F0670D1D7926DECE4AFCFD977C98ADC6A`;
- FI reference validation summary: `CA195BD8BC77FBEFB41A3E69C063E75E976BD44D4C9D75649BB242A7C126360D`.

## Step 3 — Estimate sex-specific robust reference parameters

For each of the 38 transformed biomarkers, estimate reference parameters
separately in the 547 male and 691 female FI reference members. On the locked
transformed scale calculate the first quartile, median and third quartile using
the linear sample-quantile convention, then calculate:

\[
IQR_{js}=Q_{3,js}-Q_{1,js}
\]

and:

\[
S_{js}=\frac{IQR_{js}}{1.349}.
\]

The divisor 1.349 makes the IQR-based scale comparable to a standard deviation
under a normal distribution while retaining robustness to tails. Parameters
are fitted only in same-sex healthy-reference members. No survey weights,
trimming, winsorisation, imputation or outcome information are used.

Save one row per marker and sex with the locked transformation, orientation,
adverse direction and subsystem alongside Q1, median, Q3, IQR, robust scale,
reference N, quantile method and divisor.

### Step 3 validation gate

Proceed only when:

1. there are exactly 76 unique marker-sex parameter rows;
2. marker order and scoring metadata agree exactly with the locked FI scoring
   specification within each sex;
3. every male parameter uses 547 reference observations and every female
   parameter uses 691;
4. every quartile, IQR and robust scale is finite;
5. every IQR and robust scale is strictly positive;
6. selected parameters reproduce a direct calculation from the transformed
   reference values; and
7. no participant coordinates, deficits or FI scores have yet been calculated.

The executable is `src/11_fi/estimate_fi_reference_parameters.py`. It writes
`outputs/11_fi/fi_reference_robust_parameters.csv` and
`outputs/11_fi/fi_reference_parameter_checks.csv`. Regression tests are in
`tests/test_fi_reference_parameters.py`.

Step 3 passed for all 76 marker-sex combinations. No robust scale was zero,
negative or non-finite. The observed robust-scale range was
0.004677-32.616753 in men and 0.005528-29.651594 in women. Certified SHA-256
fingerprints are:

- FI robust reference parameters: `77852C65DD36434E2806BE8F387BE0E71AD51A5512282C7D6FA0E62688839219`;
- FI reference parameter checks: `2C88BC7F557ACB55571B3400DBC520FC4883931B4BC1D6B696573F59AD270D1D`.

## Step 4 — Calculate raw-oriented robust coordinates

For participant `i`, marker `j` and sex `s`, calculate:

\[
z^{(r)}_{ijs}
=
m_j\frac{T_j(x_{ij})-M_{js}}{S_{js}},
\]

where `T_j(x_ij)` is the already transformed value, `M_js` is the same-sex
healthy-reference median, `S_js` is the same-sex healthy-reference robust scale
from Step 3, and `m_j` is the locked transformation orientation. Orientation is
`-1` for reciprocal transformations and `+1` otherwise.

Multiplication by orientation restores the original raw-value axis, so a
positive coordinate always means that the participant's raw biomarker value is
above the same-sex healthy median. These coordinates are intentionally
unbounded. Do not apply adverse direction, take absolute values, impose the
dead zone, cap values, or aggregate markers during this step.

### Step 4 validation gate

Proceed only when:

1. all 29,053 participants and 38 markers are retained in locked order;
2. all 1,104,014 participant-marker coordinates are finite;
3. the same-sex reference median is zero for every coordinate to numerical
   precision;
4. the same-sex reference IQR is 1.349 for every coordinate to numerical
   precision;
5. direct calculations reproduce saved coordinates for both orientation-
   preserving and reciprocal transformations; and
6. no adverse-tail rule, dead zone, cap, deficit or FI aggregation has yet been
   applied.

The executable is `src/11_fi/build_fi_robust_coordinates.py`. It writes
`outputs/11_fi/fi_raw_oriented_robust_coordinates.csv` and
`outputs/11_fi/fi_robust_coordinate_checks.csv`. Regression tests are in
`tests/test_fi_robust_coordinates.py`.

Step 4 passed for all 1,104,014 participant-marker coordinates. No coordinate
was missing or non-finite. All same-sex reference medians were exactly zero at
the saved numerical precision; the maximum absolute reference-IQR deviation
from 1.349 was `2.29e-14`. Certified SHA-256 fingerprints are:

- raw-oriented robust coordinates: `4F5E71E2C50A4307858AEA66C0AE5CBDAE064CAFB5AB66BED090E5194B7A028C`;
- robust-coordinate checks: `587BD596B31DF040C529C9159A0CC549EF6DD3E10D05C4F8D07FB23C9AC4BB3B`.

## Step 5 — Map coordinates to adverse excursions

Apply the locked adverse-direction rule to each raw-oriented robust coordinate
without grading or capping it. For coordinate `z`:

\[
e_{ij}=
\begin{cases}
\max(0,z^{(r)}_{ij}), & \text{high-only},\\
\max(0,-z^{(r)}_{ij}), & \text{low-only},\\
|z^{(r)}_{ij}|, & \text{two-sided}.
\end{cases}
\]

The excursion is a non-negative distance into the biomarker's adverse region,
measured in healthy-reference robust-scale units. A high-only marker therefore
has zero excursion below its healthy median, and a low-only marker has zero
excursion above its healthy median. Both tails contribute for a two-sided
marker.

Do not apply the one-scale dead zone or three-scale cap during this step. Values
greater than three must remain present so the later grading operation can be
validated independently.

### Step 5 validation gate

Proceed only when:

1. all 29,053 participants and 38 markers remain in locked order;
2. every excursion is finite and non-negative;
3. every protective-side value of a one-sided marker maps exactly to zero;
4. high-only, low-only and two-sided marker examples reproduce their defining
   equations exactly;
5. excursions greater than three remain uncapped; and
6. no deficit severity, signed deficit or FI score has yet been calculated.

The executable is `src/11_fi/build_fi_adverse_excursions.py`. It writes
`outputs/11_fi/fi_adverse_excursions.csv` and
`outputs/11_fi/fi_adverse_excursion_checks.csv`. Regression tests are in
`tests/test_fi_adverse_excursions.py`.

Step 5 passed for all 1,104,014 participant-marker excursions. The direction
rules reproduced exactly, including zero excursion for 181,556 protective-side
high-only values and 39,229 protective-side low-only values. Excursions remain
uncapped, with maxima of 17.29, 10.36 and 30.04 robust-scale units among
high-only, low-only and two-sided markers respectively. Certified SHA-256
fingerprints are:

- adverse excursions: `7D8FC71445DB1AEA89F426CF38ADDB873FC71F52B33D75A3229C5B21978548B1`;
- adverse-excursion checks: `917FE738D62999CC1153F4D28CE10BC7C5BC19DEDD7B02474BE84C94E0FF6FBB`.

## Step 6 — Grade adverse excursions into bounded deficits

**Version 2 rule (amended 5 August 2026; see Amendment 1).** The absolute deficit
severity is the adverse excursion itself, with no dead zone and no cap:

\[
d_{ij} = e_{ij}.
\]

Severity is therefore exactly proportional to displacement into the adverse tail,
measured in same-sex young-healthy robust-scale units, across the whole range. A
participant twice as far into the adverse tail scores twice the deficit, at every
magnitude. No excursion scores zero unless it is zero, and no two different
excursions receive the same severity.

Step 6 has **no free parameters**. There is nothing to tune by marker
distribution, reference-floor behaviour, age, sex or outcome, because there is
nothing to tune. Severity is non-negative and unbounded above, and contains no
directional sign.

**Superseded version 1 rule:** `d = clip((e-1)/(3-1), 0, 1)`, with dead zone
`a=1.0` and cap `c=3.0`. Removed as a specification defect; see Amendment 1. The
bounded form is retained in code as a parameterised sensitivity.

### Step 6 validation gate

Proceed only when:

1. all 29,053 participants and 38 markers remain in locked order;
2. every severity is finite and non-negative;
3. **severity equals the excursion exactly** for every participant-marker value;
4. an excursion of exactly zero maps to zero, and **every non-zero excursion maps
   to a strictly positive severity** (no dead zone);
5. **no two distinct excursions map to the same severity** (no cap), verified by
   confirming that the maximum severity equals the maximum excursion rather than
   a ceiling;
6. synthetic values, including values far beyond the former cap, reproduce the
   rule exactly; and
7. no signed deficit or participant-level FI score has yet been calculated.

The executable is `src/11_fi/build_fi_deficit_severities.py`. It writes
`outputs/11_fi/fi_absolute_deficit_severities.csv` and
`outputs/11_fi/fi_deficit_severity_checks.csv`. Regression tests are in
`tests/test_fi_deficit_severities.py`.

**SUPERSEDED — version 1 result, retained for the record.** Step 6 passed for all
1,104,014 participant-marker values. Of these, 704,288 were in the dead zone and
scored zero, 346,475 were graded linearly, and 53,251 were capped at one.
Certified SHA-256 fingerprints under version 1 were:

- absolute deficit severities: `231195D2721EC85CBCEFA42E1763F7AF3A6BCE1081D49768C0DA72290014AFF5`;
- deficit-severity checks: `9D67CF3A8FB646B90F634D5031DE3523211C2512D318578DF3A26C7653DF7CF4`.

These fingerprints are void for version 2 and must not be reused.

**Version 2 result (5 August 2026).** Step 6 passed for all 1,104,014
participant-marker values. Of these, 273,633 had zero adverse excursion and
scored zero, and 830,381 received a strictly positive severity. Maximum severity
is 30.044, exactly equal to the maximum excursion, confirming no ceiling; the
smallest positive severity is 1.03e-15, confirming no dead zone. Severities take
126,581 distinct values, against 49,663 under version 1.

Relative to version 1, **430,655 values were freed from the dead zone and 53,251
from the cap** — 484,000 in total, 43.8% of the matrix, which previously carried
no information and now do. The proportion of exactly-zero values falls from 63.8%
to 24.8%; the remainder are genuine, being participants on the non-adverse side
of a one-sided marker.

The two indices correlate 0.962 (Spearman 0.964), but **18.3% of participants
move more than ten percentile points in FI rank**, so this is a substantive
reordering rather than a rescaling. Scalar FI moves from median 0.145
(range 0–0.518) to median 0.907 (range 0.249–2.673).

Certified SHA-256 fingerprints under version 2 are:

- absolute deficit severities: `F3BA6949DE3E517F5733D85ED3DBC4194CA0EA93622E4A1821710B388834874A`;
- deficit-severity checks: `4FF5EFD90EC428279D13D139E2EC604EEAD58E3C8EBC65C3573C0D1D1A8370EE`;
- signed graded deficits: `49DEC494701ADDA311EE5E2393453EA5D8A0A4A8CD1A27CCDFD98F25DF179F26`;
- FI scores: `22A3F08F195DB035E29716883B5E5C5F771AD84004B838752EE7A718496AF26D`.

Version-1 artefacts are retained for comparison under
`outputs/11_fi/_v1_superseded/`.

## Step 7 — Construct signed graded deficits

Restore adverse-tail sign to each bounded absolute deficit severity. Let
`d_ij` be the Step-6 severity and `z_ij` the Step-4 raw-oriented coordinate:

\[
g_{ij}=
\begin{cases}
+d_{ij}, & \text{high-only},\\
-d_{ij}, & \text{low-only},\\
\operatorname{sign}(z^{(r)}_{ij})d_{ij}, & \text{two-sided}.
\end{cases}
\]

Thus positive values identify high-side deficits and negative values identify
low-side deficits. For two-sided markers the sign identifies which side of the
healthy median generated the deficit. The magnitude is unchanged:
`|g_ij| = d_ij`. Values in the dead zone remain exactly zero.

The signed matrix preserves information required for later directional-vector
analyses. It is not itself the scalar FI, which is calculated from absolute
deficit severities rather than from signed values that could cancel.

### Step 7 validation gate

Proceed only when:

1. all 29,053 participants and 38 markers remain in locked order;
2. every signed deficit is finite and bounded to `[-1,1]`;
3. `|g_ij|` equals the Step-6 absolute severity exactly for every value;
4. high-only nonzero deficits are positive and low-only nonzero deficits are
   negative;
5. two-sided deficit signs agree with their raw-oriented coordinates; and
6. no participant-level or subsystem FI aggregation has yet been performed.

The executable is `src/11_fi/build_fi_signed_deficits.py`. It writes
`outputs/11_fi/fi_signed_graded_deficits.csv` and
`outputs/11_fi/fi_signed_deficit_checks.csv`. Regression tests are in
`tests/test_fi_signed_deficits.py`.

Step 7 passed for all 1,104,014 participant-marker values. Magnitudes were
preserved exactly. The signed matrix contains 182,466 positive high-only
deficits, 41,884 negative low-only deficits, and 175,376 nonzero two-sided
deficits (77,806 below-median and 97,570 above-median). Certified SHA-256
fingerprints are:

- signed graded deficits: `5AE0B636C1B5CFDBCF287F4D05B97956B3284AABC173E4612977CCB1CC267AA4`;
- signed-deficit checks: `394579CEE91493546B0539385C86CAAF77B9F4D4C520F46F36C521AE2358C3DB`.

## Step 8 — Aggregate participant-level FI outputs

Calculate the overall sensitised linear FI as the unweighted arithmetic mean
of the 38 absolute deficit severities:

\[
FI_i=\frac{1}{38}\sum_{j=1}^{38}d_{ij}.
\]

Calculate each physiological-system FI as the unweighted mean of the absolute
deficit severities assigned to that system. This produces ten subsystem scores:
calcium, cardiovascular, immune, kidney, liver, metabolism, lipids, platelet,
red blood cell and electrolytes. Systems contribute in proportion to their
number of markers in the overall FI; no system-level reweighting is applied.

Also retain the signed-deficit vector norm:

\[
\lVert g_i\rVert_2=\sqrt{\sum_{j=1}^{38}g_{ij}^2}.
\]

The vector norm measures total directional-vector magnitude but is not the FI.
The scalar FI and every subsystem FI must use absolute severities so opposing
directions cannot cancel.

Report healthy-reference distributions and the unadjusted association with age
as outcome-blind diagnostics only. They must not be used to revise reference
membership, directionality, transformations, dead-zone or cap parameters.

### Step 8 validation gate

1. All 29,053 participants remain in locked order.
2. The overall FI reproduces the row mean of all 38 absolute deficits.
3. Every subsystem FI reproduces the row mean of its locked marker set.
4. The vector norm reproduces the Euclidean norm of the signed deficit vector.
5. Overall and subsystem FI values are finite and bounded to `[0,1]`.
6. The vector norm is finite and bounded to `[0,sqrt(38)]`.
7. Healthy-reference and age diagnostics are reported without score tuning.

The executable is `src/11_fi/calculate_fi_scores.py`. It writes
`outputs/11_fi/fi_scores.csv`, `outputs/11_fi/fi_score_checks.csv` and
`outputs/11_fi/fi_score_diagnostics.csv`. Regression tests are in
`tests/test_fi_scores.py`.

Step 8 passed for all 29,053 participants. Every aggregation identity
reproduced exactly. Overall FI ranged from 0 to 0.517581 and the maximum signed
vector norm was 4.104879. Median FI was approximately 0.0666 in each healthy
reference sex, compared with 0.1785 in men and 0.1937 in women aged 52 or older.
The outcome-blind unadjusted FI slope was +0.02116 per decade. These descriptive
diagnostics did not alter the scoring definition. Certified SHA-256
fingerprints are:

- participant FI scores: `FA4A0527141D8A0D343E5D9DF59D832B172FA67B3D7FB7B7758F2AC296946A5D`;
- FI score checks: `B32919FB67C28823DFA35D195444F297AED86805BDE05956FABF3A0ECAC60924`;
- FI score diagnostics: `407372FB52C727D3050D5034BCB82D9EAD33BE4A65C44EE4454BE0C7BCAB7A99`.

## Post-construction quality-control diagnostics

Three outcome-blind diagnostics were run after the FI definition and scores
were locked. These checks assess marker dominance, survey-cycle stability and
sex/age behaviour within physiological systems. They do not select markers,
alter transformations, tune FI thresholds or revise scores.

### QC 1 — Marker contribution and cap-frequency audit

For each marker, calculate mean and median deficit severity, nonzero fraction,
cap fraction, share of the summed 38-marker deficit burden, correlation with
overall FI and burden rank. Repeat this in the full cohort, each sex, and each
sex at age 52 or older. A conservative review flag is raised when any marker
accounts for more than 10% of total deficit burden, is capped in more than 25%
of participants, or has absolute correlation above 0.80 with the overall FI.
These thresholds are diagnostic flags rather than exclusion rules.

Waist-to-height ratio (`WHtR`) was the only flagged marker. In the full cohort
its mean severity was 0.670, 79.6% of values were nonzero, 53.2% were capped,
and it accounted for 11.47% of total deficit burden. Its correlation with the
overall FI was 0.502, so it does not determine participant ranking on its own.
The pattern was present in both sexes and was stronger at age 52 or older.

This is not a newly introduced reconstruction artefact. In the historical FI,
`WHtR` had mean severity 0.616, 40.3% capping, 10.87% of total burden and
correlation 0.516 with overall FI. The updated transformation increased its
capping but preserved the pre-existing contribution pattern. The primary FI
therefore retains `WHtR`; a future leave-one-marker-out sensitivity analysis
may quantify its influence on substantive associations.

No other marker accounted for more than 5.38% of total burden. The next largest
contributors were eGFR (5.38%) and systolic pressure (5.25%).

### QC 2 — Survey-cycle distribution stability

FI means were directly age-standardised within sex and survey cycle. Five age
bands (18-29, 30-39, 40-51, 52-64 and 65-79 years) were weighted to the pooled
age distribution of the corresponding sex. Cycle deviations were expressed in
units of the full-cohort FI standard deviation. An absolute deviation above
0.20 SD was flagged for review.

Nineteen of the twenty sex-cycle cells were within 0.20 SD. Women in 2017-2018
were the sole marginal flag: age-standardised mean FI 0.17354 versus the pooled
female mean 0.15730, a difference of +0.01624 FI units or +0.211 SD. Men in
2017-2018 showed the same direction but remained below the flag threshold
(+0.169 SD). The next largest deviation was -0.138 SD in men in 2003-2004.

The 2017-2018 elevation is also present in the historical FI (+0.220 SD in
women and +0.183 SD in men), while the historical 2003-2004 male deviation was
-0.145 SD. This close directional and numerical agreement argues against an
error introduced by the rebuilt transformations or scoring pipeline. The cycle
indicator should nevertheless remain available for downstream sensitivity
adjustment because a real secular or residual measurement-era effect cannot be
excluded by this diagnostic.

### QC 3 — Sex- and age-stratified subsystem review

For every subsystem, report mean, median, quartiles, zero fraction and one
fraction by sex and the same five age bands. Also calculate the unadjusted
within-sex linear slope per decade and Spearman age correlation. Negative age
slopes are flagged for interpretation but are not assumed to be errors: a
subsystem FI can be non-monotonic even when the overall FI increases with age.

Eighteen of twenty subsystem-sex slopes were positive. The two negative slopes
were small and occurred in women:

- red-blood-cell FI: -0.00221 per decade, Spearman correlation +0.024;
- platelet FI: -0.00119 per decade, Spearman correlation -0.020.

Both reflect shallow non-monotonicity rather than a strong inverse age trend.
Female red-blood-cell mean FI rose from 0.104 at ages 18-29 to 0.131 at ages
40-51, then was 0.095 at 52-64 and 0.099 at 65-79. Female platelet mean FI was
0.107, 0.113, 0.120, 0.103 and 0.104 across the five age bands. The same weak
patterns were already present historically: female red-blood-cell slope
-0.00202 per decade and platelet slope approximately -0.00009. They therefore
do not indicate a new orientation or reconstruction error.

The strongest positive age slopes occurred in metabolism (+0.0894 per decade
in men and +0.0857 in women), followed by kidney (+0.0505 and +0.0625) and
cardiovascular FI (+0.0287 and +0.0458). This is physiologically coherent with
the overall positive FI-age gradient.

### QC implementation and outputs

The executable is `src/11_fi/audit_fi_quality.py`, with regression tests in
`tests/test_fi_quality_audit.py`. It writes:

- `outputs/11_fi/fi_qc_marker_contributions.csv`;
- `outputs/11_fi/fi_qc_cycle_stability.csv`;
- `outputs/11_fi/fi_qc_subsystem_age_sex.csv`; and
- `outputs/11_fi/fi_qc_subsystem_age_slopes.csv`.

All three diagnostics passed structural and formula checks. The findings call
for transparent sensitivity analysis of `WHtR` and survey cycle in later
inferential models, but they do not identify a reason to revise or reject the
locked FI.

Certified SHA-256 fingerprints are:

- marker contributions: `24B59AB77F29457E19BE3720FD073AF6D8B6A5AEB9881C09CE7B43673C616608`;
- cycle stability: `7ADFE5757AF41EBEAABA15AFAD3BD3DDE01449FC7829CE8C2CE79B60412F2CCF`;
- subsystem age-sex summary: `2A2010801CC2215B7722278DD98BA974CBEDFC13B6130C6E07EA4471F61249DA`;
- subsystem age slopes: `B4A7026EAB40ADB60CFA568E47176CB8EF55A6B099652E0AD1CC45E9D0390628`.
