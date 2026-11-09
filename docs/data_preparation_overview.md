# Data preparation: living reconstruction record

> **Status:** Chunk 1 verified.  
> **Rule:** a stage is marked verified only after its output is compared with the certified release-candidate artefact.

> **Disease update (2026-08-02):** MCQ, DIQ, KIQ, BPQ, and OSQ disease
> responses have been independently reconstructed and audited. See
> [`docs/health_conditions.md`](../docs/health_conditions.md).

## Behavioural domain — smoking

### Scientific purpose

Construct a consistent three-level cigarette-smoking status across NHANES
1999–2000 through 2017–2018 for descriptive and sensitivity analyses.

### Decision rule

| Lifetime smoking (`SMQ020`) | Current smoking (`SMQ040`) | Status |
|---|---|---|
| No | Structurally skipped | Never |
| Yes | Every day or some days | Current |
| Yes | Not at all | Former |
| Missing, refused, unknown, or inconsistent | Any | Unclassified |

The raw import retains all 64,874 questionnaire participants and both source
codes. Harmonisation classifies 55,913 participants: 30,781 Never, 13,616 Former,
and 11,516 Current. The 8,961 unresolved records remain visible with explicit
classification reasons in the audit output.

### Code and outputs

- `src/01_import/import_smoking.py`
- `src/02_harmonise/harmonise_smoking.py`
- `src/10_validate/compare_smoking.py`
- `outputs/01_imported/smoking_raw.csv`
- `outputs/02_harmonised/smoking_status.csv`
- `outputs/02_harmonised/smoking_classification_audit.csv`
- `outputs/02_harmonised/smoking_by_cycle.csv`
- `provenance/smoking_participant_comparison.csv`

### Certification

The independently reconstructed output exactly matches the certified historical
output for participant membership and classification: 55,913 of 55,913 records
matched, with zero new-only records, zero historical-only records, and zero status
disagreements.

### Status

Verified.

## Behavioural domain — physical activity

### Instruments and estimands

The 1999–2006 questionnaire used `PAQ180`, a four-level description of usual
daily activity. Code 1 was classified sedentary and codes 2–4 nonsedentary. The
2007–2018 GPAQ instrument records work, transport, and recreational activity.
Domain-specific weekly MET-minutes were calculated as days multiplied by minutes
and conventional moderate (4 MET) or vigorous (8 MET) weights.

The all-cycle `physically_active` field is deliberately broad: it contrasts
sedentary with any recorded activity. It is not equivalent to meeting exercise
guidelines. `healthy_exercise` identifies at least 600 MET-min/week and is available
only in the GPAQ era.

### Historical reconstruction and correction

The raw import contains 90,816 participants. The exact historical algorithm yields
75,034 rows and 62,041 observed all-cycle activity flags. Every participant and all
four retained fields exactly match the certified historical dataframe.

The audit identified two GPAQ exclusions: 11,522 respondents lacked valid responses
to all five activity-domain gates, and 4,260 otherwise valid respondents had a
calculated dose above 10,000 MET-min/week. The historical algorithm dropped the
second group completely. This is unnecessarily destructive for the binary exposure:
even if the exact dose is implausible, these respondents clearly reported activity.

The corrected output therefore retains those participants, caps only the continuous
MET dose at 10,000, and preserves `physically_active=1` and
`healthy_exercise=1`. Both versions remain available:

- `outputs/02_harmonised/physical_activity.csv` — exact historical reconstruction;
- `outputs/02_harmonised/physical_activity_inclusion_preserving.csv` — corrected candidate.

The corrected dataframe also preserves `met_minutes_week_raw`, the numeric
`met_extreme_above_10000_flag` (1=yes, 0=no, missing before 2007), and the text field
`met_quality_flag`. Flagged values are
not declared false; they are marked `extreme_self_report_review_recommended` because
they may reflect overlapping domains, misunderstanding, or implausible duration.
Participant-level and distributional translations are saved in
`physical_activity_extreme_met_audit.csv` and
`physical_activity_extreme_met_summary.csv`.

In the Section 3 complete-case cohort, the correction restores 2,189 participants
(1,601 men and 588 women), increasing the descriptive cohort from 27,061 to 29,250.
Because restoration is sex-imbalanced, the corrected cohort must be used to rerun
the mortality displays rather than assuming the historical estimates remain valid.

### Code and audit outputs

- `src/01_import/import_physical_activity.py`
- `src/02_harmonise/harmonise_physical_activity.py`
- `src/10_validate/compare_physical_activity.py`
- `outputs/01_imported/physical_activity_raw.csv`
- `outputs/02_harmonised/physical_activity_by_cycle.csv`
- `outputs/02_harmonised/physical_activity_gpaq_gate_audit.csv`
- `provenance/physical_activity_participant_comparison.csv`

### Status

Historical reconstruction verified; inclusion-preserving correction selected for
downstream reassessment.

## Behavioural domain — Healthy Eating Index 2020

### Historical method identified

The historical R script called `heiscore::score()` with `method="simple"`,
`component="total score"`, and each cycle from 2005–2006 through 2017–2018.
Inspection of frozen `heiscore` 0.1.3/0.1.4 source and `heiscore.data` 0.0.1
established that the adult scoring standards are HEI-2020. The exposure is therefore
named `HEI_2020_total`, rather than the ambiguous historical name `score`.

### Reconstruction

The exact CRAN source archives were frozen and SHA-256 hashed. Python reads the
embedded cycle-specific FPED-derived tables without executing R. Participants require
a nonzero two-day dietary weight and the package's default demographic eligibility.
Day-one and day-two amounts are combined and converted to energy-density ratios for
all components except the unsaturated-to-saturated fatty-acid ratio. Thirteen component
scores are calculated using the frozen HEI-2020 standards and summed to 0–100.

The reconstruction retains participant, cycle, dietary weights, demographics,
component amounts, component ratios, component scores, the package-compatible total,
and a strict complete-component total. Eight historical participants aged 2–3 had
one unresolved fatty-acid component that contributed zero to the package-compatible
total; this edge behavior is explicit and does not occur in the adult study cohort.

### Certification

The reconstructed output contains 49,621 participants, exactly matching historical
membership. Every total matches within `1e-10`; there are zero value disagreements
and the maximum absolute floating-point difference is `7.11e-14`.

### Outputs

- `outputs/01_imported/hei_package_fped_inputs.csv`
- `outputs/02_harmonised/hei_2020_component_scores.csv`
- `outputs/02_harmonised/hei_2020_total.csv`
- `outputs/02_harmonised/hei_2020_scoring_audit.csv`
- `outputs/02_harmonised/hei_2020_by_cycle.csv`
- `config/hei_2020_scoring_standards.csv`
- `provenance/hei_package_manifest.csv`
- `provenance/hei_participant_comparison.csv`

### Certification boundary

Historical package reproduction is verified. The embedded FPED-derived tables remain
a frozen third-party derived input. Independent comparison with official USDA FPED
cycle releases is the remaining higher-level provenance task; this distinction does
not prevent exact reproduction of the manuscript exposure.

### Status

Verified for historical package reconstruction; official-FPED source validation pending.

## Behavioural domain — alcohol

### Scientific purpose

Construct average past-year alcoholic drinks per week on one scale across all ten
cycles while retaining abstainers and distinguishing them from missing responses.

### Harmonisation

For 1999–2016, `ALQ120Q` and `ALQ120U` were converted to drinking days per year
and multiplied by `ALQ130`, the average drinks per drinking day. For 2017–2018,
`ALQ121` categories were mapped to prespecified annual-frequency midpoints and
combined with `ALQ130`. The final calculation was:

\[
\text{drinks/week}=\frac{\min(\text{drinking days/year},365)}{52}
\times\min(\text{drinks/drinking day},15).
\]

Participants routed past the consumption questions by a valid nondrinker screener,
or reporting zero drinking occasions, were assigned zero. Refused, unknown, and
insufficient responses remained missing.

The 15-drink quantity cap harmonizes older cycles with the 2017–2018 questionnaire
top-code. It affected 152 participants. The maximum derived exposure was 105.29
drinks/week with the cap and 396 drinks/week without it. Both versions are retained
in `alcohol_derivation_audit.csv`; downstream analyses use the capped harmonized
measure and can explicitly test the uncapped alternative.

### Results and certification

The raw import contains 53,559 participants; 48,671 received a drinks-per-week
estimate, of whom 39.4% were assigned zero. Participant membership and all numeric
values exactly matched the certified historical output: zero new-only records, zero
historical-only records, zero value disagreements, and maximum absolute difference
zero.

### Code and outputs

- `src/01_import/import_alcohol.py`
- `src/02_harmonise/harmonise_alcohol.py`
- `src/10_validate/compare_alcohol.py`
- `outputs/01_imported/alcohol_raw.csv`
- `outputs/02_harmonised/alcohol_drinks_per_week.csv`
- `outputs/02_harmonised/alcohol_derivation_audit.csv`
- `outputs/02_harmonised/alcohol_by_cycle.csv`
- `outputs/02_harmonised/alcohol_cap_sensitivity.csv`
- `provenance/alcohol_participant_comparison.csv`

### Status

Verified.

## Stage 0 — Authority and scope

### Scientific purpose

Freeze the data sources and comparison targets before rebuilding participant cohorts or biomarker representations.

### Scope

- NHANES 1999–2000 through 2017–2018: ten intact two-year cycles.
- Participant key: `SEQN`.
- Raw inputs remain outside this public-facing directory.
- Frozen release-candidate outputs are comparison targets, not inputs to reconstruction.

### Current authority

| Object | Role |
|---|---|
| Raw DEMO and HSQ XPT files | Immutable reconstruction inputs |
| `release_candidate_v1/analysis/00_data_foundation/output/demographics_dataframe.csv` | Certified eligible-cohort comparison target |
| `release_candidate_v1/analysis/00_data_foundation/code/cleaning/build_eligible_sample.py` | Historical certified implementation evidence |

### Provenance result

- Twenty selected raw files were hashed: ten DEMO and ten HSQ files.
- Their combined size is 52,290,880 bytes.
- All 20 selected hashes occur in the certified release-candidate raw dependency manifest.
- The duplicate 2007–2008 DEMO copy is byte-identical; the demographics path is selected.

### Status

Verified for Chunk 1.

## Stage 1 — Demographic and recent-illness import

### Scientific purpose

Create minimally interpreted, cycle-labelled demographic and recent-illness tables before applying study eligibility rules.

### Inputs

- DEMO files: `SEQN`, `RIDAGEYR`, `RIAGENDR`, `RIDEXPRG`.
- HSQ files: `SEQN`, `HSQ500`, `HSQ510`, `HSQ520`.

### Import rules

- Restore literal SAS zero values that pandas represents as extremely small positive numbers.
- Retain source questionnaire codes at this boundary.
- Add the two-year NHANES cycle explicitly.
- Do not exclude participants during import.
- Fail on duplicated or missing `SEQN` values within the stacked table.

### Code

- `src/01_import/nhanes_reader.py`
- `src/01_import/import_demographics.py`
- `src/01_import/import_recent_illness.py`

### Outputs

- `outputs/01_imported/demographics.csv`
- `outputs/01_imported/recent_illness.csv`
- cycle-level import summaries in `outputs/01_imported/`

### Reconciliation result

| Imported table | Participants | Cycles | Missing SEQN | Duplicate SEQN rows |
|---|---:|---:|---:|---:|
| Demographics | 101,316 | 10 | 0 | 0 |
| Recent illness | 93,240 | 10 | 0 | 0 |

Of the 93,240 HSQ participants, 84,803 have non-missing responses to all three required recent-illness questions.

### Status

Verified.

## Stage 2 — Kidney measurements and derivations

### Scientific structure

The kidney block has three separate roles. Urine albumin and urine creatinine
are measured inputs used to calculate ACR. Assay-corrected serum creatinine from
the verified biochemistry block is combined with age and sex to calculate eGFR.
The questionnaire report of weak or failing kidneys is only a healthy-reference
exclusion flag and is never scored as a biomarker.

### Raw sources and harmonisation

Twenty raw files are used: ten urine laboratory files and ten kidney
questionnaire files. All twenty SHA-256 hashes match the certified dependency
manifest. The urine import contains 81,913 participants; 79,437 have both
`URXUMA` and `URXUCR`. Neither concentration is clinically clipped because
extreme albuminuria and urine-concentration variation can be genuine.

`URXUMA` is selected instead of `URXUMS`. It is present in all ten cycles and
its units, ug/mL, are numerically identical to mg/L. Using `URXUMS` would have
silently lost the first three cycles.

### ACR

\[
\mathrm{ACR\ (mg/g)} = \frac{\mathrm{URXUMA\ (mg/L)}}
{\mathrm{URXUCR\ (mg/dL)}}\times100.
\]

The factor of 100 is required by the unit conversion. The component output uses
the historical name `Alb_Cr`; the final panel uses `ACR`. All 79,437 calculated
values exactly reproduce the certified derived-urine output.

### eGFR

The rebuild uses the 2021 race-free CKD-EPI creatinine equation:

\[
\mathrm{eGFR}=142\min(S_{cr}/\kappa,1)^\alpha
\max(S_{cr}/\kappa,1)^{-1.200}0.9938^{Age}
\times(1.012\text{ if female}),
\]

where `kappa` is 0.9 for males and 0.7 for females, and `alpha` is -0.302
for males and -0.241 for females. Race is not used. Serum creatinine is the
already calibrated `LBXSCR`, avoiding duplicated assay correction.

There are 32,421 eligible participants with complete ACR and eGFR dependencies.
For every participant subsequently retained in the frozen locked sample, both
ACR and eGFR agree exactly with the certified values.

### Kidney-failure flag

`KIQ020` in 1999–2000 and `KIQ022` thereafter are harmonised to
`kidney_fail={1 yes, 0 no}`. Refused, unknown and missing responses are discarded
rather than coded as healthy. Of 55,081 questionnaire records, 54,965 have a
valid binary response; 1,767 report weak or failing kidneys.

### Reconciliation

| Check | Result |
|---|---|
| Raw source hashes | 20 of 20 certified matches |
| Complete urine pairs and ACR values | 79,437 |
| Valid kidney-failure flags | 54,965 |
| Eligible participants with ACR and eGFR | 32,421 |
| Urine measured output | Exact match |
| Urine derived output | Exact match |
| Kidney flag output | Exact match |
| Locked-sample ACR and eGFR | Exact match |
| Maximum numerical difference | 0 |
| Repository tests | 22 passed |

Exact rules are in `config/kidney_harmonisation.csv` and comparisons are in
`provenance/kidney_comparison.csv`. Later transformation and standardisation
remain separate stages.

### Status

Verified.

## Stage 3 — Complete blood count

### Scientific purpose and layer boundary

CBC is reconstructed in two deliberately separate tables. The measured layer
contains twelve analyser-reported fields. The derived layer contains absolute
differential counts, calculated red-cell indices and inflammatory ratios. This
separation makes deterministic relationships visible and prevents calculated
features from masquerading as independent measurements.

Five measured fields enter the final 38-marker panel: platelet count, mean
platelet volume, haemoglobin, mean corpuscular volume and red-cell distribution
width. Calculated MCHC (`LBXMC`) is also a final marker. The remaining measured
fields are retained because they are required to reproduce the derived layer.

### Import and cross-cycle handling

Ten CBC files, renamed from `LAB25.XPT` through `CBC_J.XPT`, contribute 92,431
unique participants. No variable rename or assay equation is applied. Cycle is
retained because NHANES documents analyser and anticoagulant changes that may
produce cycle effects later handled during standardisation.

### Error screening

The frozen broad screens retain platelet 1–999, MCV 50–130 fL, RDW 10–45%, and
each differential percentage 0–100 inclusive. WBC, RBC, haemoglobin and mean
platelet volume remain unbounded because clinically extreme values can reflect
real disease. Observed removals were five platelet values, one MCV value and two
RDW values; no differential percentages were outside their bounds.

Complete-casing across the twelve measured fields yields 82,733 participants.

### Derived equations

Absolute differential counts use:

\[
\mathrm{absolute\ count}=\mathrm{WBC}\times\mathrm{percentage}/100.
\]

The red-cell calculations are:

\[
\mathrm{HCT}=\mathrm{RBC}\times\mathrm{MCV}/10,
\quad \mathrm{MCH}=10\times\mathrm{HGB}/\mathrm{RBC},
\quad \mathrm{MCHC}=100\times\mathrm{HGB}/\mathrm{HCT}.
\]

The documented `LBDLYMNO > 10` guard made 66 derived lymphocyte counts missing.
NLR, PLR, MLR, SII, SIRI and plateletcrit are exploratory derived features;
division by zero is converted to missing. LMR is omitted because it is exactly
the reciprocal of MLR.

### Reconciliation

| Check | Result |
|---|---|
| Raw imported participants | 92,431 |
| Complete measured CBC participants | 82,733 |
| Raw source hashes | 10 of 10 certified matches |
| Measured output | Exact match |
| Derived output | Exact match |
| Maximum numerical difference | 0 |
| Repository tests | 18 passed |

Rules are in `config/cbc_harmonisation.csv`; exact comparisons are in
`provenance/cbc_comparison.csv`. Transformation and standardisation occur later.

### Status

Verified.

## Stage 4 — Cardiovascular examination

### Scientific purpose

Reconstruct resting systolic pressure, diastolic pressure and 60-second pulse
rate from the blood-pressure examination component for all ten cycles. Pulse
regularity travels with this component solely for later healthy-reference
screening; it is not one of the 38 scored biomarkers.

### Raw inputs and audit boundary

The import preserves `BPXSY1`–`BPXSY4`, `BPXDI1`–`BPXDI4`, `BPXPLS` and
`BPXPULS` for 96,766 unique examination participants. Thus, a reviewer can see
both the unused first pressure and every reading that contributes to a mean.
The ten source paths are in `config/cardiovascular_sources.csv`; their sizes and
SHA-256 hashes are in `provenance/cardiovascular_raw_manifest.csv`.

### Derivation

Reading 1 is excluded as the settling or white-coat measurement, following the
frozen analysis. Literal zero pressure readings are converted to missing before
aggregation. The two pressure biomarkers are then:

\[
\mathrm{systolic}=\operatorname{mean}(\mathrm{BPXSY2},\mathrm{BPXSY3},\mathrm{BPXSY4})
\]

and

\[
\mathrm{diastolic}=\operatorname{mean}(\mathrm{BPXDI2},\mathrm{BPXDI3},\mathrm{BPXDI4}),
\]

using every available non-missing reading. A missing fourth reading therefore
does not discard otherwise usable measurements. The reading-count audit table
shows how many values contributed to each participant's two means.

The observed zero removal affected no systolic readings and 1,874 diastolic
readings: 814 second, 906 third and 154 fourth readings.

### Pulse regularity

Raw `BPXPULS` codes are recoded as `1 regular -> 0` and `2 irregular -> 1`.
Missing values remain missing rather than being incorrectly labelled irregular.
There were 3,498 missing regularity values at raw import, but none among the
70,830 participants who passed the three-biomarker complete-case rule. This is
an observed result, not a requirement imposed on `BPXPULS`.

### Error screening and complete cases

Only broad impossibility bounds are applied here: systolic 60–290 mmHg,
diastolic 30–150 mmHg and pulse 30–220 bpm, all inclusive. These are not
healthy ranges, and clinically abnormal but plausible observations remain.
The screen removed 314 derived diastolic values and two pulse values; no derived
systolic values were removed.

Complete-casing requires the three true biomarkers only. It does not require
`BPXPULS`. All rules and their rationales are enumerated in
`config/cardiovascular_harmonisation.csv`.

### Code and outputs

- `src/01_import/import_cardiovascular.py`
- `src/02_harmonise/harmonise_cardiovascular.py`
- `src/10_validate/compare_cardiovascular.py`
- `outputs/01_imported/cardiovascular_raw.csv`
- `outputs/02_harmonised/cardiovascular.csv`
- change, reading-count and cycle summaries beside those outputs

### Reconciliation result

| Check | Result |
|---|---|
| Raw imported participants | 96,766 |
| Complete harmonised participants | 70,830 |
| Cycles | 10 |
| Raw XPT hashes matching certified manifest | 10 of 10 |
| Participant membership | Exact match |
| Column names and order | Exact match |
| Maximum numerical difference | 0 |
| Automated repository tests | 14 passed |

The reconstructed table is exactly identical to the frozen cardiovascular
output. The later skew transformation and standardisation stages have not been
applied here and remain separately auditable work.

### Status

Verified.

## Stage 5 — Eligible adult cohort

### Scientific purpose

Identify adults whose examination biomarkers are not interpreted during pregnancy or a recorded recent acute illness.

### Eligibility rules

A participant enters the certified eligible cohort when:

\[
18 \leq \mathrm{RIDAGEYR} < 80,
\]

`RIDEXPRG != 1`, and all three recent-illness responses equal `2` (“No”):

- `HSQ500` — recent head or chest cold;
- `HSQ510` — recent stomach or intestinal illness;
- `HSQ520` — recent flu, pneumonia or ear infection.

`SEQN`, sex and age must be non-missing. Missing HSQ responses do not count as “No” and therefore cannot enter the eligible cohort.

The historical pregnancy rule is reproduced exactly: only `RIDEXPRG == 1` is excluded. Structural missingness for participants who were not asked the pregnancy question is not treated as pregnancy.

### Audit fields

The participant-level eligibility table records:

- age eligibility;
- not-pregnant flag;
- recent-illness response completeness;
- no-recent-acute-illness flag;
- demographic completeness;
- final eligibility;
- primary exclusion reason.

### Certified target

| Cohort | Total | Men | Women |
|---|---:|---:|---:|
| Eligible adult base | 34,785 | 17,924 | 16,861 |

### Code

- `src/04_eligibility/build_eligible_cohort.py`

### Outputs

- `outputs/04_cohorts/eligibility_flags.csv`
- `outputs/04_cohorts/eligible_adult_cohort.csv`
- `outputs/04_cohorts/eligibility_cascade.csv`
- `outputs/04_cohorts/eligibility_by_cycle_and_sex.csv`

### Status

Verified.

### Eligibility cascade

| Step | Total retained | Men | Women | Removed at step |
|---|---:|---:|---:|---:|
| Raw demographics | 101,316 | 49,893 | 51,423 | — |
| Required demographics complete | 101,316 | 49,893 | 51,423 | 0 |
| Matched recent-illness record | 93,240 | 45,863 | 47,377 | 8,076 |
| All three illness responses complete | 84,803 | 42,028 | 42,775 | 8,437 |
| Age 18–79 | 48,123 | 23,570 | 24,553 | 36,680 |
| No recorded recent acute illness | 35,820 | 17,924 | 17,896 | 12,303 |
| Not pregnant | 34,785 | 17,924 | 16,861 | 1,035 |

The cascade order reproduces the historical builder. “Removed at step” is conditional on passing every preceding step.

### Exact comparison with the certified cohort

| Check | Result |
|---|---|
| Row count | 34,785 = 34,785 |
| `SEQN` set | Exact match |
| Sex values | 0 differences |
| Age values | 0 differences |
| Source row order | Exact match |
| CSV SHA-256 | Exact match: `1555929F237EAA9BBC2F08E133EC74C2A67222449707AD4DA19F837B80533728` |

The rebuilt cohort is byte-for-byte identical to the certified release-candidate cohort.

### Tests

Six automated tests passed, covering participant-key validation, the SAS-zero threshold, certified cohort counts, unique/non-missing `SEQN`, and agreement between eligibility flags and the materialised cohort.

### Questions and discrepancies

- No scientific discrepancy was found.
- A path-root error in the new eligibility entry point was detected before cohort creation, corrected, and recorded in the audit log. It did not affect scientific results.

## Stage 6 — Standard biochemistry harmonisation

### Scientific purpose

Construct a comparable participant-level biochemistry table across ten NHANES cycles despite changes in source filenames, variable prefixes and laboratory analysers.

This stage concerns assay comparability. It does not apply the later skew-reducing transformations used for HD/FI coordinates.

### Inputs

Ten laboratory XPT files, from `LAB18.XPT` in 1999–2000 through `BIOPRO_J.XPT` in 2017–2018. Exact paths and hashes are recorded in `provenance/standard_biochemistry_raw_manifest.csv`.

### Raw import

The minimally interpreted import contains 69,132 unique participants across ten cycles. It retains five historical 2001–2002 aliases alongside canonical variable slots so renaming remains visible at the next boundary.

### Variable-name harmonisation

Five 2001–2002 variables use an `LBD` rather than `LBX` prefix and are copied into canonical fields:

- `LBDSAPSI → LBXSAPSI`;
- `LBDSLDSI → LBXSLDSI`;
- `LBDSCR → LBXSCR`;
- `LBDSPH → LBXSPH`;
- `LBDSTB → LBXSTB`.

### Assay calibration

Published frozen crossover equations are applied only to affected cycles. Creatinine is corrected in 1999–2000 and 2005–2006. Twelve analytes are mapped from the 2017–2018 Cobas 6000 scale to the historical DxC 660i scale.

Examples:

\[
\mathrm{creatinine}_{1999-2000}^{*}=1.013x+0.147,
\]

\[
\mathrm{creatinine}_{2005-2006}^{*}=0.978x-0.016,
\]

and alkaline phosphatase uses the log-linear crossover:

\[
\mathrm{ALP}^{*}=10^{1.001\log_{10}(x)-0.04294}.
\]

All equations are stored in `config/standard_biochemistry_harmonisation.csv`; counts and before/after summaries are stored in `outputs/02_harmonised/standard_biochemistry_change_log.csv`.

### Plausibility and complete-case rules

All 23 retained analytes are strictly positive, so values at or below zero are converted to missing after assay calibration. No upper clinical tail is clipped: severe hyperglycaemia, azotaemia, transaminitis and other genuine pathology remain available.

The observed guard removed 27 total-bilirubin values recorded as zero and no non-positive values for the other 22 analytes. A historical code comment anticipated that the 2017–2018 iron correction could produce negative values, but none occurred in these frozen inputs. The audit reports the observed counts rather than treating the comment as a result.

No additional biochemistry error bounds are active. Participants are then required to have all 23 fields for this component-level output.

### Code

- `src/01_import/import_standard_biochemistry.py`
- `src/02_harmonise/harmonise_standard_biochemistry.py`
- `src/10_validate/compare_standard_biochemistry.py`

### Outputs

- `outputs/01_imported/standard_biochemistry_raw.csv`
- `outputs/02_harmonised/standard_biochemistry.csv`
- cycle, missingness, calibration and non-positive-value audit tables beside the main output.

### Reconciliation result

| Check | Result |
|---|---|
| Raw imported participants | 69,132 |
| Complete harmonised participants | 63,570 |
| Analytes | 23 |
| Cycles | 10 |
| `SEQN` membership | Exact match |
| Column names and order | Exact match |
| Maximum marker difference | 0 |
| Raw XPT hashes matching certified manifest | 10 of 10 |
| Non-positive values converted to missing | 27, all total bilirubin |

The rebuilt table is numerically identical to the certified release-candidate output for every participant and analyte.

### Interpretation boundary

The component table includes analytes that are later excluded, used only as derivation inputs, or superseded by dedicated components. Clean availability does not imply membership in the final 38-marker panel. `config/biomarker_registry.csv` records these roles explicitly.

### Status

Verified.

## Stage 7 — Anthropometry

### Scientific structure

The measured layer contains weight (`BMXWT`), standing height (`BMXHT`) and
waist circumference (`BMXWAIST`). BMI and waist-to-height ratio are calculated
in a separate derived layer because neither is an independent measurement.

The two derived variables have different roles. BMI is required to define the
young healthy reference but is deliberately excluded from biomarker scoring.
WHtR is the scored anthropometric marker. Scoring BMI after selecting the
reference on BMI would compress its reference variance and create circularity;
including both BMI and WHtR would also double-count closely related adiposity.

### Import and retention

Ten examination files contribute 96,766 participants and all ten hashes match
the certified manifest. No raw physiological bounds are applied because the
component includes children before adult eligibility is merged and extreme
adult body sizes can be genuine.

Only 1,128 examinations missing all three measurements are removed. Partial
rows remain because BMI requires weight plus height, whereas WHtR requires waist
plus height. Requiring all three prematurely would discard valid values for one
derived feature merely because the other cannot be calculated.

### Equations and naming

\[
\mathrm{BMI}=\frac{\mathrm{weight\ (kg)}}{[\mathrm{height\ (cm)}/100]^2},
\qquad
\mathrm{WHtR}=\frac{\mathrm{waist\ (cm)}}{\mathrm{height\ (cm)}}.
\]

The second variable is explicitly named WHtR. The historic abbreviation `WHR`
could imply waist-to-hip ratio, but hip circumference is not available across
all ten cycles. There are 87,799 usable BMI values and 84,473 usable WHtR values.

### Reconciliation

| Check | Result |
|---|---|
| Raw examination participants | 96,766 |
| Participants with at least one measurement | 95,638 |
| BMI values | 87,799 |
| WHtR values | 84,473 |
| Raw source hashes | 10 of 10 certified matches |
| Measured output | Exact match |
| Derived output | Exact match |
| Maximum numerical difference | 0 |
| Repository tests | 26 passed |

Rules are recorded in `config/anthropometry_harmonisation.csv`, and exact
comparison results are in `provenance/anthropometry_comparison.csv`.

### Status

Verified.

## Stage 8 — Lipids

### Scientific structure

The lipid reconstruction distinguishes two sampling frames. Total cholesterol
(`LBXTC`) and HDL cholesterol are full-sample measurements and enter the final
panel. Dedicated triglycerides (`LBXTR`) come from the morning fasting
subsample. They are reconstructed for provenance but are superseded in the
final panel by full-sample BIOPRO triglycerides (`LBXSTR`) from the already
verified standard-biochemistry block.

### Cross-cycle harmonisation

In 1999–2004, total cholesterol and HDL share one source file. From 2005 onward
they occupy separate files and are joined one-to-one by `SEQN`. HDL also changes
name from `LBDHDL` in 1999–2002, to `LBXHDD` in 2003–2004, to `LBDHDD` from
2005 onward. All are harmonised to canonical `LBDHDD`.

No clinical lipid bounds are applied during cleaning. Very high or low values
may represent real disease, and clinical reference ranges belong downstream.
Complete-casing requires both total cholesterol and HDL within this component.

### Triglyceride ownership and LDL exclusion

The dedicated fasting sources yield only 32,409 complete triglyceride values,
compared with 74,681 complete full-sample cholesterol/HDL records. Admitting a
fasting-only marker would substantially contract the shared locked cohort and
would require fasting survey weights. The final triglyceride is therefore
`LBXSTR`, the 63,570-participant full-sample BIOPRO field already calibrated and
verified in Stage 6.

Calculated LDL (`LBDLDL`) is excluded. It is Friedewald-derived from total
cholesterol, HDL and triglycerides, is therefore deterministic rather than an
independent measurement, and is systematically unavailable at high
triglycerides.

### Reconciliation

| Check | Result |
|---|---|
| Unique raw files | 27 |
| Raw hashes matching certified manifest | 27 of 27 |
| Imported TC/HDL participant pairs | 82,539 |
| Complete total cholesterol and HDL | 74,681 |
| Imported fasting records | 35,456 |
| Complete fasting triglycerides | 32,409 |
| Full-sample lipid output | Exact match |
| Fasting triglyceride output | Exact match |
| Maximum numerical difference | 0 |
| Repository tests | 29 passed |

The ownership decision and variable-name mapping are explicit in
`config/lipid_harmonisation.csv`; comparisons are in
`provenance/lipid_comparison.csv`.

### Status

Verified.

## Stage 9 — Glycated haemoglobin and RBC folate

### Glycated haemoglobin

HbA1c (`LBXGH`) is a single measured metabolic marker available in all ten
cycles. The import contains 69,132 laboratory participants. Values outside the
inclusive 3–20% impossibility bounds are made missing; this affected two
observations. Complete-casing then yields 64,950 measurements.

No cross-cycle calibration equation is applied. NHANES withdrew the proposed
2005–2006 crossover regression and retained the investigated 2007–2010 shift
unchanged. The frozen methodology therefore documents that shift and preserves
`cycle` for later adjustment rather than modifying participant values.

### RBC-folate assay eras

RBC folate requires the most involved assay harmonisation in this foundation:

- 1999–2006 uses Bio-Rad `LBDRBFSI`, serum folate `LBDFOLSI`, and analyser-
  reported `LBXHCT` in the Fazili/NHANES conversion;
- 2007–2010 uses microbiologic-assay `LBXRBFSI` directly;
- 2011–2018 uses microbiologic-assay `LBDRFOSI` directly.

The early-cycle haematocrit is read from the original CBC source files because
that measured field is the specified crossover input. It is not substituted by
the separately reconstructed calculated haematocrit.

### Four-step Bio-Rad conversion

For haematocrit fraction `h=HCT/100`, the frozen conversion calculates:

\[
WBF=RBF\,h+FOL(1-h),
\]

\[
WBF_{adj}=10^{0.2204+1.017\log_{10}(WBF)},
\]

\[
FOL_{adj}=10^{0.0188x^3-2.7109x^{-0.5}+3.8276},
\quad x=\log_{10}(FOL),
\]

and

\[
RBC\ folate=\frac{WBF_{adj}-FOL_{adj}(1-h)}{h}.
\]

Infinite, undefined and non-positive results are invalid. No finite
non-positive result occurred in the frozen inputs. After missing dependencies
and invalid logarithmic edge cases are removed, 78,147 harmonised RBC-folate
values remain. The previously omitted 2003–2004 cycle is explicitly present.

### Reconciliation

| Check | Result |
|---|---|
| Unique raw files including HCT dependencies | 24 |
| Raw hashes matching certified manifest | 24 of 24 |
| Complete HbA1c values | 64,950 |
| Harmonised RBC-folate values | 78,147 |
| HbA1c output | Exact match |
| RBC-folate output | Exact match |
| Maximum numerical difference | 0 |
| Repository tests | 32 passed |

Rules are in `config/hba1c_folate_harmonisation.csv`, raw crossover inputs are
preserved in `outputs/01_imported/folate_raw_inputs.csv`, and exact comparisons
are in `provenance/hba1c_folate_comparison.csv`.

### Status

Verified.

## Stage 10 — Locked untransformed final-38 panel

### Purpose and definition of raw

This stage creates the single authoritative input for all later transformation
and standardisation work. “Raw panel” means values remain in interpretable
biomarker units: no skew transformation, z-standardisation or winsorisation has
occurred. The values have already passed the documented source import, assay
harmonisation, physiological screening and mathematical derivation stages.

### Merge boundary

The merge starts with 34,785 eligible adults and left-joins each verified
component using unique `SEQN`. Participant count therefore remains 34,785 while
component availability accumulates. The cumulative complete counts are:

| Boundary | Complete participants |
|---|---:|
| Eligible base | 34,785 |
| Standard biochemistry | 32,707 |
| Cardiovascular | 31,438 |
| CBC derived | 31,309 |
| HbA1c | 31,265 |
| Anthropometry | 30,588 |
| Total cholesterol and HDL | 30,572 |
| Measured CBC | 30,572 |
| RBC folate | 29,211 |
| ACR and eGFR | 29,053 |

The marker-level attrition table is sequential: `newly_excluded` means first
encountered at that point in the frozen completeness order. It does not imply
that a participant has only one missing marker.

### Historical lock versus public panel

Exact cohort reconstruction requires complete data on the final 38 plus six
historical auxiliary fields. BMI and calibrated serum creatinine are screening
fields. Osmolality, calculated HCT, calculated MCH and measured RBC count were
present when the sample was locked but later pruned from scoring as redundant.

Those six fields are isolated in
`final_38_panel_screening_metadata.csv`. They do not appear in the authoritative
panel, preventing them from being mistaken for final biomarkers.

### Certified output

`final_38_biomarker_panel_raw.csv` contains:

- one unique row per participant;
- `SEQN`, `RIAGENDR` and `RIDAGEYR`;
- exactly 38 biomarkers in frozen specification order;
- 29,053 participants: 14,916 men and 14,137 women;
- no missing values;
- no transformed or standardised coordinates.

Every participant and value exactly matches the corresponding subset of the
certified harmonised table. Maximum numerical difference is zero. The file's
SHA-256 is:

`D9B2E795379571BFF781F2EE6503FA81B6AF5C04E59C47103519A0C25DC7F07A`

### Audit outputs

- `outputs/05_panel/final_38_biomarker_panel_raw.csv`
- `outputs/05_panel/final_38_panel_screening_metadata.csv`
- `outputs/05_panel/final_38_panel_merge_flow.csv`
- `outputs/05_panel/final_38_panel_missingness.csv`
- `outputs/05_panel/final_38_panel_attrition_by_marker.csv`
- `outputs/05_panel/final_38_panel_source_map.csv`
- `provenance/final_38_panel_comparison.csv`
- `provenance/final_38_panel_certification.csv`

### Status

Locked and verified.

## Stage 11 — full-cohort transformation selection

The five manuscript-era candidate formulas were reassessed prospectively on
the entire locked 29,053-person panel. Thirteen markers pass the identity gate.
Of 25 entering nonlinear testing, 13 have a fully adequate least-aggressive
choice, four have a borderline best available choice and eight remain
inadequate under all permitted formulas. No reference subset, outcome or
downstream model result entered selection.

All 38 choices were applied without offsets, missing values or infinities. The
new recipe agrees with 20 of 38 frozen manuscript transformations. Full-cohort
HD stress testing finds a full-rank, well-conditioned correlation structure but
materially heavy-tailed squared distances, so strict distributional HD
readiness is not yet supported.

The full rationale, selected recipe and stability results are in
`TRANSFORMATION_PROTOCOL.md` and `TRANSFORMATION_RESULTS.md`. The transformed,
unstandardised output SHA-256 is
`E2F404E78414F9C33B9DCB8DA5BC905C3AD25F4C47472E406F13564CB9243EDC`.

### Status

Two-step selection locked and full-cohort HD stress test complete;
reference-specific HD assessment and standardisation not begun.

## Restart checkpoint

The inclusion-preserving shared master is now complete at 29,053 participants x
89 columns (`outputs/07_master/analysis_master_raw.csv`; SHA-256
`F84DF97D19A6FAE316D5482022F97AC9CBB789A7487353353C9DD1A58E99C692`). It
contains the unchanged raw 38-marker panel plus eligibility provenance,
reference-screening covariates and mortality. Optional missing covariates have
not caused participant deletion.

The availability-aware young-healthy reference is now locked at 1,238 (547 men,
691 women). **Next task:** branch into HD sex-specific mean/SD/covariance
standardisation and sensitised-FI sex-specific median/robust-scale grading. See
`HEALTHY_REFERENCE.md` for the reference audit.

## HD-standardisation checkpoint

The transformed biomarker panel has now been standardised separately by sex
using the locked healthy-reference mean and sample SD for each marker. The
result contains all 29,053 participants, the reference flag and 38 finite
coordinates. Reference means are zero and SDs one to numerical precision.

The separate output is
`outputs/09_standardised/final_38_panel_transformed_standardised_hd.csv`, SHA-256
`D21583D09C00B49097AEF53BE7F43287CD1C233281D34FB22B7EDB8305011B78`.
Covariance/HD scoring and FI robust grading have not yet been performed.
