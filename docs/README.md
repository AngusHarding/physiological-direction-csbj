# Method notes

Written while the analysis was built, copied unchanged. Some refer to old
folder names and to checks against the original analysis; the method they
describe is the one the code applies.

| Note | Covers |
|---|---|
| `endpoint.md` | Accident-censored mortality: which deaths count, follow-up, censoring |
| `data_preparation_overview.md` | The preparation stages and the rules they share |
| `mortality_linkage.md` | Reading the 2019 linked mortality files |
| `health_conditions.md` | Rebuilding the self-reported condition flags |
| `analysis_master.md` | The 29,053-person analysis master |
| `young_healthy_reference.md` | Who is in the young-healthy reference, and how missing answers are handled |
| `transformations.md` | How each marker's transformation was chosen |
| `standardisation.md` | Standardising markers against the sex-specific reference (DM coordinates) |
| `adverse_distance.md` | The adverse-distance coordinates (called "FI" in the note) |
| `nhanes_iii.md` | The NHANES III build |

## Old script names in these notes

These notes were written for the original analysis, whose scripts lived in
`src/`. Where a note names one, this table gives its file in this code (or says
that it was a check of the original analysis and is not included).

| Named in the notes | In this code |
|---|---|
| `src/01_import/import_alcohol.py` | `1_prepare_nhanes/4_lifestyle/03_import_alcohol.py` |
| `src/01_import/import_cardiovascular.py` | `1_prepare_nhanes/3_body_and_bp/04_import_blood_pressure.py` |
| `src/01_import/import_demographics.py` | `1_prepare_nhanes/1_cohort/01_import_demographics.py` |
| `src/01_import/import_physical_activity.py` | `1_prepare_nhanes/4_lifestyle/05_import_physical_activity.py` |
| `src/01_import/import_recent_illness.py` | `1_prepare_nhanes/1_cohort/02_import_recent_illness.py` |
| `src/01_import/import_smoking.py` | `1_prepare_nhanes/4_lifestyle/01_import_smoking.py` |
| `src/01_import/import_standard_biochemistry.py` | `1_prepare_nhanes/2_blood_tests/01_import_standard_biochemistry.py` |
| `src/01_import/nhanes_reader.py` | `1_prepare_nhanes/shared/nhanes_reader.py` |
| `src/02_harmonise/harmonise_alcohol.py` | `1_prepare_nhanes/4_lifestyle/04_harmonise_alcohol.py` |
| `src/02_harmonise/harmonise_cardiovascular.py` | `1_prepare_nhanes/3_body_and_bp/05_harmonise_blood_pressure.py` |
| `src/02_harmonise/harmonise_physical_activity.py` | `1_prepare_nhanes/4_lifestyle/06_harmonise_physical_activity.py` |
| `src/02_harmonise/harmonise_smoking.py` | `1_prepare_nhanes/4_lifestyle/02_harmonise_smoking.py` |
| `src/02_harmonise/harmonise_standard_biochemistry.py` | `1_prepare_nhanes/2_blood_tests/02_harmonise_standard_biochemistry.py` |
| `src/04_eligibility/build_eligible_cohort.py` | `1_prepare_nhanes/1_cohort/03_build_eligible_cohort.py` |
| `src/10_validate/audit_disease_candidates.py` | not included (original analysis only) |
| `src/10_validate/compare_alcohol.py` | not included (original analysis only) |
| `src/10_validate/compare_cardiovascular.py` | not included (original analysis only) |
| `src/10_validate/compare_physical_activity.py` | not included (original analysis only) |
| `src/10_validate/compare_self_reported_disease.py` | not included (original analysis only) |
| `src/10_validate/compare_smoking.py` | not included (original analysis only) |
| `src/10_validate/compare_standard_biochemistry.py` | not included (original analysis only) |
| `src/10_validate/validate_fi_scoring_spec.py` | not included (original analysis only) |
| `src/11_fi/audit_fi_quality.py` | not included (original analysis only) |
| `src/11_fi/build_fi_adverse_excursions.py` | `1_prepare_nhanes/8_adverse_distance/04_build_adverse_excursions.py` |
| `src/11_fi/build_fi_deficit_severities.py` | `1_prepare_nhanes/8_adverse_distance/05_build_deficit_severities.py` |
| `src/11_fi/build_fi_robust_coordinates.py` | `1_prepare_nhanes/8_adverse_distance/03_build_robust_coordinates.py` |
| `src/11_fi/build_fi_signed_deficits.py` | `1_prepare_nhanes/8_adverse_distance/06_build_signed_deficits.py` |
| `src/11_fi/calculate_fi_scores.py` | not included (original analysis only) |
| `src/11_fi/estimate_fi_reference_parameters.py` | `1_prepare_nhanes/8_adverse_distance/02_estimate_reference_parameters.py` |
| `src/11_fi/lock_fi_reference.py` | `1_prepare_nhanes/8_adverse_distance/01_lock_reference.py` |
