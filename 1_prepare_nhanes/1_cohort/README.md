# 1. Cohort: who is in, survey design, deaths

| Step | Script | Makes |
|---|---|---|
| 1 | `01_import_demographics.py` | Age, sex, ethnicity, pregnancy status per participant (`demographics.csv`) |
| 2 | `02_import_recent_illness.py` | Recent acute illness questions, used to exclude people unwell at the exam (`recent_illness.csv`) |
| 3 | `03_build_eligible_cohort.py` | Eligibility flags: adults 18–79, not acutely ill, not pregnant. Records why each excluded person was excluded (`eligibility_flags.csv`, `eligible_adult_cohort.csv`) |
| 4 | `04_import_survey_design.py` | Survey weights, strata and clusters (`survey_design_1999_2018.csv`) |
| 5 | `05_import_socioeconomic.py` | Income-to-poverty ratio and education (`socioeconomic.csv`) |
| 6 | `06_import_mortality.py` | The linked mortality record: vital status, cause of death and both follow-up clocks (`mortality.csv`) |

Step 3 needs steps 1 and 2. The others are independent.

The mortality import keeps every death, including accidents. The paper's
endpoint, accident-censored mortality, is applied later in `3_results/data.py`
(see `docs/endpoint.md`).
