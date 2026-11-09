# 4. Lifestyle

| Step | Script | Makes |
|---|---|---|
| 1–2 | `01_import_smoking.py`, `02_harmonise_smoking.py` | Never, former or current smoker (`smoking_status.csv`) |
| 3–4 | `03_import_alcohol.py`, `04_harmonise_alcohol.py` | Average drinks per week; the questionnaire changed in 2017–2018 and both versions are converted to the same scale (`alcohol_drinks_per_week.csv`) |
| 5–6 | `05_import_physical_activity.py`, `06_harmonise_physical_activity.py` | Weekly activity in MET-minutes and an active/inactive flag (`physical_activity.csv`, `physical_activity_inclusion_preserving.csv`) |
| 7–8 | `07_import_hei_package_data.py`, `08_calculate_hei_2020.py` | Healthy Eating Index 2020 score and its 13 components, 2005–2018 only (`hei_2020_total.csv`, `hei_2020_component_scores.csv`) |

Smoking and alcohol screen the young-healthy reference (never smoked, at most
seven drinks a week). Diet never removes anyone from the cohort.

The HEI scripts need the R package data from `heiscore` and `heiscore.data`,
read in Python with `pyreadr`.
