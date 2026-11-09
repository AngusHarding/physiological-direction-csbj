# Preparing NHANES 1999–2018

Turns the raw NHANES files into one row per participant. Folders run in number
order; each one's scripts run in number order too.

| Folder | What it makes |
|---|---|
| `1_cohort/` | Who is eligible (adults 18–79, not acutely ill, not pregnant), survey weights, income, and deaths from the linked mortality file |
| `2_blood_tests/` | Biochemistry, blood count, cholesterol, HbA1c, folate, and kidney measures (eGFR, urine albumin-to-creatinine ratio) |
| `3_body_and_bp/` | Body size (waist-to-height ratio) and blood pressure |
| `4_lifestyle/` | Smoking, alcohol, physical activity, diet quality (HEI-2020) |
| `5_health_conditions/` | Self-reported conditions and functional limits, used to screen the young-healthy reference |
| `6_extra_markers/` | CRP, fasting glucose and lung function (Results 2, 4 and 7) |
| `7_marker_panel/` | Joins the 38 markers, picks each marker's transformation, builds the analysis master and the young-healthy reference, and standardises every marker against that reference |
| `8_adverse_distance/` | The coordinates the adverse distance (DA) is built from |
| `config/` | Hand-written tables: which file and variable hold each measure in each survey cycle, and the panel and scoring specifications |
| `shared/` | `nhanes_reader.py`, which reads an NHANES `.XPT` file |

Folders 2 to 6 depend only on `1_cohort/` and can run in any order among
themselves. `7_marker_panel/` needs all of them; `8_adverse_distance/` needs
`7_marker_panel/`.

**Pattern inside each subject.** `import_` scripts copy the needed variables
out of the raw files, cycle by cycle, without changing them. `harmonise_`
scripts make the cycles comparable (unit changes, instrument changes, error
screens). `derive_` scripts compute measures from others, such as eGFR from
creatinine, age and sex.

What the Results read from here: 25 of the files these scripts make, plus the
scoring table `config/fi_scoring_spec.csv`.

## Running it

First set where the raw files are, in `../data_locations.py` (the one file to
edit). Then, from this folder:

```
python -B run_all.py
```

runs all 46 scripts in order and stops at the first that fails. Each script
also runs on its own, for example `python -B 1_cohort\01_import_demographics.py`.

- **Where things are:** `shared/locations.py` names every folder from
  `../data_locations.py`: the raw NHANES files, the mortality files, the HEI
  package data, and `outputs/`, where everything is written.
- **Safety:** importing `locations.py` switches on a write guard. A script can
  write only inside `outputs/`; anywhere else in this repository or in the data
  folders, it stops with "write guard: refused to write".
- **Slow step:** `7_marker_panel/02_select_transformations.py` repeats the
  transformation choice on 100 resamples to check its stability; expect
  minutes rather than seconds.

The scripts are the originals with their folder lines changed so they find
`locations.py`; every changed line is marked `# public_code:`.
