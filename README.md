# Distance, Survival and Sex: code

Code for Harding et al. (2026), *Distance, Survival and Sex: Addressing a
Three-Way Trade-Off in Healthspan Trial Recruitment*, Computational and
Structural Biotechnology Journal.

It goes from the raw NHANES files to every number quoted in the paper's Results,
except five robustness checks (see "Not included"), and to Tables 1, 2 and 4,
and ends by checking each number against the paper. It also draws Figures 1 and
2 as published, the source tables of Figure 4, and checking versions of Figures
3 and 5–7 (see `3_results/README.md`).

```
raw NHANES 1999–2018 ──► 1_prepare_nhanes ──────┐
                                                ├──► 3_results ──► Tables 1, 2, 4, figures and
raw NHANES III ────────► 2_prepare_nhanes_iii ──┘                 the numbers quoted in the Results,
                                                                  then check_against_paper: PASS or FAIL
```

## Two ways to reproduce the paper

**Quick: start from the prepared tables.** Download `prepared_tables.zip` from
the data deposit (https://doi.org/10.5281/zenodo.24096988) and unzip it into this
folder; it fills `1_prepare_nhanes/outputs/` and `2_prepare_nhanes_iii/outputs/`.
Then:

```
cd 3_results
python -B make_all.py
```

The short run checks every point estimate and Tables 1, 2 and 4 in about 20
minutes; `python -B make_all.py --full` also recomputes every 95% interval
(several hours). See `3_results/README.md`.

**Full: start from the raw NHANES files.** Download the public files listed in
`data_locations.py`, write their folders into that file (the only file to
edit), then run the two preparation stages before the Results:

```
cd 1_prepare_nhanes
python -B run_all.py
cd ../2_prepare_nhanes_iii
python -B run_all.py
cd ../3_results
python -B make_all.py
```

The prepared tables you make this way are identical to the deposited ones;
`SHA256SUMS.txt` in the deposit lists their fingerprints.

## What is here

| Folder | What it does | Scripts | Lines |
|---|---|---|---|
| `1_prepare_nhanes/` | NHANES 1999–2018 only: raw files to one row per participant: eligibility, deaths, the 38 markers, lifestyle and conditions, the young-healthy reference, and the coordinates for each distance | 46, plus shared helpers | 4,100 |
| `2_prepare_nhanes_iii/` | NHANES III (1988–1994) only: the same, reusing the modern choices, plus the other NHANES III variables and the linked mortality file the Results use | 5, plus shared helpers | 1,400 |
| `3_results/` | Both eras: one script per Results section, shared helpers, and the check against the paper | 16 | 6,000 |
| `docs/` | Method notes: endpoint, reference, transformations, adverse distance | 10 notes | |
| `tests/` | Unit tests for the preparation scripts (`python -B -m unittest discover -s tests -p "test_*.py"`) | 28 | 1,200 |

Every analysis script was copied from the code that produced the paper. The
lines changed to make it run from this folder are marked `# public_code`.

Not included: the authors' full working archive, with the development analyses
that led to the final specification; five robustness checks quoted in the
Results, which come from that archive (Supplementary Tables S1.4–S1.6, S2.3,
S4.5 and S4.6: the all-cause and survey-design hazard ratios, the complete-panel
ratio, the alternative-reference comparison and the 32-marker score); the
Results 4 penalty search (the Results use its stored choices); the scripts the
authors used to verify the Results 3–7 intervals against the complete archived
bootstrap draws (a full run recomputes those intervals instead); and the
renderers that drew the published Figures 3–7. The paper's supplementary data
files hold the numbers behind every figure and table.

## Reading order

1. This page, then `1_prepare_nhanes/README.md`.
2. The preparation folders in number order. Each has a short README; read it
   before the scripts.
3. `2_prepare_nhanes_iii/README.md`. It reuses choices made for the modern data.
4. `3_results/README.md`, then `settings.py`, `data.py`, `measures.py`,
   `models.py`, `prediction.py`, then `s1_sex_gap.py` to `s7_functional.py`.

## Requirements

Python 3.10 with numpy, pandas, scipy, statsmodels, lifelines, patsy,
matplotlib, openpyxl, scikit-learn and pyreadr (HEI step only). Tested with Python 3.10.11,
numpy 1.26.4, pandas 2.2.3, scipy 1.15.3, statsmodels 0.14.5, lifelines 0.30.0,
patsy 1.0.2, matplotlib 3.10.8, openpyxl 3.1.5, scikit-learn 1.7.2 and pyreadr 0.5.6. Every
Results run records the versions it used.

## Names used in the code

Most code names follow the paper. Three older names remain inside the copied
scripts:

| In the code | In the paper |
|---|---|
| `fi`, "FI" | the coordinates for **adverse distance (DA)** |
| `hd`, "HD" | the standardised coordinates for **Mahalanobis distance (DM)** |
| `final_38` | the **38-marker panel** |

## Data

NHANES and the NCHS linked mortality files are public-use data from the US
National Center for Health Statistics (NCHS). The NCHS data user agreement
(https://www.cdc.gov/nchs/policy/data-user-agreement.html) applies to them and
to the prepared tables made from them: use the data for statistical reporting
and analysis only; make no attempt to learn the identity of any person;
do not link the data with individually identifiable data from other datasets;
and do not try to assess the disclosure methods that protect participants.
Please cite NCHS as the source of the data. The HEI-2020 scoring data come from
the CRAN package heiscore.data (CC0).

## How this code was written

The authors wrote, restructured and checked this code with two large language
models, Claude (Anthropic) and ChatGPT (OpenAI), working under author direction:
Claude mainly on the code and its comments, ChatGPT mainly on independent review
of the code against the Methods. Neither is an author; the authors are
responsible for the code. Commits made with Claude carry a `Co-Authored-By`
line. The corresponding author ran the code from the public NHANES files and
read it against the Methods before release (see the paper's statement on the
use of large language models).

## Licence and citation

The code is released under the MIT licence (`LICENSE`) and archived as release
v1.0.0 at https://doi.org/10.5281/zenodo.24097390. The prepared tables in the
data deposit (https://doi.org/10.5281/zenodo.24096988) are released under
CC BY 4.0. If you use either, please cite the paper (`CITATION.cff`).
