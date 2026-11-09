# Reproduction pipeline

This folder makes the paper's Results tables, figures and quoted numbers from
the prepared NHANES data, and checks each number against the reviewed analysis
and the manuscript as printed.

It is the short version of the analysis. The full record of how every number
was first produced and reviewed (every sensitivity analysis, audit and run log)
is the authors' working archive, which is not part of this code.

**Quick start with the prepared tables.** Instead of running the two
preparation stages, you can download the prepared tables from the data
deposit (see the README at the top of this code), unzip them into this code's
top folder, and run the Results straight away.

## Run it

From this folder:

```
python -B make_all.py --dry-run     # list the commands; runs nothing
python -B make_all.py               # SHORT run (the default): everything except the long resampling
python -B make_all.py --full        # FULL run: every bootstrap and all 472 panel evaluations (hours)
python -B check_against_paper.py --all   # the check on its own (seconds; reads saved outputs only)
```

**The short run** computes every point estimate, table and figure, plus a few
bootstrap draws (and, for Results 4, five panel evaluations) in each section, and
compares those draws one by one with the archived draws (stored in
`expected/`). It took about 20 minutes on a 12-core desktop (6 Nov 2026). It
ends with `check_against_paper.py --all`; a pass looks like

```
s1: 38 numbers, 1 table(s): PASS [20 intervals pending a full run; 4 draws or evaluations match the archive]
...
s4: 45 numbers, 0 table(s): PASS [... 62 draws or evaluations match the archive]
```

followed by one allowed formatting difference (Table 1, "17.0" for the
archive's "17") and the known slips in the manuscript (below). Exit code 0
means everything passed. Every point estimate is checked; the intervals are
reported as pending, because they need every draw.

**The full run** recomputes every interval. Approximate times on 12 cores:
Results 1 a few minutes; 2 about an hour; 3 about 1.5 hours; 4 about 30
minutes plus the benchmark bootstrap; 5, 6 and 7 under an hour each (Results
1 and 2 timed on 6 Nov 2026; the others estimated from the archived runs).

Each section script also runs on its own; its header lists its options
(`python -B s4_panel.py --help`). Unknown options stop a script before it
does any work.

**Requirements.** Python 3.10 with numpy, pandas, scipy, statsmodels,
lifelines, patsy, matplotlib, openpyxl and scikit-learn. Tested with Python
3.10.11, numpy 1.26.4, pandas 2.2.3, scipy 1.15.3, statsmodels 0.14.5,
lifelines 0.30.0, patsy 1.0.2, matplotlib 3.10.8, openpyxl 3.1.5,
scikit-learn 1.7.2. Each run records the versions it used in its `run.json`.

**Before the first run**, run `../1_prepare_nhanes/` and then
`../2_prepare_nhanes_iii/` (see their READMEs).

## Known slips in the manuscript

The checker lists these as known errors, not failures; the pipeline agrees with
the reviewed analysis in each case.

| Where | Printed | Reviewed value | Rounds to |
|---|---|---|---|
| Results 2, paragraph [39], SS adjusted female HR, lower limit | 0.82 | 0.8149 | 0.81 |
| Results 4, paragraph [68], female HR with the full-panel score | 0.93 | 0.9248 | 0.92 |
| Results 4, paragraph [68], female HR with the locally refitted inputs | 0.66 | 0.6550 | 0.65 |

## What each file does

| File | Role |
|---|---|
| `settings.py` | Folders, seeds, age limit, draw counts, and the write guard. |
| `runs.py` | Each section's run record (`outputs/<section>/run.json`) and the rules for reusing stored bootstrap draws. |
| `data.py` | Loads the two analysis cohorts (NHANES 1999–2018, NHANES III): participants, endpoint, follow-up and marker coordinates. |
| `measures.py` | The distances from the young-healthy reference (DA, DM, DE) and the three reference choices (SS, SC, PC). |
| `models.py` | Cox regression (ordinary and ridge-penalised), the women-versus-men contrast, the cross-fitted score q, the bootstraps and percentile intervals. Records every fit warning. |
| `prediction.py` | Held-out prediction: concordance, five-year risk and its error (Table 2, Phenotypic Age comparison). |
| `displays.py` | Shared figure style; writes tables as CSV and as all-text Excel; loads the shared Figure 1/2 renderer. |
| `s1_sex_gap.py` … `s7_functional.py` | One script per Results section. Each header states its question, population, endpoint, estimands, split, uncertainty, assumptions and permitted interpretation. |
| `make_all.py` | Runs every section in order (short or full), then the check. |
| `check_against_paper.py` | Compares every number with the reviewed value and the printed text; tables as exact text; draws one by one. |
| `inputs/` | Fixed analysis choices too small to be data: clinical and young-healthy limits and the 20 cycle divisions (Results 3); the elastic-net penalty chosen for each panel selection (Results 4, below); the 32 markers common to both eras with their transformations (Results 5, 6); the reviewed marker ranking that fixes Results 7's ten markers. |
| `render_figures.py` | Draws the accepted Figures 1 and 2 from the results the pipeline hands it (copied unchanged from the archive's `displays/current_figures/`). |
| `expected/` | The reviewed values and the printed values; the first archived draws of each section; per-column draw tolerances and allowed table-format differences, each with its reason; `interval_verification.csv` (below). |
| `outputs/` | Everything the pipeline writes, one folder per section. |

## The one step that reads stored choices

Choosing the elastic-net penalty for each of Results 4's 472 panel selections
took about 177 CPU-hours (an inner grid search). `inputs/panel_selection_tuning.csv`
holds the choices made in the reviewed run; `s4_panel.py` repeats every step
after them (ranking, refitting, evaluation), which reproduces the reviewed
selections. The search itself is in the archive and is not part of the pipeline.

## Where the inputs come from

The inputs are the outputs of the two preparation stages, which run first:
`../1_prepare_nhanes/outputs/` (with its `config/`) and
`../2_prepare_nhanes_iii/outputs/`. They turn the raw NHANES files into one
row per participant: markers transformed towards a normal shape, the
young-healthy reference chosen, and each marker expressed relative to that
reference. The prepared CRP, fasting glucose and spirometry extracts used by
Results 2, 4 and 7 (`settings.EXTRACTS`) are made by
`1_prepare_nhanes/6_extra_markers/`. The NHANES III questionnaire and
examination variables and the linked mortality file are copied into prepared
tables by `2_prepare_nhanes_iii/05_extract_results_variables.py`; the Results
read no raw file.

One input is a stored result: `inputs/full_development_ranking.csv`, the
marker ranking of the reviewed Results 4 run, which fixes Results 7's ten
markers (`settings.LOCKED_TEN`). A full run of `s4_panel.py` remakes it as
`outputs/s4/locked_ten_modern.csv`.

The prepared files are not in git because of their size; the preparation
stages make them. Every file a run read, with its SHA-256, is listed in that
run's `run.json`.

## Safety: the write guard, run records and stored draws

**Write guard.** Importing `settings.py` switches on a guard. Inside the input
repository and inside the repository that holds this pipeline, a script may
write only in `outputs/` (and the system's temporary folder); any other write
stops it. Writes outside both repositories (e.g. matplotlib's font cache) are
allowed. `NHANES_ALLOW_WRITES=1` (exactly "1") switches it off; a normal run
never needs that.

**Run records.** Each section writes `outputs/<section>/run.json`: the command,
the check draws it wrote, the SHA-256 of every input file it or any of its
worker processes read (each process lists its reads in
`outputs/<section>/process_logs/`), a
fingerprint of the code (comments and docstrings excluded), the package
versions, and every model fit that reported a warning or failed to converge.
A fit with non-finite coefficients stops the run.

**Stored draws.** A full bootstrap's draws are stored with a manifest
(`<file>.manifest.json`). Reuse always needs the file unchanged since its
manifest (SHA-256), the same seed, and the same data (no listed input changed,
no new input read). Then:

* draws made by a full run of the pipeline (`"provenance": "computed"`) are
  reused only if the code is unchanged. After a code change the script stops
  and the draws must be remade by a full run; a few matching check draws are
  never taken as evidence for the rest.
* draws registered as an archived bootstrap (`"provenance":
  "archive_registered"`) are reused and reported as the archived draws, not as
  draws made by the current code. The authors' registration tool (not part of
  this code) registers a file only if every draw and every archived quantity
  equals the complete archived bootstrap. A fresh copy of this code has no
  stored draws: a full run makes them.

Check draws made in the same run must agree with the stored draws; otherwise
the script stops and says how to remake them.

## How the check works

`expected/expected_values.csv` holds one row per number quoted in the Results:
the reviewed value at full precision, its source, and the value as printed. A
section passes when every number equals the reviewed value (to 1e-8, scaled,
unless a looser tolerance is recorded with its reason) and reads exactly as
printed, every table matches the reviewed table as text, and every check draw
matches the archived draw.

Intervals need every draw, and some draws take hours. Two checks stand in for a
full run:

1. **Draw by draw.** The short run's draws must equal the archived draws one by
   one (same random numbers, same estimates). Every archived quantity must be
   present; quantities without archived draws are listed as unverified.
2. **The interval rule on the archived draws** (done by the authors; the
   archive is too large to publish). For each section, a verifier applied the
   pipeline's own interval code to the complete archived draws and compared the
   result with every reviewed and printed interval. All 163 intervals reproduce,
   and for Results 4 so do the 31 selection summaries (from the 5,664 archived
   evaluations) and Figure 4's two source tables. A verifier writes its receipt
   only when everything passes (otherwise it exits with status 1):
   `expected/interval_verification.csv` holds each verified bound or summary
   value and the number of draws or evaluations used;
   `expected/interval_verification_sources.json` holds, per section, the code
   fingerprints of the pipeline modules and verifier scripts used and the
   SHA-256 of every archived file read. The checker counts an interval as
   "verified on the archived results" only while those fingerprints still hold
   and the receipt's numbers still equal the reviewed and printed values with
   a full run's draw count. Without the archive (as in a fresh copy of this
   code) the receipt cannot be confirmed, and its intervals are reported as
   pending.

Neither is the same as a full run, which is how anyone else recomputes the
intervals (`--full`).

## Figures

* Figures 1 and 2 are drawn by the accepted shared renderer
  (`render_figures.py`) from the pipeline's computed results.
* Figure 4: the pipeline writes its two source tables
  (`outputs/s4/figure4_*_source_data.csv`) in the layout the accepted renderer
  reads; they equal the accepted tables. The drawing itself is not yet wired in.
* Figures 3, 5, 6 and 7 are drawn by the section scripts as checking displays,
  with interval bars once the bootstrap exists. The published versions were
  drawn from the same results by renderers in the authors' archive.

## Verified by the authors (6 November 2026)

Run from the prepared tables made by this code's two preparation stages:

| Section | Short check | Intervals |
|---|---|---|
| Results 1 | Passes; 2 draws per era match the archive | Full run (500 draws per era): every interval equals the reviewed and printed value |
| Results 2 | Passes; 2 draws of each bootstrap match | Full run (both bootstraps): every interval, Table 2 included, equals the reviewed and printed value |
| Results 3 | Passes; 2 draws match (all 16 quantities) | 16 verified on the archived draws |
| Results 4 | Passes; 5 evaluations (incl. NHANES III fold F09) and 2 benchmark draws match | 8 intervals and 31 summaries verified on the archive |
| Results 5 | Passes; 3 draws match | 16 verified on the archived draws |
| Results 6 | Passes; 2 draws per era match (261 quantities) | 55 verified on the archived draws |
| Results 7 | Passes; 1 draw per era matches (both gains) | 68 verified on the archived draws |
