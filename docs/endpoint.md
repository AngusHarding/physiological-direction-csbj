# Mortality endpoint specification

This module holds the one definition of the mortality endpoint used by every
analysis in the repository. No other module derives a death event from the
cause-of-death code; a test enforces that.

## Source

NCHS 2019 public-use Linked Mortality Files for NHANES 1999-2018 and for
NHANES III, with follow-up to 31 December 2019. Participants are those eligible
for mortality linkage. Follow-up is measured from the examination.

## What the public-use files release

The underlying cause of death is released as the ten-category leading-cause
recode `UCOD_LEADING`:

| Code | Cause |
|---|---|
| 1 | Diseases of heart |
| 2 | Malignant neoplasms |
| 3 | Chronic lower respiratory diseases |
| 4 | Accidents (unintentional injuries) |
| 5 | Cerebrovascular diseases |
| 6 | Alzheimer's disease |
| 7 | Diabetes mellitus |
| 8 | Influenza and pneumonia |
| 9 | Nephritis, nephrotic syndrome and nephrosis |
| 10 | All other causes (residual) |

For NHANES 2015-2016 and 2017-2018 only codes 1, 2 and 10 are released. Deaths
from every other cause, including accidents, fall in the residual category and
cannot be separated from it. For NHANES 1999-2014 and for NHANES III all ten
codes are released. The NCHS codebook for the 2019 public-use files states the
restriction ("For 2015-2018 NHANES only, UCOD_LEADING=001, 002, and 010 are
available"), and the constructor checks it against the data on every call.

A few deaths have no cause code at all: 2 among eligible NHANES 1999-2018
participants (none in the analysis cohorts) and 13 among eligible NHANES III
participants.

The two multiple-cause flags, diabetes and hypertension recorded anywhere on the
death certificate, are released for every cycle.

## Endpoints

| Column | Role | Event | Defined for |
|---|---|---|---|
| `event_submitted_v1` | Principal | A death with a recorded cause, leaving out deaths identified as accidental (code 4) | Every participant. In 2015-2018 no death can be identified as accidental, so all coded deaths count |
| `event_allcause` | Check 1 | Any linked death, including a death with no cause code | Every participant |
| `event_nonaccidental` | Check 2 | A death with a recorded cause other than code 4 | NHANES 1999-2014 and NHANES III. Missing for 2015-2018 |

**Principal endpoint: accident-censored mortality.** Death from any recorded
cause, with accidental deaths censored at the time of death wherever the data
identify them. This is the
endpoint of the analysis as first reported; the column name records that. It
was retained as principal on 2026-10-08 by decision of the manuscript lead, to
remove accidental deaths to the extent the data allow and to keep the reported
analysis intact.

Its limitation is that the definition differs between survey cycles. In
1999-2014, identified accidental deaths are censored at death. In 2015-2018 they
cannot be identified and are counted as events. The survey cycles are also the
cross-validation folds, so the event is not defined identically across folds.
In the direction cohort this concerns the 92 deaths of 2015-2018, of which 41
are in the residual category that would contain any accident.

Two checks address the limitation:

- **Check 1, all-cause mortality.** One definition in every cycle. The main
  analyses, including panel selection, are repeated on it.
- **Check 2, non-accidental mortality in 1999-2014.** The principal definition
  applied only to the cycles where it is exact. The headline analyses are
  repeated on it. Panel selection is not, because its partition design needs
  all ten cycles.

In NHANES III the principal endpoint and check 2 coincide.

Under the principal endpoint and under check 2, deaths coded as accidents, and
deaths with no cause code, are censored at the time of death, so a Cox model
estimates a cause-specific hazard. Where a probability of death is wanted,
accidental death is a competing event.

A death with no cause code is censored under this endpoint because its cause is
unknown. It is not an identified non-accidental death, and it is not an
identified accident. This is a limitation of the sensitivity analysis.

Accident status is never inferred from the residual category.

## Cause-specific analyses

An analysis of death from one named cause selects deaths by that cause's code.
It is the only kind of analysis outside this module that reads the cause code,
and it must restrict its cohort to the cycles in which the cause is released as
its own category, using `cause_identifiable`:

| Cause | Cycles analysed |
|---|---|
| Diseases of heart (1), malignant neoplasms (2) | All ten, 1999-2018 |
| Every other code, 3 to 10 | 1999-2014 |
| Multiple-cause diabetes and hypertension flags | All ten, 1999-2018 |

The participants are restricted, not only the events. A 2015-2018 participant
who died of stroke has the death recorded under the residual category; left in a
stroke analysis, that death would count as censoring by another cause.

The residual category itself is restricted to 1999-2014. In 2015-2018 it holds
every cause except heart disease and cancer, so it is a different category from
the residual of earlier cycles.

Results for the restricted causes therefore describe 1999-2014 and are reported
with the number of participants and the cycles they rest on. NHANES III releases
all ten categories throughout and needs no restriction.

## Configuration

`code/endpoints.py` holds one constant, `ENDPOINT`, naming the endpoint that
every module analyses. Analysis code reads the column `endpoint_event`, which is
a copy of the configured endpoint. In this repository it is always the principal
endpoint.

A check is run in a separate copy of the analysis tree in which that one
constant is changed, so that check results cannot be written into the folders
that hold the principal results. To reproduce check 1: copy the repository, set
`ENDPOINT = "event_allcause"`, and rerun the modules.

To reproduce check 2: copy the repository, set
`ENDPOINT = "event_nonaccidental"`, and rerun the headline scripts of module 05
(`01_build_cohort.py`, `02_fit_loco_directions.py`, `04_sex_contrasts.py`,
`05_mortality_associations.py`, `09_bootstrap_uncertainty.py`). Nothing else has
to be changed: `endpoint_cycles` in `code/endpoints.py` says on which cycles an
endpoint is defined, and module 05 builds its cohort and its
leave-one-cycle-out folds from those cycles. Under the non-accidental endpoint
that is the eight cycles of 1999-2014: 9,362 of the 11,497 participants aged
52 to 79, with 2,402 events. Check 2 covers those five scripts and nothing
else. Module 05's split-stability script and panel selection (module 06) are
defined for ten cycles and stop with an error under this setting.

Check 2 differs from the principal analysis in more than the event definition.
It leaves out 2,135 participants and their 92 deaths, uses eight
cross-validation folds where the principal analysis uses ten, and leaves out
the two cycles with the shortest follow-up. Its result describes 1999-2014. A
difference from the ten-cycle result cannot be put down to the endpoint alone.
The marker coordinates are not rebuilt: they stay standardised to the same
young healthy reference as in the principal analysis.

Changing `ENDPOINT` is an outcome-relevant change. Modules with resumable
computations include the endpoint in their cache keys, so results computed under
one endpoint cannot be resumed under another.

## Tests

```
python -B -m unittest discover -s analysis/00_mortality_endpoint/tests
```

`test_endpoints.py` pins the truth table of each endpoint, the cycles in which
each cause is identifiable, and the event counts in the locked mortality import.
`test_no_local_endpoint_derivation.py` scans the analysis code and fails if a
module defines an event for itself, or selects deaths by cause code in
continuous NHANES without restricting its cohort to identifiable cycles.
