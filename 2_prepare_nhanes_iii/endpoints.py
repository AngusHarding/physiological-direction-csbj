# =============================================================================
# READ ME FIRST: 2_prepare_nhanes_iii/endpoints.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  States, in one place, the endpoint definitions and which cause
#               codes each cycle releases. The Results compute the principal
#               endpoint with the same rule in 3_results/data.py
#               (principal_endpoint); the release's run does not import this
#               file, which is kept as the documented definition (see
#               docs/endpoint.md).
#
# METHODS       "Endpoints": "The principal endpoint counted deaths with a
#               recorded underlying cause other than an identified accident.
#               Identified accidental deaths and deaths without a recorded cause
#               were censored at death. Accidents cannot be identified in the
#               public-use files for 2015-2018, where all deaths with a recorded
#               cause were counted".
#
# DEFINITIONS   event_submitted_v1 (PRINCIPAL) = died AND cause recorded AND
#                 cause != 4 (accident)
#               event_allcause (check) = died
#               event_nonaccidental (check) = the principal rule, 1999-2014
#                 and NHANES III only
# =============================================================================
"""The single definition of the mortality endpoint used by every analysis module.

Concern addressed
    Accidental deaths cannot be identified in the 2015-2018 public-use linked
    mortality files, so an endpoint described as "deaths excluding accidents in
    all cycles" was not what the submitted analysis computed. The cause
    classification behind the endpoint also needed to be stated explicitly.
    This module states it once, in code.

Audit objective (no statistical hypothesis applies)
    Make the event indicator mean the same thing in every survey cycle, and make
    it impossible for two modules to disagree about what a death event is.

Why a shared module, and what it replaces
    Before the revision the event was re-derived from `ucod_leading` in about
    seventeen separate places. Each copy was correct against the others, but a
    definition that lives in seventeen places cannot be changed safely or shown
    to be consistent. Every module now calls `attach_endpoints` instead, and a
    test (tests/test_no_local_endpoint_derivation.py) fails if any module
    derives an event from the cause code itself.

Population, endpoint, time origin
    Unchanged by this module: participants eligible for mortality linkage in the
    NCHS 2019 public-use linked mortality files, followed from the examination
    (`permth_exm`) to death or 31 December 2019. This module defines only WHO
    COUNTS AS AN EVENT, never who is in the cohort or how time is measured.

The three endpoints, and the role of each
    event_submitted_v1    PRINCIPAL. Death with a recorded underlying cause,
                          leaving out deaths identified as accidental. Accidents
                          are removed wherever the public-use files identify
                          them: NHANES 1999-2014 and NHANES III. In NHANES
                          2015-2018 they cannot be identified, so there every
                          death with a recorded cause counts. This is the
                          endpoint of the analysis as first reported, which is
                          where its name comes from. It stays the principal
                          endpoint by decision of the manuscript lead
                          (2026-10-08): it removes accidental deaths to the
                          extent the data allow and keeps the reported analysis
                          intact. Its limitation is stated, not hidden: the
                          definition is not the same in the last two cycles.
                          The two checks below exist to show whether that
                          matters.
    event_allcause        CHECK 1. Death from any cause, including deaths with
                          no cause code. Defined identically in every cycle.
                          The main analyses are repeated on it to see whether
                          the findings hold on an endpoint with one definition
                          throughout.
    event_nonaccidental   CHECK 2. The principal definition applied only where
                          it is exact: death with a recorded underlying cause
                          other than accidents (UCOD_LEADING 4), in NHANES
                          1999-2014 and NHANES III. It is missing (pd.NA), not
                          zero, for every 2015-2018 participant, so that a model
                          which asks for it across all ten cycles fails loudly
                          on the missing values instead of running quietly.
                          An analysis of this endpoint must therefore leave the
                          2015-2018 participants out of its cohort;
                          `endpoint_cycles` names the cycles it may keep.

    In NHANES III the principal endpoint and check 2 are the same thing,
    because accidents are identifiable for every participant.

Why accidents cannot be recovered in 2015-2018
    For those two cycles NCHS releases only three underlying-cause categories:
    heart disease (1), malignant neoplasms (2) and a residual (10). Accidental
    deaths are inside the residual and cannot be separated from it. Nothing
    here, or anywhere downstream, may infer accident status from code 10.

Interpretation of the non-accidental endpoint
    Coding accidental deaths as non-events censors them at the time of death.
    A Cox model on that indicator estimates a cause-specific hazard. It does not
    give the probability of non-accidental death: for that, accidental death is
    a competing event and a cumulative-incidence estimator is required.

    A death with no cause code is also censored under this endpoint. Its cause
    is unknown, so it is not an identified non-accidental death and cannot be
    counted as one. That is a limitation of the sensitivity, not a
    classification: 2 such deaths among NHANES 1999-2014 eligible participants
    and 13 in NHANES III.

Cause-specific analyses
    An analysis of death from one named cause may select deaths by testing the
    cause code for equality with that cause. That is the one place a module
    other than this one reads the cause code. It must restrict its COHORT, not
    only its events, to the cycles in which that cause is released as its own
    category, using `cause_identifiable`:
        heart disease (1), malignant neoplasms (2)   every cycle
        every other code, 3 to 10                    1999-2014 and NHANES III
    Without the restriction a 2015-2018 death from, say, stroke sits in the
    residual category, is counted as a non-event for the stroke analysis, and
    is censored as if the person had died of something else. The residual code
    10 is restricted as well, because in 2015-2018 it holds every cause except
    heart disease and cancer and is not the same category as in earlier cycles.
    The multiple-cause diabetes and hypertension flags are released for every
    cycle and need no restriction.

Assumptions and limits
    - `mortstat` is 1 for a linked death; anything else, including a missing
      value for a participant without linkage, is not an event. Participants
      without linkage must be removed by the caller's follow-up filter, exactly
      as before the revision.
    - Cause codes are the ten NCHS public-use leading-cause recodes, 1-10.
    - Which categories are released for which survey is taken from the NCHS
      codebook for the 2019 public-use linked mortality files ("For 2015-2018
      NHANES only, UCOD_LEADING=001, 002, and 010 are available"), and is
      checked against the data every time the endpoints are attached.

Outputs (columns added to the caller's frame)
    event_allcause, event_submitted_v1   int, 0/1
    event_nonaccidental                  nullable Int64, 0/1/NA
    accidental_death                     int, 1 for a death coded 4. Can only
                                         be 1 where accidents are identifiable.
    endpoint                             name of the configured endpoint
    endpoint_event                       copy of that endpoint's column; the one
                                         column analysis code should read

Decision record
    docs/endpoint.md

Usage
    from endpoints import ENDPOINT, attach_endpoints, expected_events
    from endpoints import cause_identifiable   # cause-specific analyses only
    from endpoints import endpoint_cycles      # cycles the endpoint is defined on
    from endpoints import endpoint_label, endpoint_events_label   # names for text
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


# The endpoint every module analyses unless a script names another explicitly.
#
# It is deliberately a constant in version-controlled code and not an
# environment variable or command-line flag: an endpoint chosen by hidden
# session state is how a rerun ends up on the wrong outcome without anyone
# noticing. Changing this line is an outcome-relevant amendment and requires the
# affected modules to be rerun from empty output directories.
#
# It is `event_submitted_v1`, the principal endpoint (see the module docstring
# for what that is and why). The two check endpoints are never configured here:
# a check is run in a separate copy of the analysis tree, so that its results
# cannot land in the folders that hold the principal results.
ENDPOINT = "event_submitted_v1"

EVENT_COLUMN = "endpoint_event"
ENDPOINTS = ("event_allcause", "event_nonaccidental", "event_submitted_v1")

ENDPOINT_LABELS = {
    "event_allcause": "all-cause mortality",
    "event_nonaccidental": "non-accidental mortality, 1999-2014",
    # The manuscript's term. Accidental deaths are censored at death wherever
    # the public-use files identify them (1999-2014, NHANES III); in 2015-2018
    # they cannot be identified and are counted. Figures and generated legends
    # take the term from here, so it is changed in one place.
    "event_submitted_v1": "accident-censored mortality",
}
# The same endpoints as a name for a COUNT of events, for table headings, axis
# annotations and generated legends ("2,494 accident-censored mortality
# events"). Kept beside the labels so the two cannot drift apart.
ENDPOINT_EVENT_LABELS = {
    "event_allcause": "deaths from any cause",
    "event_nonaccidental": "non-accidental deaths, 1999-2014",
    "event_submitted_v1": "accident-censored mortality events",
}
ENDPOINT_ROLES = {
    "event_submitted_v1": "principal",
    "event_allcause": "check 1: one definition in every cycle",
    "event_nonaccidental": "check 2: principal definition where it is exact",
}

ACCIDENT_CODE = 4
VALID_CAUSE_CODES = frozenset(range(1, 11))
# The only cause categories NCHS releases for NHANES 2015-2018.
RESTRICTED_RELEASE_CODES = frozenset({1, 2, 10})
# The categories that mean the same thing in every cycle. Code 10 is released in
# 2015-2018 but is not among them: there it is "everything except heart disease
# and cancer", which is a different category from the residual of earlier cycles.
CAUSES_IDENTIFIABLE_IN_ALL_CYCLES = frozenset({1, 2})

ACCIDENT_IDENTIFIABLE_CYCLES = frozenset({
    "1999-2000", "2001-2002", "2003-2004", "2005-2006",
    "2007-2008", "2009-2010", "2011-2012", "2013-2014",
})
CAUSE_RESTRICTED_CYCLES = frozenset({"2015-2016", "2017-2018"})
ALL_CYCLES = ACCIDENT_IDENTIFIABLE_CYCLES | CAUSE_RESTRICTED_CYCLES
# The same eight cycles under the name that fits a cause-specific analysis: the
# cycles for which all ten cause categories are released.
FULL_CAUSE_RELEASE_CYCLES = ACCIDENT_IDENTIFIABLE_CYCLES

# This file's own path, so that modules with resumable checkpoints can include
# it in their cache keys. A checkpoint computed under one endpoint must never be
# resumed under another.
ENDPOINT_MODULE_PATH = Path(__file__).resolve()

# Event counts that certify the two development cohorts. They are facts about
# the locked inputs, recorded so that a silent change in either the data or the
# definition stops the pipeline.
#   older_52_79: 11,497 panel-complete participants aged 52-79 with follow-up.
#                2,539 deaths, of which 45 are identifiable accidents.
#   adult_linked: 29,023 panel-complete adults with positive follow-up.
#                3,028 deaths, of which 106 are identifiable accidents.
#   older_52_79_1999_2014: the 9,362 members of older_52_79 examined in
#                1999-2014, the cohort of check 2. 2,447 deaths, of which 45 are
#                identifiable accidents. It is certified for the non-accidental
#                endpoint only, which is the one endpoint analysed on it.
_EXPECTED_EVENTS = {
    "older_52_79": {"event_allcause": 2_539, "event_submitted_v1": 2_494},
    "adult_linked": {"event_allcause": 3_028, "event_submitted_v1": 2_922},
    "older_52_79_1999_2014": {"event_nonaccidental": 2_402},
}


def canonical_cycle(values: pd.Series) -> pd.Series:
    """Return cycle labels as '1999-2000'; the import step writes '1999_2000'."""
    return values.astype("string").str.replace("_", "-", regex=False)


def expected_events(cohort: str, endpoint: str = ENDPOINT) -> int:
    """Certified event count for a named development cohort under an endpoint.

    Raises KeyError for a cohort/endpoint pair that has not been certified, so a
    new combination cannot pass a count check by default.
    """
    return _EXPECTED_EVENTS[cohort][endpoint]


def _first_letter_upper(text: str) -> str:
    return text[:1].upper() + text[1:]


def endpoint_label(endpoint: str = ENDPOINT, *, capitalised: bool = False) -> str:
    """The name of an endpoint as generated text should give it.

    "accident-censored mortality" for the principal endpoint. Every figure
    label, table heading and generated legend that names the endpoint takes it
    from here, and no script writes the name itself. Two things follow: the
    term is changed in one place, and an analysis run on another endpoint (a
    check) labels its outputs with that endpoint's name without any edit.

    `capitalised=True` gives the form for the start of a heading or a table
    row ("Accident-censored mortality").
    """
    label = ENDPOINT_LABELS[endpoint]
    return _first_letter_upper(label) if capitalised else label


def endpoint_events_label(endpoint: str = ENDPOINT, *, capitalised: bool = False) -> str:
    """The name for a count of an endpoint's events.

    "accident-censored mortality events" for the principal endpoint: the
    participants who died of any recorded cause, without the deaths identified
    as accidental. Use it wherever generated text counts events; use
    `endpoint_label` where it names the outcome.
    """
    label = ENDPOINT_EVENT_LABELS[endpoint]
    return _first_letter_upper(label) if capitalised else label


def endpoint_cycles(endpoint: str = ENDPOINT) -> tuple[str, ...]:
    """Survey cycles on which `endpoint` is defined, in chronological order.

    All ten cycles for the principal endpoint and for all-cause mortality. The
    eight cycles of 1999-2014 for the non-accidental endpoint, which does not
    exist where accidental deaths cannot be identified.

    This says where an endpoint is defined. It does not build a cohort: a
    module that analyses the configured endpoint restricts its own cohort to
    these cycles, so that repeating an analysis on the non-accidental endpoint
    (check 2) needs the one configured line changed and nothing else. For the
    two endpoints defined everywhere the restriction removes nobody.

    Removing two cycles changes more than the event definition: participants,
    events, cross-validation folds and the mix of follow-up times change with
    it. A result from the restricted cohort is a result about 1999-2014, and a
    difference from the ten-cycle result cannot be put down to the endpoint.
    """
    if endpoint not in ENDPOINTS:
        raise ValueError(f"unknown endpoint: {endpoint!r}")
    cycles = (ACCIDENT_IDENTIFIABLE_CYCLES if endpoint == "event_nonaccidental"
              else ALL_CYCLES)
    return tuple(sorted(cycles))


def _accident_identifiable(frame: pd.DataFrame, cycle_col: str) -> pd.Series:
    """True where the participant's cycle releases all ten cause categories.

    Those are the cycles in which accidents, and every other cause apart from
    heart disease and cancer, have a category of their own.
    """
    cycle = canonical_cycle(frame[cycle_col])
    linked = frame["mortstat"].notna()
    unknown = linked & ~cycle.isin(list(ALL_CYCLES)).fillna(False)
    if unknown.any():
        examples = sorted(cycle[unknown].dropna().unique())[:5]
        raise AssertionError(
            "mortality-linked participants have an unrecognised or missing "
            f"survey cycle: {examples}"
        )
    return cycle.isin(list(ACCIDENT_IDENTIFIABLE_CYCLES)).fillna(False)


def cause_identifiable(
    frame: pd.DataFrame, cause_code: int, *, cycle_col: str = "cycle"
) -> pd.Series:
    """True for participants whose cycle releases `cause_code` as its own category.

    A cause-specific analysis must keep only these participants. Restricting
    the events alone is not enough: a participant from a cycle that does not
    release the cause would stay in the risk set, and if they died of that
    cause the death would be recorded under the residual code and treated as
    censoring by another cause. The estimate would then mix cycles where the
    event can be seen with cycles where it cannot.

    Heart disease (1) and malignant neoplasms (2) are identifiable in every
    cycle. Every other code, including the residual 10, is identifiable only
    where all ten categories are released (see the module docstring).

    The frame must contain `mortstat` and the cycle column. This is for
    continuous NHANES; NHANES III releases all ten categories throughout and
    needs no restriction.
    """
    code = int(cause_code)
    if code not in VALID_CAUSE_CODES:
        raise ValueError(f"not a UCOD_LEADING category: {cause_code!r}")
    # Called for every code so that an unrecognised cycle is always an error,
    # including for the two causes that need no restriction.
    all_categories_released = _accident_identifiable(frame, cycle_col)
    if code in CAUSES_IDENTIFIABLE_IN_ALL_CYCLES:
        return pd.Series(True, index=frame.index)
    return all_categories_released


def attach_endpoints(
    frame: pd.DataFrame,
    *,
    cycle_col: str = "cycle",
    accident_identifiable: bool | pd.Series | None = None,
    endpoint: str = ENDPOINT,
) -> pd.DataFrame:
    """Add the endpoint columns to `frame` in place and return it.

    Parameters
    ----------
    frame
        Must contain `mortstat` (1 = linked death) and `ucod_leading` (numeric
        NCHS leading-cause recode 1-10, missing for survivors and for deaths
        with no recorded cause).
    cycle_col
        Column holding the NHANES cycle. Used to decide where accidents are
        identifiable. Ignored when `accident_identifiable` is given.
    accident_identifiable
        Pass True for a cohort whose mortality file releases all ten cause
        categories for every participant (NHANES III). A death there can still
        have no cause code: it is an all-cause event and is not a
        non-accidental event. Leave as None for continuous NHANES so the
        decision is made per cycle.
    endpoint
        The endpoint copied into `endpoint_event`. Defaults to the configured
        principal endpoint; a sensitivity analysis names another explicitly.

    Fails if a cause code is outside 1-10, if a 2015-2018 death carries a code
    NCHS does not release for those cycles, or if `endpoint` is undefined for
    any row (which is what happens when the non-accidental endpoint is requested
    for a cohort that still contains 2015-2018 participants).
    """
    if endpoint not in ENDPOINTS:
        raise ValueError(f"unknown endpoint: {endpoint!r}")

    # `.fillna(False)` makes the comparisons behave the same for NumPy float
    # columns (missing = NaN) and pandas nullable integers (missing = pd.NA).
    died = frame["mortstat"].eq(1).fillna(False).astype(bool)
    code = frame["ucod_leading"]
    has_code = code.notna().astype(bool)
    accident_code = code.eq(ACCIDENT_CODE).fillna(False).astype(bool)

    bad_code = has_code & ~code.isin(list(VALID_CAUSE_CODES))
    if bad_code.any():
        raise AssertionError(
            f"unexpected UCOD_LEADING values: {sorted(code[bad_code].unique())[:5]}"
        )

    if accident_identifiable is None:
        identifiable = _accident_identifiable(frame, cycle_col)
    elif isinstance(accident_identifiable, pd.Series):
        identifiable = accident_identifiable.astype(bool)
    else:
        identifiable = pd.Series(bool(accident_identifiable), index=frame.index)

    # A data-definition guard, not a transformation: if NCHS ever released fuller
    # cause coding for the restricted cycles, the premise of this module would
    # have changed and the definition would need revisiting.
    restricted_death = died & has_code & ~identifiable
    leaked = restricted_death & ~code.isin(list(RESTRICTED_RELEASE_CODES))
    if leaked.any():
        raise AssertionError(
            "a cause-restricted cycle contains a cause code outside {1, 2, 10}; "
            "the public-use release has changed and AD-001 must be revisited"
        )

    non_accidental = died & has_code & ~accident_code

    frame["event_allcause"] = died.astype(int)
    # Same truth table as the pre-revision `disease_death`:
    # mortstat == 1 & ucod_leading recorded & ucod_leading != 4.
    frame["event_submitted_v1"] = non_accidental.astype(int)
    frame["event_nonaccidental"] = (
        non_accidental.astype("Int64").where(identifiable, pd.NA)
    )
    frame["accidental_death"] = (died & accident_code).astype(int)

    chosen = frame[endpoint]
    if chosen.isna().any():
        raise AssertionError(
            f"{endpoint} is undefined for {int(chosen.isna().sum())} rows. "
            "Restrict the cohort to cycles where accidents are identifiable "
            "before requesting the non-accidental endpoint."
        )
    frame["endpoint"] = endpoint
    frame[EVENT_COLUMN] = chosen.astype(int)
    return frame
