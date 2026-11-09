"""Build the NHANES III young-healthy reference cohort.

SPECIFICATION: SPECIFICATION.md Decision 8.

Outcome-blind. No mortality variable is read.

The reference defines the origin and unit of every coordinate downstream, so it
is built from questionnaire and examination data only, and separately by sex.
Two criteria of the main analysis cannot be replicated and are documented in
Decision 8: no kidney-failure item exists (immaterial, 0.06% of under-30s have
eGFR < 60), and alcohol is recorded as frequency rather than volume.

Structural non-coverage is not disqualifying. The functional-limitation items are
age-targeted -- answered by 53% of under-60s against 97% of over-60s -- so an
unanswered item means the question was not put, not that the participant failed
it. This mirrors the main analysis, whose availability-aware reference of 1,238
differs from its permissive 1,241 for exactly this reason.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from nhanes3_locations import OUTPUTS, RAW_NHANES3  # public_code: folders and write guard
from nhanes3_io import read_fixed_width, read_layout  # noqa: E402

OUT = OUTPUTS  # public_code

# HAC1 battery, matching the main analysis's thirteen-condition list.
# HAC1E (asthma) is NOT here: the main reference retains asthmatics, and testing
# showed exclusion moves every immune reference parameter by under 1%.
# HAC1B is a follow-up on arthritis type, HAC1H hay fever and HAC1I cataracts are
# not age-related chronic disease in the main sense.
CONDITIONS = ["HAC1A", "HAC1C", "HAC1D", "HAC1F", "HAC1G",
              "HAC1J", "HAC1K", "HAC1L", "HAC1M", "HAC1N", "HAC1O"]
ALCOHOL = ["HAN6HS", "HAN6IS", "HAN6JS"]
FUNCTIONAL = ["HAH1", "HAH2", "HAH8"]
ALCOHOL_MAX_OCCASIONS_PER_MONTH = 30


def main() -> None:
    layout = read_layout(RAW_NHANES3 / "adult.sas")  # public_code
    wanted = ["SEQN", "HAR1", "HAD1", "HAF10", "HAG11"] + CONDITIONS + ALCOHOL + FUNCTIONAL
    adult = read_fixed_width(RAW_NHANES3 / "adult.dat", layout, wanted)  # public_code

    panel = pd.read_csv(OUT / "nhanes3_panel_raw.csv")
    d = panel.merge(adult, on="SEQN", how="left", validate="one_to_one")
    occasions = d[ALCOHOL].fillna(0).sum(axis=1)

    steps, mask = [], d.HSAGEIR.lt(30)
    def record(label):
        steps.append({"criterion": label, "remaining": int(mask.sum())})
    record("age < 30")

    mask &= d.bmi.between(18.5, 24.9);           record("BMI 18.5-24.9")
    mask &= d.HAR1.eq(2);                        record("never smoker")
    mask &= occasions.lt(ALCOHOL_MAX_OCCASIONS_PER_MONTH)
    record(f"< {ALCOHOL_MAX_OCCASIONS_PER_MONTH} drinking occasions/month")
    mask &= d.HAD1.ne(1);                        record("no diabetes")
    mask &= d.HAF10.ne(1);                       record("no myocardial infarction")
    mask &= d.HAG11.ne(1);                       record("no osteoporosis")
    for c in CONDITIONS:
        mask &= d[c].ne(1)
    record("none of 11 chronic conditions")
    # Availability-aware: an unanswered functional item is a question not asked.
    impaired = (d[FUNCTIONAL].gt(1).any(axis=1)).fillna(False)
    mask &= ~impaired;                           record("no functional impairment")

    reference = d.loc[mask].copy()
    reference["is_reference"] = 1

    cascade = pd.DataFrame(steps)
    cascade["removed"] = cascade.remaining.shift().sub(cascade.remaining).fillna(0).astype(int)
    cascade.to_csv(OUT / "reference_cascade.csv", index=False)
    reference[["SEQN", "HSSEX", "HSAGEIR", "is_reference"]].to_csv(
        OUT / "reference_members.csv", index=False)

    print(cascade.to_string(index=False))
    men = int((reference.HSSEX == 1).sum())
    women = int((reference.HSSEX == 2).sum())
    print(f"\nyoung-healthy reference: {len(reference)}  ({men} men, {women} women)")
    print(f"main analysis comparator: 1,238  (547 men, 691 women)")
    markers = pd.read_csv(OUT / "nhanes3_panel_markers.csv").marker.tolist()
    print(f"participants per marker : {len(reference)/len(markers):.1f} "
          f"(main analysis 32.6)")
    if min(men, women) <= len(markers):
        raise AssertionError(
            "reference too small in one sex to invert a covariance matrix")
    print(f"\nwritten -> {OUT}")


if __name__ == "__main__":
    main()
