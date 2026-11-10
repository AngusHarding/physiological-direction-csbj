# =============================================================================
# READ ME FIRST: 2_prepare_nhanes_iii/05_extract_results_variables.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Copies the NHANES III questionnaire and examination variables
#               the Results use (Table 1, Results 5-7), and the NCHS linked
#               mortality file, into prepared tables. Derives nothing; NHANES
#               III's 8- and 9-fill codes become missing.
#
# METHODS       Data import only. The endpoint is applied in 3_results/data.py.
# =============================================================================
"""Copy the NHANES III variables the Results use, and the linked mortality file, into three tables.

Steps 1-4 build the marker panel. The Results also need some NHANES III
questionnaire and examination variables (education, income, smoking, physical
activity, conditions, BMI, alcohol, blood pressure, spirometry, bone density,
bioimpedance, balance) and the NCHS linked mortality file. This step copies
them out of the raw fixed-width files, so the Results read prepared tables
only, like every other input.

Nothing is derived here. Each value is read as the raw file holds it, with one
rule: NHANES III fills a missing number with 8s or 9s across the whole field
(88888 in a five-character field), and those values are set to missing.

Makes, in outputs/:
  results_adult_variables.csv   SEQN and 28 adult-questionnaire variables
  results_exam_variables.csv    SEQN and 17 examination variables
  mortality_linked.csv          the NCHS 2019 linked mortality file: SEQN, eligibility,
                                vital status, leading underlying cause, months of
                                follow-up from the examination
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from nhanes3_locations import OUTPUTS, RAW_NHANES3  # folders and write guard

# Which section uses each variable is given on the right.
ADULT = (["HSAGEIR", "HSSEX", "DMARETHN", "HFA8R", "DMPPIR", "HAR1", "HAR3", "HAT1S"]      # Results 1, 6, 7
         + ["HAT2", "HAT4", "HAT6", "HAT8", "HAT10", "HAT12", "HAT14", "HAT16", "HAT18"]   # Results 1: activity
         + ["HAC1A", "HAC1C", "HAC1D", "HAC1F", "HAC1G", "HAC1J", "HAC1K", "HAC1L",
            "HAC1M", "HAC1N", "HAC1O"])                                                     # Results 1: conditions
EXAM = (["BMPBMI", "MAPE1", "MAPE2", "MAPE3S", "MAPE4"]                                     # Results 1: BMI, alcohol
        + ["PEP6H1", "PEP6I1", "PEP6H3", "PEP6I3"]                                          # blood pressure (data.py)
        + ["BMPHT", "SPPFEV1", "SPPFVC", "BDPFNBMD", "BDPTOBMC", "PEP12A1", "PEP12B1",
           "PFPTTAND"])                                                                     # Results 7


def read_layout(sas_file: Path) -> dict[str, tuple[int, int]]:
    """Column positions {NAME: (first, last)} from a NHANES III SAS layout file.

    A NHANES III .dat file is plain text with one participant per line and no
    separators; the SAS file lists which characters hold which variable.
    """
    text = sas_file.read_text(encoding="latin-1")
    upper = text.upper()
    block = text[upper.find("INPUT") + 5: upper.find("LABEL")]
    return {m.group(1).upper(): (int(m.group(2)), int(m.group(3) or m.group(2)))
            for m in re.finditer(r"^\s*([A-Z0-9_]+)\s+(\d+)(?:-(\d+))?\s*$", block, re.M | re.I)}


def read_nhanes3(file: str, variables: list[str]) -> pd.DataFrame:
    """Selected variables from a NHANES III .dat file, with 8- and 9-fills set to missing."""
    layout = read_layout(RAW_NHANES3 / f"{file}.sas")
    columns = {v: [] for v in variables}
    with open(RAW_NHANES3 / f"{file}.dat", encoding="latin-1") as handle:
        for line in handle:
            for v in variables:
                first, last = layout[v]
                columns[v].append(line[first - 1:last].strip())
    frame = pd.DataFrame(columns)
    for v in variables:
        first, last = layout[v]
        values = pd.to_numeric(frame[v], errors="coerce")
        for fill in ("8", "9"):
            values = values.mask(values.eq(float(fill * (last - first + 1))))
        frame[v] = values
    return frame


def read_mortality() -> pd.DataFrame:
    """The NCHS 2019 linked mortality file for NHANES III, every row (eligible or not)."""
    rows = []
    with open(RAW_NHANES3 / "NHANES_III_MORT_2019_PUBLIC.dat", encoding="latin-1") as fh:
        for line in fh:
            rows.append({"SEQN": int(line[0:5]), "eligstat": int(line[14:15]),
                         "mortstat": line[15:16].strip(), "ucod_leading": line[16:19].strip(),
                         "permth_exm": line[45:48].strip()})
    m = pd.DataFrame(rows)
    for column in ("mortstat", "ucod_leading", "permth_exm"):
        m[column] = pd.to_numeric(m[column], errors="coerce")
    return m


def main() -> None:
    for file, variables in (("adult", ADULT), ("exam", EXAM)):
        frame = read_nhanes3(file, ["SEQN"] + variables)
        assert frame.SEQN.is_unique, f"{file}: participant numbers repeat"
        frame.to_csv(OUTPUTS / f"results_{file}_variables.csv", index=False)
        print(f"results_{file}_variables.csv: {len(frame):,} participants, {len(variables)} variables")
    mortality = read_mortality()
    assert mortality.SEQN.is_unique, "mortality file: participant numbers repeat"
    mortality.to_csv(OUTPUTS / "mortality_linked.csv", index=False)
    print(f"mortality_linked.csv: {len(mortality):,} rows, {int(mortality.eligstat.eq(1).sum()):,} eligible for linkage")


if __name__ == "__main__":
    main()
