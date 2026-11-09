"""Parse, harmonise, derive and QC the NHANES III biomarker panel.

SPECIFICATION: SPECIFICATION.md sections 4, 6 and 6a (Decisions 1, 6, 7).

Outcome-blind. No mortality variable is read anywhere in this script.

Order of operations is fixed by Decision 6:
    parse -> clear fill codes -> plausibility -> distributions
    -> side-by-side against continuous NHANES -> transformation selection
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from nhanes3_locations import OUTPUTS, RAW_NHANES3  # public_code: folders and write guard
from nhanes3_io import (  # noqa: E402
    apply_plausibility, distribution_report, load_marker_map, range_report,
    read_fixed_width, read_layout,
)

OUT = OUTPUTS  # public_code
AGE_MIN, AGE_MAX = 18, 79


def selvin_standardise_creatinine(scr: pd.Series) -> pd.Series:
    """Calibrate NHANES III creatinine to an IDMS-traceable scale.

    NHANES III creatinine predates assay standardisation and reads high. CKD-EPI
    2021 assumes standardised creatinine, and eGFR is NONLINEAR in creatinine, so
    the correction cannot be absorbed by later reference-relative standardisation
    the way a pure location shift would be.

    standardised = 0.960 * NHANES_III - 0.184

    CITATION MUST BE VERIFIED before write-up (Selvin et al., calibration of
    NHANES III serum creatinine). Recorded here from working knowledge; the
    coefficients should not be trusted without checking the source. An
    uncalibrated eGFR is retained alongside as a sensitivity.
    """
    return 0.960 * scr - 0.184


def egfr_ckdepi_2021(creatinine, age, sex):
    """Race-free CKD-EPI 2021. Identical to the main pipeline's implementation
    (data_preparation/src/03_derive/derive_kidney_features.py). NHANES sex codes
    1 = male, 2 = female."""
    creatinine = np.asarray(creatinine, dtype=float)
    female = np.asarray(sex) == 2
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    ratio = creatinine / kappa
    estimate = (142.0 * np.minimum(ratio, 1.0) ** alpha
                * np.maximum(ratio, 1.0) ** -1.200
                * 0.9938 ** np.asarray(age, dtype=float))
    return np.where(female, estimate * 1.012, estimate)


def load_source(spec: pd.DataFrame) -> pd.DataFrame:
    """Parse lab and exam, clear fill codes, apply plausibility, join on SEQN."""
    joined = None
    audit = []
    for src, sas, dat in (("lab", "lab.sas", "lab.dat"),
                          ("exam", "exam.sas", "exam.dat")):
        sub = spec[spec.source_file.eq(src)]
        layout = read_layout(RAW_NHANES3 / sas)  # public_code
        wanted = [v for v in dict.fromkeys(
            ["SEQN", "HSAGEIR", "HSSEX", "SDPPHASE", "SDPPSU6", "SDPSTRA6"]
            + sub.nhanes3_var.tolist()) if v in layout]
        frame = read_fixed_width(RAW_NHANES3 / dat, layout, wanted)  # public_code
        frame = frame[frame.HSAGEIR.between(AGE_MIN, AGE_MAX)].copy()
        frame, set_missing = apply_plausibility(frame, sub)
        if len(set_missing):
            set_missing.insert(0, "source_file", src)
            audit.append(set_missing)
        report = range_report(frame, sub)
        report.to_csv(OUT / f"qc_range_{src}.csv", index=False)
        distribution_report(frame, sub).to_csv(
            OUT / f"qc_distribution_raw_{src}.csv", index=False)
        if not report.passes.all():
            raise AssertionError(
                f"{src}: range assertions failed for "
                f"{report.loc[~report.passes, 'marker'].tolist()} "
                "(Decision 6: a range failure halts the build)")
        print(f"  {src:5s} {len(frame):6d} adults | {int(report.passes.sum())}/"
              f"{len(report)} markers pass range assertions")
        keep = [c for c in frame.columns if c not in ("HSAGEIR", "HSSEX",
                "SDPPHASE", "SDPPSU6", "SDPSTRA6")] if joined is not None else list(frame.columns)
        joined = frame[keep] if joined is None else joined.merge(
            frame[keep], on="SEQN", how="inner", validate="one_to_one")
    if audit:
        pd.concat(audit).to_csv(OUT / "qc_implausible_set_missing.csv", index=False)
    return joined


def derive(panel: pd.DataFrame, spec: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Rename to marker names and build the three derived markers."""
    rename = dict(zip(spec.nhanes3_var, spec.marker))
    d = panel.rename(columns=rename)

    # ACR (mg/g) = urine albumin (ug/mL = mg/L) / urine creatinine (g/L),
    # where creatinine mg/dL -> g/L is a divide by 100. Same convention as the
    # main pipeline's Alb_Cr.
    d["ACR"] = d.urine_albumin / (d.urine_creatinine / 100.0)

    d["creatinine_std"] = selvin_standardise_creatinine(d.creatinine)
    d["eGFR"] = egfr_ckdepi_2021(d.creatinine_std, d.HSAGEIR, d.HSSEX)
    d["eGFR_uncalibrated"] = egfr_ckdepi_2021(d.creatinine, d.HSAGEIR, d.HSSEX)

    d["WHtR"] = d.waist / d.height

    direct = spec.loc[spec.in_primary_set.eq(1)
                      & ~spec.marker.isin(["urine_albumin", "urine_creatinine",
                                           "creatinine", "waist", "height"]),
                      "marker"].tolist()
    markers = direct + ["ACR", "eGFR", "WHtR"]
    return d, markers


def main() -> None:
    OUT.mkdir(exist_ok=True)
    spec = load_marker_map()
    print("parsing source files")
    panel = load_source(spec)
    panel, markers = derive(panel, spec)
    print(f"\npanel markers: {len(markers)}")

    complete = panel[panel[markers].notna().all(axis=1)].copy()
    print(f"adults {AGE_MIN}-{AGE_MAX}          : {len(panel):,}")
    print(f"complete-case on {len(markers)} markers: {len(complete):,}")
    print(f"  age >=52                : {int(complete.HSAGEIR.ge(52).sum()):,}")

    complete.to_csv(OUT / "nhanes3_panel_raw.csv", index=False)
    pd.DataFrame({"marker": markers}).to_csv(OUT / "nhanes3_panel_markers.csv", index=False)

    dist = []
    from scipy import stats
    for m in markers:
        x = complete[m].dropna()
        k2, p = stats.normaltest(x)
        dist.append({"marker": m, "n": len(x), "mean": x.mean(), "median": x.median(),
                     "sd": x.std(ddof=1), "min": x.min(), "max": x.max(),
                     "skew": stats.skew(x), "excess_kurtosis": stats.kurtosis(x),
                     "dagostino_k2": k2, "dagostino_p": p})
    pd.DataFrame(dist).round(4).to_csv(OUT / "qc_distribution_panel_raw.csv", index=False)
    print(f"\nwritten -> {OUT}")


if __name__ == "__main__":
    main()
