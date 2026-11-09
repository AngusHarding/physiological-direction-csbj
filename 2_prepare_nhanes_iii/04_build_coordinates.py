"""Transform, standardise, and build HD and FI coordinates for NHANES III.

SPECIFICATION: SPECIFICATION.md sections 6, 6a and 6b (Decisions 7, 8, 9).

Outcome-blind. No mortality variable is read.

Two coordinate systems, exactly as in the main analysis:

  HD  moment-based.  z = (T(x) - mu_s) / sigma_s, with mu and sigma the mean and
      sample SD of the SAME-SEX young-healthy reference. HD is then Mahalanobis
      distance from that centroid using the reference covariance.

  FI  robust.  z_r = m * (T(x) - M_s) / S_s, with M the reference median and
      S = IQR / 1.349 the normal-consistent robust scale. The adverse excursion
      follows the locked directionality; severity IS the excursion, with no dead
      zone and no cap (main FI_METHODS.md Amendment 1).

DECISION 9 in practice: `adverse_direction` is a clinical prior and is imported
unchanged from the main analysis. `orientation` is NOT imported -- it encodes
whether the SELECTED TRANSFORMATION reverses the raw-value axis, and NHANES III
re-selects transformations, so it is recomputed here. Importing it would silently
invert every marker whose transformation differs between eras, and ten do.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from nhanes3_locations import CONFIG, MODERN_CONFIG, OUTPUTS  # public_code: folders and write guard
OUT = OUTPUTS  # public_code

ROBUST_SCALE_DIVISOR = 1.349   # Phi^-1(0.75) - Phi^-1(0.25)


def transform(x: np.ndarray, kind: str) -> np.ndarray:
    if kind == "identity":
        return x
    if kind == "sqrt":
        return np.sqrt(x)
    if kind == "cbrt":
        return np.cbrt(x)
    if kind == "log":
        return np.log(x)
    if kind == "reciprocal":
        return -1.0 / x
    raise ValueError(kind)


def orientation_for(kind: str) -> int:
    """+1 unless the transformation reverses the raw-value axis.

    None of the transformations above reverses it: the reciprocal is computed as
    -1/x precisely so that it rises with the raw value. The submitted version
    returned -1 for the reciprocal, which reversed the frailty-index coordinate of
    AST, HbA1c, RDW and systolic pressure a second time, so the adverse-tail rule
    scored their LOW tail. Corrected 15 October 2026 (AD-009 amendment 1; defect
    and effect measured in revision note N3a, reviewed R013). HD coordinates do
    not use this multiplier and are unchanged.
    """
    if kind not in ("identity", "sqrt", "cbrt", "log", "reciprocal"):
        raise ValueError(kind)
    return 1


def main() -> None:
    panel = pd.read_csv(OUT / "nhanes3_panel_raw.csv")
    markers = pd.read_csv(OUT / "nhanes3_panel_markers.csv").marker.tolist()
    selected = pd.read_csv(OUT / "transformation_selected.csv")
    reference = pd.read_csv(OUT / "reference_members.csv")
    mapping = pd.read_csv(CONFIG / "nhanes3_marker_map.csv")  # public_code

    # adverse_direction imported unchanged from the main analysis (Decision 9)
    fi_spec = pd.read_csv(MODERN_CONFIG / "fi_scoring_spec.csv")  # public_code
    link = dict(zip(mapping.marker, mapping.main_panel_var))
    link.update({"ACR": "ACR", "eGFR": "eGFR", "WHtR": "WHtR"})
    adverse = dict(zip(fi_spec.marker, fi_spec.adverse_direction))
    direction = {m: adverse.get(link.get(m)) for m in markers}
    unmapped = [m for m, v in direction.items() if v is None]
    if unmapped:
        raise AssertionError(f"no adverse_direction for: {unmapped}")

    tf = dict(zip(selected.marker, selected.selected_transform))
    panel["is_reference"] = panel.SEQN.isin(reference.SEQN).astype(int)
    print(f"cohort {len(panel):,} | reference {int(panel.is_reference.sum()):,} "
          f"| markers {len(markers)}")

    # ---- transform ---------------------------------------------------------
    T = pd.DataFrame({"SEQN": panel.SEQN, "HSSEX": panel.HSSEX,
                      "HSAGEIR": panel.HSAGEIR, "is_reference": panel.is_reference})
    for m in markers:
        T[m] = transform(panel[m].to_numpy(float), tf[m])
    if not np.isfinite(T[markers].to_numpy(float)).all():
        raise AssertionError("non-finite value after transformation")

    # ---- sex-specific reference parameters ---------------------------------
    params, hd_z, fi_def = [], {}, {}
    hd_frames, fi_frames = [], []
    for sex_code, sex in ((1, "men"), (2, "women")):
        ref = T[(T.HSSEX == sex_code) & T.is_reference.eq(1)]
        grp = T[T.HSSEX == sex_code]
        mu = ref[markers].mean()
        sd = ref[markers].std(ddof=1)
        med = ref[markers].median()
        iqr = ref[markers].quantile(.75) - ref[markers].quantile(.25)
        robust = iqr / ROBUST_SCALE_DIVISOR
        if (sd <= 0).any() or (robust <= 0).any():
            raise AssertionError(f"degenerate reference scale in {sex}")

        z = (grp[markers] - mu) / sd                       # HD coordinates
        zr = pd.DataFrame(index=grp.index)                  # FI robust coordinates
        for m in markers:
            zr[m] = orientation_for(tf[m]) * (grp[m] - med[m]) / robust[m]

        deficits = pd.DataFrame(index=grp.index)
        for m in markers:
            d = direction[m]
            if d == "high_only":
                e = zr[m].clip(lower=0)
            elif d == "low_only":
                e = (-zr[m]).clip(lower=0)
            elif d == "two_sided":
                e = zr[m].abs()
            else:
                raise ValueError(f"{m}: unknown adverse_direction {d!r}")
            deficits[m] = e                                 # severity IS the excursion

        hd_frames.append(z.assign(SEQN=grp.SEQN.values, HSSEX=sex_code))
        fi_frames.append(deficits.assign(SEQN=grp.SEQN.values, HSSEX=sex_code,
                                         **{f"signed_{m}": np.sign(zr[m]) * deficits[m]
                                            for m in markers}))
        for m in markers:
            params.append({"sex": sex, "marker": m, "transform": tf[m],
                           "orientation": orientation_for(tf[m]),
                           "adverse_direction": direction[m],
                           "ref_mean": mu[m], "ref_sd": sd[m],
                           "ref_median": med[m], "ref_robust_scale": robust[m]})

    hd = pd.concat(hd_frames).sort_values("SEQN")
    fi = pd.concat(fi_frames).sort_values("SEQN")
    pd.DataFrame(params).round(6).to_csv(OUT / "reference_parameters.csv", index=False)

    # ---- scores ------------------------------------------------------------
    scores = pd.DataFrame({"SEQN": hd.SEQN.values, "HSSEX": hd.HSSEX.values})
    scores["FI"] = fi[markers].mean(axis=1).values
    hd_values = []
    for sex_code in (1, 2):
        sel = hd[hd.HSSEX == sex_code]
        ref_idx = T.set_index("SEQN").loc[sel.SEQN, "is_reference"].to_numpy() == 1
        Z = sel[markers].to_numpy(float)
        cov = np.cov(Z[ref_idx], rowvar=False, ddof=1)
        inv = np.linalg.inv(cov)
        cond = np.linalg.cond(cov)
        quad = np.einsum("ij,jk,ik->i", Z, inv, Z)
        hd_values.append(pd.DataFrame({"SEQN": sel.SEQN.values,
                                       "HD": np.sqrt(np.maximum(quad, 0))}))
        print(f"  {['','men','women'][sex_code]:5s} reference covariance: "
              f"condition number {cond:,.0f}")
    scores = scores.merge(pd.concat(hd_values), on="SEQN")

    hd.to_csv(OUT / "coordinates_hd.csv", index=False)
    fi.to_csv(OUT / "coordinates_fi.csv", index=False)
    scores.to_csv(OUT / "scores.csv", index=False)

    print(f"\nHD  median {scores.HD.median():.2f}  IQR "
          f"{scores.HD.quantile(.25):.2f}-{scores.HD.quantile(.75):.2f}")
    print(f"FI  median {scores.FI.median():.3f}  IQR "
          f"{scores.FI.quantile(.25):.3f}-{scores.FI.quantile(.75):.3f}")
    print(f"  main analysis FI for comparison: median 1.07, IQR 0.90-1.26")
    ref_scores = scores[scores.SEQN.isin(reference.SEQN)]
    print(f"\nreference participants sit lower, as they must:")
    print(f"  HD  reference median {ref_scores.HD.median():.2f} vs cohort {scores.HD.median():.2f}")
    print(f"  FI  reference median {ref_scores.FI.median():.3f} vs cohort {scores.FI.median():.3f}")
    print(f"\nwritten -> {OUT}")


if __name__ == "__main__":
    main()
