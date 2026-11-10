# =============================================================================
# READ ME FIRST: measures.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Calculates physiological distance (DA, DM, DE) for each
#               participant, and the coordinates under the three references
#               (SS, SC, PC) compared in Results 2. Nothing here uses mortality.
#
# METHODS       "Physiological distance and the mortality score":
#               "Young-healthy reference" and "Distances".
#
# USED BY       Results 1 (s1_sex_gap.py): DA, DM, DE; Fig. 1c,d; Supp. Table
#                 S1.3.
#               Results 2 (s2_reference.py): SS/SC/PC coordinates (Fig. 2c,d;
#                 Supp. Table S2.2) and DA under each reference (Table 2).
#               Results 3 (s3_healthy_range.py): DA as a single score; Supp.
#                 Table S3.6.
#
# PSEUDOCODE
#   DA(person) = mean over markers of |x|       x = adverse-only excursion,
#                                                built in data preparation
#   DM(person) = sqrt( z' S^-1 z )               S = covariance of z among the
#                                                same-sex young-healthy members
#   DE(person) = sqrt( sum of z^2 )
#   coordinate = (transformed value - centre) / scale, where
#       SS: centre = own-sex mean,  scale = own-sex SD        (primary)
#       SC: centre = own-sex mean,  scale = SD pooled over both sexes
#       PC: centre = pooled mean,   scale = pooled SD
# =============================================================================
"""Physiological distance from the young-healthy reference: three ways to measure it.

None of these is fitted to mortality. Each summarises, in one number per
participant, how far that participant's markers lie from their own sex's
young-healthy reference.

  adverse distance (DA)      the average, over markers, of how far each marker
                             lies beyond its young-healthy median in the
                             direction that is bad for health (x coordinates).
  Mahalanobis distance (DM)  distance in z coordinates that allows for how
                             markers vary together in the healthy reference: two
                             markers that usually move together count less than
                             twice when both are high.
  Euclidean distance (DE)    the straight-line length of the z profile.

The references are fixed by data preparation (the young-healthy members of
each sex); DM uses each sex's reference covariance. The three distances are
descriptions of displacement, not risk scores; Results 1 and 2 test how they
relate to mortality. SS, SC and PC (Results 2) are the three ways of placing
women and men on the z scale.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# -- CALCULATION: adverse distance (DA) --------------------------------------
# Methods, "Distances": "Adverse distance (DA): the mean excursion beyond the
#   reference median in each marker's prespecified adverse direction, on the
#   transformed scale and in units of the reference interquartile range divided
#   by 1.349."
# The excursions x are built in data preparation
# (1_prepare_nhanes/8_adverse_distance/), so here DA is only their mean.
# Results 1: "women's adverse distance exceeded men's by 0.19 standard
#   deviations (SD) in both eras".
def adverse_distance(x: pd.DataFrame) -> np.ndarray:
    """DA: mean absolute adverse excursion across markers.

    The x coordinates are signed (a marker too high is positive, too low is
    negative, and zero on the healthy side), so the size of the excursion is
    the absolute value.
    """
    return np.abs(x.to_numpy(float)).mean(axis=1)


# -- CALCULATION: Mahalanobis distance (DM) ----------------------------------
# Methods, "Distances": "Mahalanobis distance (DM): the basis of homeostatic
#   dysregulation [21, 46]. It discounts excursions along marker combinations
#   that vary together in the reference."
# How: DM = sqrt(z' S^-1 z). S is the covariance matrix of z among the
#   young-healthy members of the participant's own sex, so markers that usually
#   rise together in healthy people count less than twice when both are high.
# Results 1: "Mahalanobis distance did not differ by sex in NHANES 1999-2018
#   (0.01 SD; -0.02 to 0.05) and was 0.11 SD lower in women in NHANES III".
def mahalanobis_distance(z: pd.DataFrame, female: pd.Series, ref_z: dict) -> np.ndarray:
    """DM: sqrt(z' S^-1 z), with S the covariance of the z coordinates of the
    participant's own sex's young-healthy reference members."""
    out = np.full(len(z), np.nan)
    for sex, is_female in ((1, 0), (2, 1)):
        precision = np.linalg.inv(np.cov(ref_z[sex], rowvar=False, ddof=1))
        rows = female.eq(is_female).to_numpy()
        v = z.to_numpy(float)[rows]
        out[rows] = np.sqrt(np.einsum("ij,jk,ik->i", v, precision, v))
    return out


# -- CALCULATION: Euclidean distance (DE) ------------------------------------
# Methods, "Distances": "Euclidean distance (DE): the straight-line length of
#   the standardised profile."
# How: DE = sqrt(z1^2 + z2^2 + ... + zk^2). DE is also the total that q and R
#   split: DE^2 = q^2 + R^2 (Methods, "Mortality score").
# Results 1 (s1_sex_gap.py) calculates DE alongside DA and DM.
def euclidean_distance(z: pd.DataFrame) -> np.ndarray:
    """DE: length of the z profile."""
    return np.linalg.norm(z.to_numpy(float), axis=1)


# ---------------------------------------------------------------------------
# References: where each marker starts from, and what unit it is measured in
# ---------------------------------------------------------------------------
# A reference turns a transformed marker value t into a coordinate
#     (t - centre) / scale
# The paper compares three choices (Results 2):
#   SS  sex-specific centre and scale (primary): each sex's own young-healthy mean and SD
#   SC  sex-specific centre, common scale: own mean, one SD pooled across the sexes
#   PC  pooled centre and scale: one mean and one SD for everyone
REFERENCES = ("SS", "SC", "PC")


# -- CALCULATION: coordinates under the SS, SC and PC references -------------
# Methods, "Young-healthy reference": "each transformed marker was centred and
#   scaled by the mean and SD of the same-sex young-healthy reference" (SS, the
#   primary choice).
# Results 2 compares this with two alternatives: "With sex-specific starting
#   points and pooled units" (SC) and "a fully pooled reference" (PC).
# How: coordinate = (value - centre) / scale; the centre and scale for each
#   reference are listed in the READ ME block at the top of this file.
def reference_coordinates(values: pd.DataFrame, sex: np.ndarray, ref: pd.DataFrame, kind: str) -> np.ndarray:
    """Coordinates under reference `kind`, with mean/SD taken from the reference members `ref`.

    `values` are transformed marker values (rows = participants), `sex` their
    sex codes (1 = man, 2 = woman); `ref` holds the members' transformed values
    and their `sex`.
    """
    markers = list(values.columns)
    r = ref[markers].to_numpy(float)
    v = values.to_numpy(float)
    pooled_mean, pooled_sd = r.mean(0), r.std(0, ddof=1)
    out = np.empty_like(v)
    for s in (1, 2):
        own = r[ref.sex.to_numpy() == s]
        rows = sex == s
        centre = pooled_mean if kind == "PC" else own.mean(0)
        scale = own.std(0, ddof=1) if kind == "SS" else pooled_sd
        out[rows] = (v[rows] - centre) / scale
    return out


# -- CALCULATION: adverse distance (DA) under the SS, SC and PC references ---
# The same DA rule as above (median centre, IQR/1.349 scale, adverse side
# only), rebuilt here so that Table 2 (Results 2) can recompute DA from a
# reference re-estimated without the held-out survey cycle, under each of
# SS, SC and PC. These are the "Adverse distance, DA" rows of Table 2.
def robust_adverse_distance(values: pd.DataFrame, sex: np.ndarray, ref: pd.DataFrame, kind: str,
                            orientation: np.ndarray, adverse_side: np.ndarray) -> np.ndarray:
    """DA rebuilt from a reference: the same rule the data preparation used.

    Each marker is centred on the reference median and scaled by its robust SD
    (interquartile range / 1.349), then oriented so that positive follows an
    increasing raw value. The prespecified adverse side is applied afterwards:
    "high_only" keeps positive excursions, "low_only" keeps negative ones (as
    positive numbers), and "two_sided" keeps both. DA is the mean over markers.
    """
    markers = list(values.columns)
    r = ref[markers].to_numpy(float)
    v = values.to_numpy(float)

    def median_and_scale(x):
        q1, q2, q3 = np.quantile(x, [.25, .5, .75], axis=0, method="linear")
        return q2, (q3 - q1) / 1.349

    pooled = median_and_scale(r)
    z = np.empty_like(v)
    for s in (1, 2):
        own = median_and_scale(r[ref.sex.to_numpy() == s])
        rows = sex == s
        centre = pooled[0] if kind == "PC" else own[0]
        scale = own[1] if kind == "SS" else pooled[1]
        z[rows] = orientation * (v[rows] - centre) / scale
    excursion = np.where(adverse_side == "high_only", np.maximum(z, 0),
                         np.where(adverse_side == "low_only", np.maximum(-z, 0), np.abs(z)))
    return excursion.mean(axis=1)


# -- CALCULATION: rescaling to cohort SD -------------------------------------
# Every distance and score is divided by its cohort SD before it enters a model,
# so that differences and hazard ratios read "per SD", as reported in the
# Results ("women's adverse distance exceeded men's by 0.19 standard
# deviations").
def standardise(values: np.ndarray) -> np.ndarray:
    """Rescale to mean 0 and SD 1 in the cohort, so effects read 'per SD'."""
    return (values - values.mean()) / values.std(ddof=1)
