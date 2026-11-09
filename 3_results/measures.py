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


def adverse_distance(x: pd.DataFrame) -> np.ndarray:
    """DA: mean absolute adverse excursion across markers.

    The x coordinates are signed (a marker too high is positive, too low is
    negative, and zero on the healthy side), so the size of the excursion is
    the absolute value.
    """
    return np.abs(x.to_numpy(float)).mean(axis=1)


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


def robust_adverse_distance(values: pd.DataFrame, sex: np.ndarray, ref: pd.DataFrame, kind: str,
                            orientation: np.ndarray, adverse_side: np.ndarray) -> np.ndarray:
    """DA rebuilt from a reference: the same rule the data preparation used.

    Each marker is centred on the reference median and scaled by its robust SD
    (interquartile range / 1.349), oriented so that positive is the adverse
    direction where a marker has one, and only the adverse side counts:
    "high_only" keeps positive excursions, "low_only" negative ones (as
    positive numbers), "two_sided" both. DA is the mean over markers.
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


def standardise(values: np.ndarray) -> np.ndarray:
    """Rescale to mean 0 and SD 1 in the cohort, so effects read 'per SD'."""
    return (values - values.mean()) / values.std(ddof=1)
