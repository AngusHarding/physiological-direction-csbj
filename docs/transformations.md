# Full-cohort biomarker transformation protocol

## Decision being made

This is the prospective protocol for rebuilding the manuscript analysis. It
selects one transformation for each of the 38 biomarkers using the complete,
locked 29,053-person panel. The earlier young-healthy selector and the frozen
manuscript recipe are retained only for comparison after the new choices have
been locked. They cannot influence selection.

Transformation selection is unsupervised. Age, sex, healthy-reference status,
mortality, model fit and downstream scientific results are not selection
variables.

## Locked input

The sole primary input is
`outputs/05_panel/final_38_biomarker_panel_raw.csv`, certified as:

- 29,053 unique participants (14,916 men; 14,137 women);
- `SEQN`, sex, age, and the ordered 38-marker registry;
- no missing biomarker values;
- SHA-256
  `D9B2E795379571BFF781F2EE6503FA81B6AF5C04E59C47103519A0C25DC7F07A`.

The script must stop rather than operate on an input that fails any lock.

## Candidate family

Only the five transformations used by the manuscript-era selector are tested:

| Candidate | Formula | Valid raw domain | Orientation |
|---|---:|---|---:|
| identity | \(x\) | all finite values | +1 |
| cube root | \(\sqrt[3]{x}\) | all finite values | +1 |
| square root | \(\sqrt{x}\) | \(x\geq0\) | +1 |
| natural log | \(\log(x)\) | \(x>0\) | +1 |
| reciprocal | \(1/x\) | \(x>0\) | -1 |

No fitted transformation, trimming, participant deletion, winsorisation,
capping, offset, epsilon or manual shift is permitted. In particular, log and
reciprocal are invalid if a marker contains a genuine zero. This prevents
arbitrary constants from determining the spacing of observations near zero.
Box-Cox, Yeo-Johnson, log10, log1p and square are outside the candidate family.

## Candidate diagnostics

Every valid marker-candidate pair is assessed on all 29,053 participants. Four
complementary quantities are minimised:

1. absolute bias-corrected Fisher-Pearson skewness;
2. absolute bias-corrected Fisher excess kurtosis;
3. normal Q-Q deviation, defined as \(1-r^2\), where \(r\) is the correlation
   between the ordered transformed observations and their expected normal
   order statistics;
4. robust tail asymmetry,
   \(|\log[(q_{.99}-q_{.50})/(q_{.50}-q_{.01})]|\).

If a quantile gap is zero, tail asymmetry is infinite. Constant or non-finite
outputs are invalid. The ledger also records the raw minimum and maximum,
zeros, unique values, transformed range and monotonic orientation.

## Two-step selection rule

First apply the practical identity gate in `TRANSFORMATION_GO_NO_GO.md`. Identity
is locked when absolute skew is at most 1, absolute excess kurtosis is at most
3, and Q-Q R-squared is at least 0.98. A marker failing any criterion proceeds
to nonlinear testing.

For a marker requiring transformation, classify each valid nonlinear candidate
using those same three thresholds. If one or more pass all three, select the
least aggressive adequate candidate in this predeclared order:

`sqrt`, `cbrt`, `log`, `reciprocal`.

If none passes all three, retain candidates passing the greatest number of
criteria, then select the lowest equal-weight rank sum across absolute skew,
absolute excess kurtosis and Q-Q deviation. Resolve any remaining tie using the
same parsimony order and label the result `best_available_not_adequate`.

There are no manual overrides. A scientifically surprising or unstable winner
is flagged for review; it is not silently replaced. Any future protocol change
requires a new version and a complete rerun of all 38 markers.

## Stability analysis

Selection stability is descriptive and cannot alter the primary winner:

- repeat the complete selection separately in men and women;
- perform 100 deterministic sex-stratified bootstrap resamples of the full
  cohort using random seed 20260801;
- report each candidate's bootstrap selection frequency;
- flag a marker if men and women select different candidates, or if the primary
  winner is recovered in fewer than 70% of bootstrap samples.

The full-cohort winner remains authoritative even when flagged. Flags identify
where interpretation and sensitivity analysis are needed.

## Outputs and audit boundary

The executable stage writes:

- `outputs/06_transform/transformation_candidate_ledger_v2.csv`;
- `outputs/06_transform/transformation_selection_v2.csv`;
- `outputs/06_transform/transformation_sex_stability_v2.csv`;
- `outputs/06_transform/transformation_bootstrap_stability_v2.csv`;
- `outputs/06_transform/final_38_panel_transformed_unstandardised_v2.csv`;
- `provenance/transformation_manuscript_comparison_v2.csv`.

The transformed table preserves `SEQN`, sex, age, participant order and marker
order. It is not standardised. Standardisation is a later, separately audited
operation.

## Historical comparison

The manuscript-era code selected minimum absolute skew in a pooled
young-healthy reference and allowed positivity shifts. It reproduced 36 of the
38 frozen manuscript choices; eosinophils and HbA1c differed. That reconstruction
explains the old analysis but is not the rule above. Comparison with
`config/final_38_panel_spec.csv` occurs only after the new full-cohort choices
have been generated.
