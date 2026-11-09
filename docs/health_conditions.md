# Self-reported disease reconstruction and definition audit

## Status

The historical disease variables have now been independently reconstructed from
the raw NHANES questionnaires. All 14 legacy fields—12 Medical Conditions
Questionnaire fields, diabetes, and kidney failure—agree with the certified
historical outputs with **zero participant-level value mismatches**.

This audit separates two questions:

1. Can the manuscript-era 13-condition definition be reproduced? **Yes.**
2. Was that definition the only or best disease set for the sex comparison?
   **No; several defensible conditions were omitted, but they differ materially
   in cycle coverage.**

## Reproducible pipeline

```text
Raw MCQ + DIQ + KIQ + BPQ + OSQ files
                  |
                  v
import_self_reported_disease.py
  - preserves questionnaire response codes
  - outer-merges instruments by SEQN
  - preserves non-fielded instruments as missing
                  |
                  v
harmonise_self_reported_disease.py
  - 1 = yes; 2 = no
  - refused / don't know / missing = unresolved
  - diabetes "borderline" = unresolved
  - harmonises thyroid and kidney variable renames
                  |
                  v
self_reported_disease_flags.csv
  - disease flags and condition_observed flags
  - historical 13-condition complete-case flag and count
                  |
                  +--> historical value comparison
                  +--> condition-specific >=52-year audit
```

### Code, configuration, and outputs

- Source map: [`1_prepare_nhanes/config/disease_sources.csv`](../1_prepare_nhanes/config/disease_sources.csv)
- Candidate registry: `config/disease_candidate_registry.csv` (original analysis; not included)
- Raw import: [`1_prepare_nhanes/5_health_conditions/01_import_self_reported_disease.py`](../1_prepare_nhanes/5_health_conditions/01_import_self_reported_disease.py)
- Harmonisation: [`1_prepare_nhanes/5_health_conditions/02_harmonise_self_reported_disease.py`](../1_prepare_nhanes/5_health_conditions/02_harmonise_self_reported_disease.py)
- Historical comparison: `src/10_validate/compare_self_reported_disease.py` (original analysis; not included)
- Candidate audit: `src/10_validate/audit_disease_candidates.py` (original analysis; not included)
- Raw responses: `1_prepare_nhanes/outputs/01_imported/self_reported_disease_raw.csv` (made by the scripts)
- Harmonised flags: `1_prepare_nhanes/outputs/02_harmonised/self_reported_disease_flags.csv` (made by the scripts)
- Historical primary set: `1_prepare_nhanes/outputs/02_harmonised/self_reported_disease_primary_13.csv` (made by the scripts)
- Cycle audit: `1_prepare_nhanes/outputs/02_harmonised/self_reported_disease_by_cycle.csv` (made by the scripts)
- Historical comparison: `provenance/self_reported_disease_comparison.csv` (original analysis; not included)
- Cohort audit: `outputs/10_validation/disease_candidate_audit_ge52.csv` (original analysis; not included)

## Historical primary 13

The historical set comprised arthritis, congestive heart failure, coronary heart
disease, angina, heart attack, stroke, emphysema, thyroid disease, chronic
bronchitis, liver condition, cancer, diabetes, and kidney failure. Blood
transfusion was present in the old MCQ output but correctly excluded because it
is a treatment/exposure history, not a diagnosis.

The primary-13 count is emitted only when all 13 answers are observed. Missing,
refused, unknown, borderline diabetes, and structurally unasked responses are
never counted as zero.

## Omitted-condition audit

The table uses the corrected broad Table 1 cohort restricted to age at least 52.
Percentages are among respondents to that condition. Mortality risk ratios are
unadjusted diagnostics, not causal estimates or final Supplementary Table models.

| Item | Cycle coverage | Answered n (%) | Men % | Women % | Women:men prevalence ratio | Exposed disease deaths | Unadjusted disease-death RR |
|---|---|---:|---:|---:|---:|---:|---:|
| Hypertension | All 10 | 12,810 (99.8) | 51.7 | 56.1 | 1.09 | 1,784 | 1.38 |
| Asthma | All 10 | 12,830 (99.9) | 9.9 | 14.3 | 1.45 | 356 | 1.02 |
| Osteoporosis | 8 of 10; eligibility changed | 10,064 (78.4) | 2.9 | 17.2 | 5.98 | 281 | 1.11 |
| Gout | 2007–2018 only | 8,639 (67.3) | 10.8 | 4.5 | 0.41 | 135 | 1.58 |
| COPD | 2013–2018 only | 4,333 (33.7) | 8.3 | 5.3 | 0.64 | 59 | 3.61 |

The corresponding biomarker cohort has 11,509 participants aged at least 52.
Its condition-specific counts and estimates are in the audit CSV.

## Decision matrix

| Candidate | Scientifically relevant? | Safe in an all-cycle burden count? | Recommended use |
|---|---|---|---|
| Hypertension | Yes | Yes | Add to expanded condition-specific analysis; evaluate sensitivity count |
| Asthma | Chronic, but not age-related | Technically yes | Retain for provenance only; exclude from the age-related disease analysis |
| Osteoporosis | Yes; especially important to sex comparison | No | Add as condition-specific row with explicit cycle restriction/adjustment |
| Gout | Yes | No | Add as a 2007–2018 secondary row |
| COPD | Yes, but overlaps respiratory diagnoses | No | Add as a 2013–2018 secondary row; do not double-count blindly |
| Blood transfusion | No (not a disease) | Not applicable | Keep excluded |

## Locked recommendation for Supplementary Table 1

1. Preserve the historical 13 rows as the exact-reproduction analysis.
2. Add hypertension, osteoporosis, gout, and COPD as an explicitly
   labelled expanded audit panel.
3. Use condition-specific answered samples; never code non-fielded cycles as no.
4. For osteoporosis, gout, and COPD, fit models only in fielded cycles and include
   cycle in adjusted models.
5. Do not create a pooled 18-condition count: it would mix observation windows,
   and COPD overlaps existing respiratory diagnoses.
6. Treat the historical 13-condition count as a sensitivity construct, not proof
   that the definition exhausts morbidity.

This preserves manuscript reproducibility while making the previously omitted
osteoporosis signal visible and auditable.

### Why asthma remains excluded

Asthma can begin at any age, but it is not an ageing-associated disease in the
sense required by this analysis. In the corrected broad cohort, lifetime asthma
prevalence was 16.1%, 13.4%, 11.5%, 12.4%, 12.3%, and 10.7% across the successive
18–29 through 70–79-year age bands (unweighted diagnostic summaries). The
participant-level age correlation was −0.041. Thus asthma is not literally
age-invariant, but it does not show the increasing age gradient that motivates
the age-related disease panel. Its exclusion is retained as an a priori
construct decision rather than a significance-based selection.
