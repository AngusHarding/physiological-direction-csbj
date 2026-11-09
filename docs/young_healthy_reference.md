# Young-healthy reference reconstruction

## Locked primary definition

The reference is defined inside the fixed 29,053-person analysis master. It
requires age under 30, BMI 18.5-24.9 inclusive, no affirmative chronic disease,
never smoking, no more than seven alcoholic drinks per week, and no affirmative
functional impairment.

Missingness is handled by survey applicability:

- genuine missing responses in an applicable screen are `unknown` and exclude;
- chronic-disease questionnaire non-administration at ages 18-19 is explicitly
  not applicable and does not exclude in the absence of an affirmative flag;
- PFQ051/054/057 were absent from the 1999-2002 source files and are explicitly
  not collected rather than silently coded healthy;
- those functional questions are not applicable below age 20;
- missing BMI, smoking or alcohol always excludes reference membership.

The primary availability-aware reference contains 1,238 participants: 547 men
and 691 women. This provides 14.4 male and 18.2 female reference participants
per biomarker for the 38-marker panel.

## Selection cascade

| Step | Total | Men | Women |
|---|---:|---:|---:|
| Analysis master | 29,053 | 14,916 | 14,137 |
| Age under 30 | 6,811 | 3,606 | 3,205 |
| BMI 18.5-24.9 | 3,036 | 1,632 | 1,404 |
| Disease-free or structurally not applicable | 2,864 | 1,560 | 1,304 |
| Never smoker | 1,371 | 638 | 733 |
| Alcohol <=7 drinks/week | 1,272 | 565 | 707 |
| No functional impairment or structurally unavailable | 1,238 | 547 | 691 |

## Missingness sensitivity

A naïve complete-screen definition yields only 846 participants because it
systematically removes cycles and ages where questions were not administered.
The historical permissive definition yields 1,241 and exactly reproduces every
stored historical membership. The primary definition differs from it by only
three participants with genuine applicable disease-item missingness:

- one missing diabetes response;
- one missing thyroid-disease response;
- one missing liver-condition response.

These three are excluded rather than assumed healthy. The historical membership
was used only after both new definitions had been generated.

## Outputs

- `outputs/08_reference/healthy_reference_members.csv`: primary membership;
- `outputs/08_reference/healthy_reference_flags.csv`: all criterion states and
  sensitivity definitions;
- `outputs/08_reference/healthy_reference_cascade.csv`: cumulative counts;
- `outputs/08_reference/healthy_reference_criterion_counts.csv`: pass/fail/
  unknown counts;
- `provenance/healthy_reference_historical_comparison.csv`: post-freeze forensic
  comparison.

## Next boundary

Reference membership is now fixed. Next apply the locked transformations to
the master and separately estimate:

- HD sex-specific reference means, SDs and covariance matrices;
- sensitised-FI sex-specific reference medians and robust scales.

Neither branch may redefine reference membership or transformation choices.

The locked membership file SHA-256 is:

`9CFDA42A07F3DB275547A7014EFF7C11EBA8BD2F59DF93332B30CB218FEC15B2`
