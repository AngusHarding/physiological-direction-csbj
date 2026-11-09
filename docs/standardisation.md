# HD biomarker standardisation

## Method

The locked transformed panel is standardised separately by sex using the locked
availability-aware healthy reference. For marker `j` and participant sex `s`:

`z = (transformed value - healthy-sex reference mean) / healthy-sex reference SD`

The SD is the sample SD (`ddof=1`). The same two parameters are applied to every
participant of that sex, not only reference members. Transformations are not
reselected. No trimming, winsorisation, capping, offsets or missing-value
imputation occurs.

This is an HD-branch derivative. It is not an FI input and it does not yet
contain covariance matrices or Mahalanobis-distance scores.

## Outputs

- `final_38_panel_transformed_standardised_hd.csv`: identifiers, sex, age,
  reference flag and 38 standardised coordinates;
- `hd_reference_standardisation_parameters.csv`: 38 markers x 2 sexes, with
  transformation, orientation, reference mean, reference SD and reference n;
- `hd_reference_standardisation_checks.csv`: numerical mean-zero/SD-one checks.

The reference contains 547 men and 691 women. Within each reference sex, every
saved marker has mean zero and sample SD one to numerical precision.

## Certified hashes

- Standardised panel:
  `D21583D09C00B49097AEF53BE7F43287CD1C233281D34FB22B7EDB8305011B78`
- Reference parameter table:
  `324139B2DB3410A763EF09734E57B20536B1ADBC1080059E9A65049F1B1EE0EB`

Five tests lock participant identity, reference counts, parameter structure,
mean-zero/SD-one behaviour, direct formula reproduction and finite values.
