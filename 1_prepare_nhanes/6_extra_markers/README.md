# 6. Extra markers: CRP, fasting glucose, lung function

| Step | Script | Makes |
|---|---|---|
| 1 | `01_import_crp_glucose_spirometry.py` | CRP, fasting glucose and spirometry (including FEV1) for each cycle that measured them (`candidate_crp_raw.csv`, `candidate_glucose_raw.csv`, `candidate_spirometry_raw.csv`) |

These are not in the 38-marker panel. Results 2 and 4 use them in the
Phenotypic Age comparison; Results 7 tests whether lung function adds
prediction beyond the panel.

This script came from the original analysis stage, not from data preparation,
so its file paths follow a different pattern from the other folders.
