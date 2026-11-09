# 5. Health conditions

| Step | Script | Makes |
|---|---|---|
| 1–2 | `01_import_self_reported_disease.py`, `02_harmonise_self_reported_disease.py` | Flags for each self-reported condition: heart disease, stroke, diabetes, cancer, lung, liver, thyroid and kidney disease, arthritis (`self_reported_disease_flags.csv`) |
| 3 | `03_import_functional_limits.py` | The three functional-limitation questions (PFQ051, PFQ054, PFQ057; not asked in 1999–2002 or below age 20) (`functional_impairment.csv`) |

These screen who can join the young-healthy reference (`7_marker_panel/`).
Where a cycle did not ask a question, the answer stays missing; it is never
counted as "no".
