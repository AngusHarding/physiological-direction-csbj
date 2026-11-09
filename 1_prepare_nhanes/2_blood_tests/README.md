# 2. Blood tests

Five subjects, each imported, then harmonised across cycles, then (for two of
them) used to derive further measures.

| Step | Script | Makes |
|---|---|---|
| 1–2 | `01_import_standard_biochemistry.py`, `02_harmonise_standard_biochemistry.py` | The standard biochemistry profile: liver enzymes, electrolytes, urea, uric acid, albumin, triglycerides and others (`standard_biochemistry.csv`) |
| 3–5 | `03_import_cbc.py`, `04_harmonise_cbc.py`, `05_derive_cbc_features.py` | Complete blood count (`whole_blood_count.csv`) and the white-cell counts derived from it (`cbc_derived_features.csv`) |
| 6–7 | `06_import_lipids.py`, `07_harmonise_lipids.py` | Total and HDL cholesterol (`cholesterol_hdl.csv`) |
| 8–9 | `08_import_hba1c_folate.py`, `09_harmonise_hba1c_folate.py` | HbA1c and red-cell folate (`glycated_haemoglobin.csv`, `folate.csv`) |
| 10–12 | `10_import_kidney.py`, `11_harmonise_kidney.py`, `12_derive_kidney_features.py` | Urine albumin-to-creatinine ratio and eGFR (`urine_albumin_creatinine.csv`, `kidney_features_eligible.csv`) |

Step 12 also needs `1_cohort/` (age and sex for eGFR) and step 2 (serum
creatinine).

Harmonising matters most for biochemistry: laboratories and analysers changed
between cycles, so step 2 applies fixed assay corrections and renames before
keeping complete rows.
