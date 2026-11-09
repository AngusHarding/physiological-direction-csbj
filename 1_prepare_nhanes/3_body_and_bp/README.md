# 3. Body size and blood pressure

| Step | Script | Makes |
|---|---|---|
| 1–3 | `01_import_anthropometry.py`, `02_harmonise_anthropometry.py`, `03_derive_anthropometry_features.py` | Height, weight and waist; then BMI and waist-to-height ratio (`anthropometry.csv`, `anthropometry_derived_features.csv`) |
| 4–5 | `04_import_blood_pressure.py`, `05_harmonise_blood_pressure.py` | All recorded blood-pressure readings and pulse; systolic and diastolic are the mean of readings 2 to 4 (`cardiovascular.csv`) |

The blood-pressure scripts were called `cardiovascular` in the original code;
the files they make still carry that name.
