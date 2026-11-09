# Preparing NHANES III (1988–1994)

NHANES III is the second survey era. It is prepared the same way as the modern
data, so the Results can run every analysis in both eras.

| Step | Script | Makes |
|---|---|---|
| 1 | `01_build_panel.py` | Reads the raw adult, examination and laboratory files and builds the marker panel, mapped to the modern marker names (`nhanes3_panel_raw.csv`, `nhanes3_panel_markers.csv`) |
| 2 | `02_select_transformations.py` | Each marker's transformation (`transformation_selected.csv`) |
| 3 | `03_build_young_healthy_reference.py` | The NHANES III young-healthy reference (`reference_members.csv`) |
| 4 | `04_build_coordinates.py` | Coordinates for the adverse distance and for the Mahalanobis distance, and the reference parameters (`coordinates_fi.csv`, `coordinates_hd.csv`, `reference_parameters.csv`) |
| 5 | `05_extract_results_variables.py` | Copies the other NHANES III variables the Results use (education, income, smoking, activity, conditions, BMI, alcohol, blood pressure, spirometry, bone density, bioimpedance, balance) and the linked mortality file into tables (`results_adult_variables.csv`, `results_exam_variables.csv`, `mortality_linked.csv`) |

**Run it** from this folder, in order, after `1_prepare_nhanes/`:

```
python -B 01_build_panel.py
python -B 02_select_transformations.py
python -B 03_build_young_healthy_reference.py
python -B 04_build_coordinates.py
python -B 05_extract_results_variables.py
```

or all five at once with `python -B run_all.py`.

`nhanes3_locations.py` names the folders (the raw NHANES III folder is set in
`../data_locations.py`) and switches on the same write guard as
`1_prepare_nhanes/`: scripts write only in `outputs/`.

`nhanes3_io.py` reads the NHANES III fixed-width `.dat` files using their
`.sas` layout files. `config/nhanes3_marker_map.csv` says which NHANES III
variable matches each modern marker. `endpoints.py` defines the mortality
endpoints; `nhanes3_io.py` loads it, so every step needs it.

**Links to the modern data.** Steps 2 and 4 read two tables from
`1_prepare_nhanes/config/`: the modern transformation choices and the
adverse-distance scoring rules. Step 2 also reuses the modern transformation
code (`1_prepare_nhanes/7_marker_panel/02_select_transformations.py`), so both
eras choose transformations by the same rules. Run `1_prepare_nhanes/` first.

After step 5 the Results read only prepared tables, never a raw NHANES III
file. See `docs/nhanes_iii.md`.
