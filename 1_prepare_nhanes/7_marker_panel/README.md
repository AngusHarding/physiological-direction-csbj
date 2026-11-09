# 7. The 38-marker panel and the young-healthy reference

This is where the paper's measurement is set up.

| Step | Script | Makes |
|---|---|---|
| 1 | `01_build_38_marker_panel.py` | Joins the 38 markers from folders 2 and 3 into one table, one row per eligible participant with every marker present (`final_38_biomarker_panel_raw.csv`) |
| 2 | `02_select_transformations.py` | Chooses each marker's transformation (none, log, reciprocal and so on) to bring it closer to a normal shape, and applies it (`final_38_panel_transformed_unstandardised_v2.csv`) |
| 3 | `03_build_analysis_master.py` | Adds demographics, lifestyle, conditions and mortality to the panel, without dropping anyone (`analysis_master_raw.csv`) |
| 4 | `04_build_young_healthy_reference.py` | Picks the young-healthy reference: under 30, BMI 18.5–24.9, never smoked, at most seven drinks a week, no reported chronic disease or functional limit. 1,238 people: 547 men, 691 women (`healthy_reference_members.csv`) |
| 5 | `05_standardise_markers.py` | Expresses every marker relative to the reference of the participant's own sex: the coordinates for Mahalanobis distance (DM) and the mortality score q (`final_38_panel_transformed_standardised_hd.csv`) |

Which markers go in, and from where, is set in `../config/panel_merge_spec.csv`
and `../config/final_38_panel_spec.csv`. See `docs/transformations.md` and
`docs/young_healthy_reference.md`.

**Known issue.** Step 3 reads the smoking, alcohol and activity tables from the
original analysis, not the rebuilt ones in `4_lifestyle/`. To fix before
release. Step 4 ends by comparing its result with an old file; that comparison
will be removed.
