# 8. Coordinates for the adverse distance (DA)

The adverse distance counts how far each marker has moved beyond the
young-healthy range in its adverse direction only. These scripts build the
per-marker values it is made from. The distance itself is summed in
`3_results/measures.py`.

| Step | Script | Makes |
|---|---|---|
| 1 | `01_lock_reference.py` | Fixes the reference members used here, the same young-healthy reference as folder 7 (`fi_reference_membership.csv`) |
| 2 | `02_estimate_reference_parameters.py` | Each marker's reference median and robust spread (IQR / 1.349), by sex, on the transformed scale, so a few extreme reference values cannot move them (`fi_reference_robust_parameters.csv`) |
| 3 | `03_build_robust_coordinates.py` | Every participant's markers on that scale: (value − reference median) / spread (`fi_raw_oriented_robust_coordinates.csv`) |
| 4 | `04_build_adverse_excursions.py` | Applies each marker's adverse direction and keeps only movement into the adverse tail, so every excursion is zero or positive (`fi_adverse_excursions.csv`) |
| 5 | `05_build_deficit_severities.py` | Severity of each excursion. In the final version severity equals the excursion: no dead zone, no cap (`fi_absolute_deficit_severities.csv`) |
| 6 | `06_build_signed_deficits.py` | The final per-marker values, with the sign showing which tail (`fi_signed_graded_deficits.csv`), read by the Results |

The rules (which direction is adverse for each marker) are in
`../config/fi_scoring_spec.csv`. "FI" is the original internal name; the paper
calls this the adverse distance. See `docs/adverse_distance.md`.
