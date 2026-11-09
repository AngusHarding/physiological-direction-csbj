# Shared analysis master cohort

## Purpose

`analysis_master_raw.csv` is the participant-level spine shared by HD and FI.
It is anchored on the certified 29,053-person complete 38-biomarker panel and
adds demographic eligibility provenance, screening fields, reference-screening
covariates and mortality without deleting participants.

## Ordering

1. Demographic, age, recent-illness and pregnancy eligibility produced 34,785
   adults.
2. Complete final-panel availability produced 29,053 participants.
3. Additional fields were left-joined to those 29,053 by unique `SEQN`.
4. Missing optional covariates remain missing; they do not define membership.
5. Healthy-reference membership is deliberately not inherited from historical
   outputs and will be reconstructed inside this master next.

## Source status

The biomarker panel, eligibility provenance, screening metadata, functional
items and mortality are public reconstruction outputs. Smoking, alcohol,
activity and chronic-disease flags remain certified historical dependencies
pending their own public source-level reconstruction. Exact hashes and coverage are recorded in
`provenance/analysis_master_source_manifest.csv`.

## Branch point

After reference reconstruction:

- HD will use sex-specific healthy means, SDs and covariance;
- sensitised FI will use sex-specific healthy medians and robust scales plus
  marker directionality, dead zone and bounded severity;
- both branches will use the same locked transformations and participant spine.

No nonlinear transformation, standardisation, FI grading or reference flag is
stored in this raw master.

## Certified result

The master contains 29,053 participants and 89 columns. All 29,053 have complete
shared eligibility decisions and screening metadata. Optional coverage is
reported rather than enforced: smoking is unmatched for 1,709, alcohol for
1,563, physical activity for 2,150, diabetes for 545, kidney-failure history
for 2,232, and mortality for 30 participants.

The file SHA-256 is:

`F84DF97D19A6FAE316D5482022F97AC9CBB789A7487353353C9DD1A58E99C692`

Five tests lock participant membership, raw biomarker identity, eligibility,
optional missingness preservation and the absence of a silently inherited
reference flag.

## Work checkpoint — reference complete

The young-healthy reference has now been independently reconstructed and locked
at 1,238 participants: 547 men and 691 women. Survey non-administration is
distinguished from genuine applicable missingness. Full logic and sensitivity
results are in `HEALTHY_REFERENCE.md`.

The pipeline can now branch:

- **HD:** apply the already locked transformations, then calculate sex-specific
  healthy-reference means, SDs and covariance matrices;
- **sensitised FI:** apply the same transformations, then calculate sex-specific
  healthy-reference medians and robust scales before directional graded-deficit
  construction.

Do not rerun transformation selection. Transformations were selected on the
complete shared cohort by design.
