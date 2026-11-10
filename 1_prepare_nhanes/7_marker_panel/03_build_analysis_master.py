# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/7_marker_panel/03_build_analysis_master.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Attaches questionnaire, lifestyle, socioeconomic and mortality
#               fields to the 29,053 complete-panel adults, one row each.
#
# METHODS       Supp. Section Zero, "Complete panel": "Questionnaire,
#               socioeconomic and mortality fields were then attached without
#               removing anyone, and each analysis applied its own further
#               requirements."
#
# PSEUDOCODE
#   start from the certified 29,053-row panel
#   left-join each listed table by SEQN; stop if anyone is added or lost
# =============================================================================
"""Build the inclusion-preserving participant master table.

Beginner's map
--------------
The certified 29,053-person raw biomarker panel is the anchor. Every additional
table is attached by SEQN with a left join, which may add missing cells but may
never remove or duplicate a participant. No reference group is defined here.
"""

from hashlib import sha256
from pathlib import Path

import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
BASE = ROOT / "outputs/05_panel/final_38_biomarker_panel_raw.csv"
SPEC = ROOT / "config/analysis_master_sources.csv"
PANEL_SPEC = ROOT / "config/final_38_panel_spec.csv"
OUT = ROOT / "outputs/07_master"
PROVENANCE = ROOT / "outputs/provenance"  # public_code: records go with the outputs
EXPECTED_BASE_SHA256 = "D9B2E795379571BFF781F2EE6503FA81B6AF5C04E59C47103519A0C25DC7F07A"


def fingerprint(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def resolve_source(row) -> Path:
    if row.authority == "reconstructed":
        return ROOT / row.relative_path
    # public_code: the three "historical_certified" tables (smoking, alcohol, physical
    # activity) are now read from this run's outputs; the rebuilt tables matched them
    # value for value (6 Nov 2026).
    raise ValueError(f"unknown authority: {row.authority}")


def main() -> None:
    if fingerprint(BASE) != EXPECTED_BASE_SHA256:
        raise AssertionError("certified raw-panel hash changed")
    master = pd.read_csv(BASE)
    base_columns = list(master.columns)
    base_seqn = master.SEQN.copy()
    if master.shape != (29053, 41) or master.SEQN.duplicated().any():
        raise AssertionError("unexpected base-panel structure")

    spec = pd.read_csv(SPEC)
    audit, manifest = [], []
    for row in spec.itertuples(index=False):
        path = resolve_source(row)
        columns = row.columns.split("|")
        source = pd.read_csv(path, usecols=["SEQN"] + columns)
        if source.SEQN.duplicated().any():
            raise AssertionError(f"{row.source_key} contains duplicate SEQN")
        matched = int(master.SEQN.isin(source.SEQN).sum())
        before_rows, before_columns = master.shape
        master = master.merge(source, on="SEQN", how="left", validate="one_to_one")
        audit.append({"source_key": row.source_key, "authority": row.authority,
                      "master_rows_before": before_rows, "source_rows": len(source),
                      "matched_master_rows": matched, "unmatched_master_rows": len(master) - matched,
                      "columns_added": len(columns), "master_rows_after": len(master),
                      "master_columns_after": master.shape[1], "missingness_policy": row.missingness_policy})
        manifest.append({"source_key": row.source_key, "authority": row.authority,
                         "path": str(path), "sha256": fingerprint(path),
                         "source_rows": len(source), "source_columns_used": len(columns)})
        if before_rows != len(master) or master.shape[1] != before_columns + len(columns):
            raise AssertionError(f"{row.source_key} changed rows or added unexpected columns")

    # Eligibility provenance and screening metadata are locks, not optional covariates.
    screening_required = spec.loc[spec.missingness_policy.eq("must_be_complete"), "columns"].str.cat(sep="|").split("|")
    eligibility_required = ["cycle", "HSQ500", "HSQ510", "HSQ520",
                            "eligible_demographic_complete", "eligible_recent_illness_matched",
                            "eligible_recent_illness_complete", "eligible_age_18_79",
                            "eligible_not_acute", "eligible_not_pregnant", "eligible_adult",
                            "primary_exclusion_reason"]
    required = screening_required + eligibility_required
    if master[required].isna().any().any() or not master.eligible_adult.all():
        raise AssertionError("shared eligibility or screening fields are incomplete")
    if not master.SEQN.equals(base_seqn) or not master[base_columns].equals(pd.read_csv(BASE)):
        raise AssertionError("master changed participant order or raw biomarker values")

    missingness = pd.DataFrame({
        "column": master.columns,
        "missing_n": [int(master[c].isna().sum()) for c in master.columns],
        "missing_pct": [float(master[c].isna().mean() * 100) for c in master.columns],
        "unique_nonmissing": [int(master[c].nunique(dropna=True)) for c in master.columns],
    })
    OUT.mkdir(parents=True, exist_ok=True)
    master.to_csv(OUT / "analysis_master_raw.csv", index=False)
    pd.DataFrame(audit).to_csv(OUT / "analysis_master_merge_audit.csv", index=False)
    missingness.to_csv(OUT / "analysis_master_missingness.csv", index=False)
    pd.DataFrame(manifest).to_csv(PROVENANCE / "analysis_master_source_manifest.csv", index=False)
    print(f"analysis master: {len(master):,} participants x {master.shape[1]} columns")
    print(f"mortality linked: {master.mortstat.notna().sum():,}; deaths: {master.mortstat.eq(1).sum():,}")
    print("reference membership: not yet defined")


if __name__ == "__main__":
    main()
