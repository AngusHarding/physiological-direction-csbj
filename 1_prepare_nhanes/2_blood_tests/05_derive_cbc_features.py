# =============================================================================
# READ ME FIRST: 1_prepare_nhanes/2_blood_tests/05_derive_cbc_features.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Derives absolute leukocyte counts (total x percentage / 100),
#               calculated haematocrit and MCHC from the cleaned blood count.
#
# METHODS       Supp. Section Zero, "Blood cells".
# =============================================================================
"""Build transparent CBC derivatives from the cleaned measured layer.

These columns are deterministic functions, not additional independent blood
tests. They are therefore saved separately. Only ``LBXMC`` is a member of the
final 38-marker panel; the others support reconstruction and sensitivity work.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
SOURCE = ROOT / "outputs" / "02_harmonised" / "whole_blood_count.csv"
OUTPUT = ROOT / "outputs" / "03_derived"


def main() -> None:
    raw = pd.read_csv(SOURCE)
    derived = raw[["SEQN", "cycle"]].copy()
    # Absolute counts use measured total WBC and measured differential percent.
    derived["LBDLYMNO"] = raw.LBXWBCSI * raw.LBXLYPCT / 100
    derived["LBDMONO"] = raw.LBXWBCSI * raw.LBXMOPCT / 100
    derived["LBDNENO"] = raw.LBXWBCSI * raw.LBXNEPCT / 100
    derived["LBDEONO"] = raw.LBXWBCSI * raw.LBXEOPCT / 100
    derived["LBDBANO"] = raw.LBXWBCSI * raw.LBXBAPCT / 100
    unreliable = derived.LBDLYMNO.gt(10)
    derived.loc[unreliable, "LBDLYMNO"] = np.nan

    # These reproduce analyser-calculated red-cell indices from fundamental
    # measured inputs, preventing duplicate information in the clean layer.
    derived["LBXHCT"] = raw.LBXRBCSI * raw.LBXMCVSI / 10
    derived["LBXMCHSI"] = raw.LBXHGB * 10 / raw.LBXRBCSI
    derived["LBXMC"] = raw.LBXHGB * 100 / derived.LBXHCT
    derived["NLR"] = derived.LBDNENO / derived.LBDLYMNO
    derived["PLR"] = raw.LBXPLTSI / derived.LBDLYMNO
    derived["MLR"] = derived.LBDMONO / derived.LBDLYMNO
    derived["SII"] = raw.LBXPLTSI * derived.LBDNENO / derived.LBDLYMNO
    derived["SIRI"] = derived.LBDNENO * derived.LBDMONO / derived.LBDLYMNO
    derived["PCT"] = raw.LBXPLTSI * raw.LBXMPSI / 10000
    derived = derived.replace([np.inf, -np.inf], np.nan)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    derived.to_csv(OUTPUT / "cbc_derived_features.csv", index=False)
    pd.DataFrame([{"rule": "LBDLYMNO > 10 -> missing",
                   "values_to_missing": int(unreliable.sum())}]).to_csv(
                       OUTPUT / "cbc_derivation_qc_log.csv", index=False)
    print(f"CBC derived features: {len(derived):,} rows x {derived.shape[1]} columns")


if __name__ == "__main__":
    main()
