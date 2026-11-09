"""Screen HbA1c and place all RBC-folate cycles on the MA assay scale."""

from pathlib import Path
import numpy as np
import pandas as pd

import sys  # public_code: needed for the line below
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)


def transform_br_to_ma(rbf, fol, hct):
    """Frozen four-step Fazili/NHANES Bio-Rad to microbiologic conversion."""
    with np.errstate(invalid="ignore",divide="ignore"):
        h=hct/100; wbf=rbf*h+fol*(1-h)
        wbf_adj=10**(0.2204+1.017*np.log10(wbf)); x=np.log10(fol)
        fol_adj=10**(0.0188*x**3-2.7109*x**(-0.5)+3.8276)
        return (wbf_adj-fol_adj*(1-h))/h


def main() -> None:
    inp=ROOT/"outputs/01_imported"; out=ROOT/"outputs/02_harmonised"; out.mkdir(parents=True,exist_ok=True)
    ghb=pd.read_csv(inp/"hba1c_raw.csv"); invalid=ghb.LBXGH.notna() & ~ghb.LBXGH.between(3,20)
    ghb.loc[invalid,"LBXGH"]=np.nan; ghb=ghb.dropna(subset=["LBXGH"])[["SEQN","cycle","LBXGH"]].reset_index(drop=True)
    fol=pd.read_csv(inp/"folate_raw_inputs.csv"); br=fol.assay.eq("Bio-Rad")
    fol["RBC_folate"]=fol.raw_rbc_folate
    fol.loc[br,"RBC_folate"]=transform_br_to_ma(fol.loc[br,"raw_rbc_folate"],fol.loc[br,"raw_serum_folate"],fol.loc[br,"LBXHCT"])
    nonpositive=fol.RBC_folate.le(0); fol.loc[nonpositive,"RBC_folate"]=np.nan
    fol=fol.replace([np.inf,-np.inf],np.nan)
    fol_clean=fol.dropna(subset=["RBC_folate"])[["SEQN","cycle","RBC_folate"]].reset_index(drop=True)
    ghb.to_csv(out/"glycated_haemoglobin.csv",index=False); fol_clean.to_csv(out/"folate.csv",index=False)
    pd.DataFrame([{"rule":"HbA1c outside 3-20 -> missing","values_to_missing":int(invalid.sum())},
                  {"rule":"RBC folate non-positive -> missing","values_to_missing":int(nonpositive.sum())}]).to_csv(out/"hba1c_folate_change_log.csv",index=False)
    print(f"HbA1c/folate harmonisation: {len(ghb):,} HbA1c; {len(fol_clean):,} RBC folate")


if __name__ == "__main__": main()
