"""Import HbA1c and all raw inputs required for RBC-folate harmonisation."""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))  # public_code: find shared helpers
from locations import ROOT  # public_code: folders and write guard (shared/locations.py)
sys.path.insert(0, str(HERE))
from nhanes_reader import raw_data_root, read_nhanes_xpt, validate_participant_key


def main() -> None:
    sources = pd.read_csv(ROOT / "config/hba1c_folate_sources.csv")
    base = raw_data_root(); hba1c_frames=[]; folate_frames=[]; summaries=[]
    for source in sources.itertuples(index=False):
        raw = read_nhanes_xpt(base / Path(source.relative_path))
        if source.component == "hba1c":
            frame = raw[["SEQN", "LBXGH"]].copy()
            frame.insert(1, "cycle", source.cycle); hba1c_frames.append(frame)
        elif source.component == "folate_br":
            hct = read_nhanes_xpt(base / Path(source.hct_relative_path))[["SEQN", "LBXHCT"]]
            frame = raw[["SEQN", source.value_variable, source.serum_variable]].merge(hct, on="SEQN", how="inner", validate="one_to_one")
            frame = frame.rename(columns={source.value_variable:"raw_rbc_folate", source.serum_variable:"raw_serum_folate"})
            frame.insert(1, "cycle", source.cycle); frame["assay"]="Bio-Rad"; folate_frames.append(frame)
        else:
            frame = raw[["SEQN", source.value_variable]].rename(columns={source.value_variable:"raw_rbc_folate"})
            frame.insert(1, "cycle", source.cycle)
            # Explicit float-NaN columns keep the cross-cycle schema stable and
            # avoid pandas having to guess a dtype for entirely empty columns.
            frame["raw_serum_folate"] = np.full(len(frame), np.nan, dtype=float)
            frame["LBXHCT"] = np.full(len(frame), np.nan, dtype=float)
            frame["assay"]="microbiologic"; folate_frames.append(frame)
        validate_participant_key(frame, f"{source.component} {source.cycle}")
        summaries.append({"component":source.component,"cycle":source.cycle,"source_rows":len(raw),"imported_rows":len(frame)})
    hba1c=pd.concat(hba1c_frames,ignore_index=True); folate=pd.concat(folate_frames,ignore_index=True)
    validate_participant_key(hba1c,"stacked HbA1c"); validate_participant_key(folate,"stacked folate inputs")
    out=ROOT/"outputs/01_imported"; out.mkdir(parents=True,exist_ok=True)
    hba1c.to_csv(out/"hba1c_raw.csv",index=False); folate.to_csv(out/"folate_raw_inputs.csv",index=False)
    pd.DataFrame(summaries).to_csv(out/"hba1c_folate_import_summary.csv",index=False)
    print(f"HbA1c/folate import: {len(hba1c):,} HbA1c; {len(folate):,} folate records")


if __name__ == "__main__": main()
