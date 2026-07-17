#!/usr/bin/env python3
"""Summarize the stratified all-variant CONUS stress (out_fvs/conus_*_b0.csv):
per-variant mean default vs calibrated AGB at decades + cal/def ratio."""
import glob, os, pandas as pd
OUT="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_fvs"
rows=[]
for f in sorted(glob.glob(OUT+"/conus_*_b0.csv")):
    try:
        d=pd.read_csv(f)
        if "CONFIG" not in d or "AGB_TONS_AC" not in d: continue
        v=d["VARIANT"].iloc[0] if "VARIANT" in d else os.path.basename(f).split("_")[1].upper()
        n=d["STAND_CN"].nunique() if "STAND_CN" in d else d.iloc[:,0].nunique()
        p=d[d["PROJ_YEAR"].isin([100])].groupby("CONFIG")["AGB_TONS_AC"].mean()
        if {"default","calibrated"}.issubset(p.index):
            rows.append({"variant":v,"stands":n,"def_y100":round(p["default"],1),
                         "cal_y100":round(p["calibrated"],1),
                         "cal/def":round(p["calibrated"]/p["default"],3)})
    except Exception as e:
        rows.append({"variant":os.path.basename(f),"stands":"ERR","def_y100":str(e)[:40]})
df=pd.DataFrame(rows)
print(f"variants summarized: {len(df)}")
print(df.to_string(index=False))
if "cal/def" in df: print("\nmean cal/def at yr100:", round(df["cal/def"].dropna().mean(),3))
