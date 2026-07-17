#!/usr/bin/env python3
"""Probe: run ONE NE stand, dump FVS treelist columns + check for stable tree id."""
import os, sys
import numpy as np, pandas as pd
PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
import perseus_100yr_projection as P
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import run_conus_task_fvstreeinit as R

FIPS = R.FIPS
SI_DIR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant"
TI_DIR="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit"
variant="ne"
si=pd.read_csv(os.path.join(SI_DIR,f"standinit_{variant.upper()}.csv"),low_memory=False)
si["STAND_CN"]=si["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
done=False
for state_fips,grp in si.groupby("STATE"):
    if done: break
    try: state=FIPS[int(float(state_fips))]
    except Exception: continue
    tfile=os.path.join(TI_DIR,f"{state}_FVS_TREEINIT_PLOT.csv")
    if not os.path.exists(tfile): continue
    tt=pd.read_csv(tfile,low_memory=False)
    tt["STAND_CN"]=tt["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
    by_cn={k:v for k,v in tt.groupby("STAND_CN")}
    for _,stand in grp.iterrows():
        cn=stand["STAND_CN"]; fr_rows=by_cn.get(cn)
        if fr_rows is None or fr_rows.empty: continue
        sid=f"S{cn}"; iy=int(float(stand.get("INV_YEAR") or 2010))
        pd_={"INVYR":iy,"LAT":stand.get("LATITUDE"),"LON":stand.get("LONGITUDE"),
             "ELEV":stand.get("ELEVFT") or 500,"SLOPE":stand.get("SLOPE") or 10,
             "ASPECT":stand.get("ASPECT") or 180,"STDAGE":stand.get("AGE") or 50}
        sdf=P.build_fvs_standinit(pd_,sid,variant)
        tdf=R.treeinit_for_stand(fr_rows,sid)
        if tdf.empty: continue
        fr=P.run_fvs_projection(sdf,tdf,sid,variant,config_version=None,num_cycles=20,cycle_length=5)
        tls=fr["treelists"]
        print("STAND",cn,"n_treelist_years",len(tls),"exit",fr["exit_code"])
        yrs=sorted(tls.keys())
        print("YEARS",yrs[:6],"...")
        df=tls[yrs[0]]
        print("COLUMNS:",list(df.columns))
        idcols=[c for c in df.columns if c.lower() in ("treeid","treeindex","tree","ptindex","plot")]
        print("ID-LIKE COLS:",idcols)
        for c in idcols:
            print(f"  {c}: head={df[c].head(4).tolist()}  nunique={df[c].nunique()} nrows={len(df)}")
        # check id stability across two cycles
        if len(yrs)>=2 and idcols:
            c=idcols[0]
            s0=set(df[c].tolist()); s1=set(tls[yrs[1]][c].tolist())
            print(f"  overlap {c} yr{yrs[0]}->yr{yrs[1]}: {len(s0&s1)} of {len(s0)} / {len(s1)}")
        print("SAMPLE ROW0:",df.iloc[0].to_dict())
        done=True; break
