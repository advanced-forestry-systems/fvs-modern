#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
PROJECT_ROOT=os.environ.get("FVS_PROJECT_ROOT",os.path.expanduser("~/fvs-modern"))
sys.path.insert(0,PROJECT_ROOT); sys.path.insert(0,os.path.join(PROJECT_ROOT,"calibration","python"))
import perseus_100yr_projection as P
sys.path.insert(0,"/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import run_conus_task_fvstreeinit as R
FIPS=R.FIPS
SI_DIR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant"
TI_DIR="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit"
variant="ne"
si=pd.read_csv(os.path.join(SI_DIR,f"standinit_{variant.upper()}.csv"),low_memory=False)
si["STAND_CN"]=si["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
checked=0; fully=0
for state_fips,grp in si.groupby("STATE"):
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
        try:
            sdf=P.build_fvs_standinit(pd_,sid,variant); tdf=R.treeinit_for_stand(fr_rows,sid)
        except Exception: continue
        if tdf.empty: continue
        fr=P.run_fvs_projection(sdf,tdf,sid,variant,config_version=None,num_cycles=20,cycle_length=5)
        tls=fr["treelists"]; yrs=sorted(tls.keys())
        checked+=1
        print(f"stand {cn}: n_trees_init={len(tdf)} n_years={len(yrs)} yrs={yrs} exit={fr['exit_code']}")
        if len(yrs)>=5:
            fully+=1
            # TreeId stability across first->last cycle
            a=tls[yrs[0]]; b=tls[yrs[-1]]
            sa=set(a["TreeId"].astype(str)); sb=set(b["TreeId"].astype(str))
            print(f"   TreeId overlap yr{yrs[0]}({len(sa)}) -> yr{yrs[-1]}({len(sb)}): {len(sa&sb)} persist")
            print(f"   yr{yrs[0]} TPA range [{a['TPA'].min():.3g},{a['TPA'].max():.3g}] sumTPA={a['TPA'].sum():.1f}")
            print(f"   yr{yrs[-1]} sumTPA={b['TPA'].sum():.1f}  cols PctCr present={'PctCr' in b.columns}")
        if checked>=8: break
    if checked>=8: break
print(f"DONE checked={checked} fully_projecting={fully}")
