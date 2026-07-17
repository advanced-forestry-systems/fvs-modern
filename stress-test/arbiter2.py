#!/usr/bin/env python3
import os, glob, numpy as np, pandas as pd
SCR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
TREEDIR="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h"
EQDIR=os.path.join(SCR,"conus_eq_proj","out_conus_eq")
ENG=os.path.join(SCR,"out_conus_engine_v2")
BA_CONST=np.pi/(4.0*144.0)
def eng_states(v):
    s=set()
    for f in glob.glob(os.path.join(ENG,f"conus_{v.lower()}_b*.csv")):
        try: s|=set(pd.read_csv(f,usecols=["STATE"])["STATE"].dropna().unique())
        except Exception: pass
    return sorted(s)
def observed_ba(states):
    out={}
    for st in states:
        f=os.path.join(TREEDIR,f"{st}_FVS_TREEINIT_PLOT.csv")
        if not os.path.exists(f): continue
        d=pd.read_csv(f,usecols=["STAND_CN","TREE_COUNT","DIAMETER"],low_memory=False)
        d["STAND_CN"]=d["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
        d=d[(d["DIAMETER"]>=1.0)&(d["TREE_COUNT"]>0)]
        d=d.assign(ba=BA_CONST*d["DIAMETER"]**2*d["TREE_COUNT"])
        for cn,ba in d.groupby("STAND_CN")["ba"].sum().items(): out[cn]=out.get(cn,0.0)+ba
    return out
def eq_year0(v):
    # eq metrics CSV: conus_eq_<v>_conus_b2_metrics.csv
    f=os.path.join(EQDIR,f"conus_eq_{v.lower()}_conus_b2_metrics.csv")
    if not os.path.exists(f): return {}
    d=pd.read_csv(f)
    d=d[d["PROJ_YEAR"]==0]
    d["STAND_CN"]=d["STAND_CN"].astype(str)
    return dict(zip(d["STAND_CN"],d["BA_FT2AC"]))
for v in ["CR","PN","WC"]:
    sts=eng_states(v); obs=observed_ba(sts); eq=eq_year0(v)
    shared=[cn for cn in eq if cn in obs and obs[cn]>0]
    if not shared:
        print(f"{v}: eq metrics rows={len(eq)} obs={len(obs)} NO shared"); continue
    o=np.array([obs[cn] for cn in shared]); e=np.array([eq[cn] for cn in shared])
    diff=np.abs(e-o)/np.where(o>0,o,np.nan)
    print(f"{v}: n_shared(eq vs obs)={len(shared)}  obs_mean={o.mean():.2f}  eq_mean={e.mean():.2f}  median|rel diff|={np.nanmedian(diff)*100:.3f}%  frac within1%={np.mean(diff<0.01)*100:.1f}%")
