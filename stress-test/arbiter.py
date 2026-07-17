#!/usr/bin/env python3
import os, glob, numpy as np, pandas as pd
SCR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
TREEDIR="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h"
ENG=os.path.join(SCR,"out_conus_engine_v2")
BA_CONST=np.pi/(4.0*144.0)
VAR_STATES={"CR":["CA","CO","AZ","NV","NM","UT"],"PN":["OR","WA"],"WC":["OR","WA"]}
# Actually determine states per variant from engine outputs (only states actually present)
def eng_states(v):
    s=set()
    for f in glob.glob(os.path.join(ENG,f"conus_{v.lower()}_b*.csv")):
        try:
            d=pd.read_csv(f,usecols=["STATE"])
            s|=set(d["STATE"].dropna().unique())
        except Exception: pass
    return sorted(s)

def observed_ba_by_cn(states):
    """Observed stand BA from treeinit: sum(pi/4*DBH^2 * TREE_COUNT)/144, grouped by STAND_CN.
       TREE_COUNT is the per-acre TPA (= TPA_UNADJ). This equals the eq arm's year-0 BA exactly."""
    out={}
    for st in states:
        f=os.path.join(TREEDIR,f"{st}_FVS_TREEINIT_PLOT.csv")
        if not os.path.exists(f): continue
        d=pd.read_csv(f,usecols=["STAND_CN","TREE_COUNT","DIAMETER"],low_memory=False)
        d["STAND_CN"]=d["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
        d=d[(d["DIAMETER"]>=1.0)&(d["TREE_COUNT"]>0)].copy()
        d["ba"]=BA_CONST*(d["DIAMETER"]**2)*d["TREE_COUNT"]
        g=d.groupby("STAND_CN")["ba"].sum()
        for cn,ba in g.items(): out[cn]=out.get(cn,0.0)+ba
    return out

def engine_year0_ba(v):
    """Engine year-0 BA per STAND_CN (default config)."""
    rows=[]
    for f in glob.glob(os.path.join(ENG,f"conus_{v.lower()}_b*.csv")):
        try:
            d=pd.read_csv(f)
        except Exception: continue
        d=d[(d["PROJ_YEAR"]==0)&(d["CONFIG"]=="default")]
        rows.append(d[["STAND_CN","BA_FT2AC"]])
    if not rows: return {}
    d=pd.concat(rows,ignore_index=True)
    d["STAND_CN"]=d["STAND_CN"].astype(str)
    return dict(zip(d["STAND_CN"],d["BA_FT2AC"]))

for v in ["CR","PN","WC"]:
    sts=eng_states(v)
    obs=observed_ba_by_cn(sts)            # == eq arm year-0 BA
    eng=engine_year0_ba(v)
    shared=[cn for cn in eng if cn in obs and obs[cn]>0]
    if not shared:
        print(f"{v}: NO shared stands (states={sts}, obs={len(obs)}, eng={len(eng)})"); continue
    obs_a=np.array([obs[cn] for cn in shared])
    eng_a=np.array([eng[cn] for cn in shared])
    ratio=eng_a/np.where(obs_a>0,obs_a,np.nan)
    eng_zero=np.mean(eng_a==0)*100
    print(f"=== {v} (states={sts}) n_shared={len(shared)} ===")
    print(f"  OBSERVED/EQ year0 BA  mean={obs_a.mean():.2f} median={np.median(obs_a):.2f} ft2/ac")
    print(f"  ENGINE   year0 BA     mean={eng_a.mean():.2f} median={np.median(eng_a):.2f} ft2/ac")
    print(f"  ENGINE/OBS ratio      mean={np.nanmean(ratio):.3f} median={np.nanmedian(ratio):.3f}")
    print(f"  ENGINE BA=0 fraction  {eng_zero:.1f}%  (of shared)")
    # on nonzero engine subset
    nz=eng_a>0
    if nz.sum():
        r2=eng_a[nz]/obs_a[nz]
        print(f"  ENGINE/OBS ratio (engine>0 only) mean={r2.mean():.3f} median={np.median(r2):.3f} n={nz.sum()}")
