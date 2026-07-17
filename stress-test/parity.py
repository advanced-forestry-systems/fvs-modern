#!/usr/bin/env python3
import os,glob,numpy as np,pandas as pd
SCR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
TREEDIR="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h"
V3=os.path.join(SCR,"out_conus_engine_v3")
EQDIR=os.path.join(SCR,"conus_eq_proj","out_conus_eq")
BA_CONST=np.pi/(4.0*144.0)
def observed_ba_for_cns(states,cns):
    cns=set(cns); out={}
    for st in states:
        f=os.path.join(TREEDIR,f"{st}_FVS_TREEINIT_PLOT.csv")
        if not os.path.exists(f): continue
        d=pd.read_csv(f,usecols=["STAND_CN","TREE_COUNT","DIAMETER"],low_memory=False)
        d["STAND_CN"]=d["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
        d=d[d["STAND_CN"].isin(cns)]
        d=d[(d["DIAMETER"]>=1.0)&(d["TREE_COUNT"]>0)]
        d=d.assign(ba=BA_CONST*d["DIAMETER"]**2*d["TREE_COUNT"],
                   d2t=d["DIAMETER"]**2*d["TREE_COUNT"], tpa=d["TREE_COUNT"])
        for cn,sub in d.groupby("STAND_CN"):
            ba=sub["ba"].sum(); st_=sub["tpa"].sum(); d2t=sub["d2t"].sum()
            qmd=np.sqrt(d2t/st_) if st_>0 else 0; tph=st_*2.4710538
            o=out.get(cn,[0,0,0]); out[cn]=[o[0]+ba,o[1]+st_,d2t]  # store sums
    # recompute QMD/TPH from sums
    res={}
    for cn,(ba,tpa,d2t) in out.items():
        res[cn]={"BA":ba,"QMD":np.sqrt(d2t/tpa) if tpa>0 else 0,"TPH":tpa*2.4710538}
    return res
def v3_year0(v):
    rows=[]
    for f in glob.glob(os.path.join(V3,f"conus_{v.lower()}_b*.csv")):
        try:d=pd.read_csv(f)
        except:continue
        if d.empty: continue
        d=d[(d["PROJ_YEAR"]==0)&(d["CONFIG"]=="default")]
        rows.append(d[["STAND_CN","BA_FT2AC","QMD_IN","TPH"]])
    if not rows:return pd.DataFrame()
    d=pd.concat(rows,ignore_index=True); d["STAND_CN"]=d["STAND_CN"].astype(str)
    return d
def eq_year0(v):
    f=os.path.join(EQDIR,f"conus_eq_{v.lower()}_conus_b2_metrics.csv")
    if not os.path.exists(f):return pd.DataFrame()
    d=pd.read_csv(f); d=d[d["PROJ_YEAR"]==0]; d["STAND_CN"]=d["STAND_CN"].astype(str)
    return d[["STAND_CN","BA_FT2AC","QMD_IN","TPH"]]
VSTATES={"CR":["AZ","CA","CO","KS","NE","NM","OK","SD","TX","WY"],"PN":["OR","WA"],"WC":["OR","WA"]}
print(f"{'VAR':<4}{'n':>7}{'obsBA':>9}{'engBA':>9}{'eqBA':>9}{'eng/obs':>9}{'eng/eq':>9}{'within5%':>10}{'BA0%':>7}")
for v in ["CR","PN","WC"]:
    eng=v3_year0(v)
    if eng.empty: print(f"{v}: no v3 output yet"); continue
    cns=eng["STAND_CN"].tolist()
    obs=observed_ba_for_cns(VSTATES[v],cns)
    eq=eq_year0(v); eqm=dict(zip(eq["STAND_CN"],eq["BA_FT2AC"])) if not eq.empty else {}
    shared=[cn for cn in eng["STAND_CN"] if cn in obs and obs[cn]["BA"]>0]
    em=dict(zip(eng["STAND_CN"],eng["BA_FT2AC"]))
    o=np.array([obs[cn]["BA"] for cn in shared]); e=np.array([em[cn] for cn in shared])
    q=np.array([eqm.get(cn,np.nan) for cn in shared])
    within5=np.mean(np.abs(e-o)/o<=0.05)*100
    ba0=np.mean(np.array([em[cn] for cn in eng["STAND_CN"]])==0)*100
    eng_eq=np.nanmean(e/np.where(q>0,q,np.nan)) if np.isfinite(q).any() else np.nan
    print(f"{v:<4}{len(shared):>7}{o.mean():>9.2f}{e.mean():>9.2f}{np.nanmean(q):>9.2f}{(e/o).mean():>9.3f}{eng_eq:>9.3f}{within5:>9.1f}%{ba0:>6.1f}%")
