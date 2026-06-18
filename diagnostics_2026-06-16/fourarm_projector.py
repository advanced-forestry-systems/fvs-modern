#!/usr/bin/env python3
"""Projector-side arms of the four-arm comparison, disturbance-clean, with bootstrap CIs (2026-06-18).

Arms A' (projector default FVS equations) and C (fvs-conus species-free equations) live in the fvs-conus
standalone projector, which already wrote per-condition default and calibrated predictions to
validation_data.csv. This script filters those conditions to COND-undisturbed (matching the engine arms'
basis), computes BA/TPA/QMD/merch-volume bias for both arms vs observed t2, and adds percentile-bootstrap
95% CIs. Combined with fourarm_engine.py (engine arms A and B), this gives both halves of the four-arm on
the disturbance-clean basis; the synthesis doc reconciles them via within-framework deltas.
Env: VD (validation_data.csv path), STATES (comma FIA abbrevs), OUT.
"""
import os, sys, csv
import pandas as pd, numpy as np
from pathlib import Path
FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"
VD=os.environ.get("VD","/users/PUOM0008/crsfaaron/fvs-conus/output/comparisons_overstory_NEonly/intermediate/validation_data.csv")
STATES=os.environ.get("STATES","CT,ME,MA,NH,NY,VT").split(",")
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/fourarm_projector.csv"))

def cond_class(states):
    cls={}
    for s in states:
        f=Path(FIA)/f"{s}_COND.csv"
        if not f.exists(): continue
        c=pd.read_csv(f,usecols=lambda x:x in ("PLT_CN","TRTCD1","TRTCD2","TRTCD3","DSTRBCD1","DSTRBCD2","DSTRBCD3"),low_memory=False)
        for cn,g in c.groupby("PLT_CN"):
            try: cn=int(cn)
            except: continue
            trt=g[["TRTCD1","TRTCD2","TRTCD3"]].fillna(0).values; dst=g[["DSTRBCD1","DSTRBCD2","DSTRBCD3"]].fillna(0).values
            cls[cn]="harvested" if (trt==10).any() else ("disturbed" if (dst>0).any() else "undisturbed")
    return cls

def bias(pred,obs):
    m=(obs>0)&pred.notna()&obs.notna()
    return float("nan") if m.sum()<3 else 100*(pred[m]-obs[m]).sum()/obs[m].sum()
def boot_ci(pred,obs,nb=2000,seed=11):
    m=(obs>0)&pred.notna()&obs.notna(); p=pred[m].values.astype(float); o=obs[m].values.astype(float); n=len(o)
    if n<3: return (float("nan"),float("nan"))
    rng=np.random.default_rng(seed); out=[]
    for _ in range(nb):
        i=rng.integers(0,n,n); out.append(100*(p[i]-o[i]).sum()/o[i].sum())
    return (round(float(np.percentile(out,2.5)),1),round(float(np.percentile(out,97.5)),1))

df=pd.read_csv(VD,low_memory=False)
cls=cond_class(STATES)
df["cls"]=df.PLT_CN_t2.apply(lambda x: cls.get(int(x)) if pd.notna(x) else None)
und=df[df.cls=="undisturbed"].copy()
print("conditions: total %d, undisturbed %d"%(len(df),len(und)))
# metric -> (observed col, default pred col, calib/fvs-conus pred col)
METR={"BA":("BA_t2","BA_pred_default","BA_pred_calib"),
      "TPA":("TPA_t2","TPA_pred_default","TPA_pred_calib"),
      "QMD":("QMD_t2","QMD_pred_default","QMD_pred_calib"),
      "VOL":("VOL_CFNET_t2","VOL_CFNET_pred_default","VOL_CFNET_pred_calib")}
rows=[]
for var in sorted(und.VARIANT.dropna().unique()):
    d=und[und.VARIANT==var]
    if len(d)<10: continue
    row={"variant":str(var).lower(),"n":len(d)}
    for m,(oc,dc,cc) in METR.items():
        if oc in d and dc in d and cc in d:
            row["Aproj_%s"%m]=round(bias(d[dc],d[oc]),1); row["C_%s"%m]=round(bias(d[cc],d[oc]),1)
            row["Aproj_%s_ci"%m]="[%+.1f,%+.1f]"%boot_ci(d[dc],d[oc]); row["C_%s_ci"%m]="[%+.1f,%+.1f]"%boot_ci(d[cc],d[oc])
    rows.append(row)
    print("%-4s n=%-5d | proj-default->fvs-conus  BA %+.1f>%+.1f  TPA %+.1f>%+.1f  QMD %+.1f>%+.1f  VOL %+.1f>%+.1f"%(
        row["variant"],row["n"],row.get("Aproj_BA",float("nan")),row.get("C_BA",float("nan")),
        row.get("Aproj_TPA",float("nan")),row.get("C_TPA",float("nan")),
        row.get("Aproj_QMD",float("nan")),row.get("C_QMD",float("nan")),
        row.get("Aproj_VOL",float("nan")),row.get("C_VOL",float("nan"))))
if rows:
    with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_PROJECTOR")
