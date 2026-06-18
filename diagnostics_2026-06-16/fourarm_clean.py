#!/usr/bin/env python3
"""Disturbance-clean fvs-conus vs default (NE), by filtering the existing per-condition predictions
(2026-06-17). Join validation_data.csv conditions to FIA COND status, filter to undisturbed, and compute
BA/TPA/QMD/volume bias for default FVS and the fvs-conus calibrated equations on the same disturbance-clean
conditions. This is arms A (default) and C (fvs-conus) of the four-arm comparison; arm B (fvs-modern
keyword-calibrated) comes from calib_final (separate sample, noted)."""
import os,sys
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"
import pandas as pd, numpy as np
from pathlib import Path
VD="/users/PUOM0008/crsfaaron/fvs-conus/output/comparisons_overstory_NEonly/intermediate/validation_data.csv"
NE_STATES=["CT","ME","MA","NH","NY","VT"]
def load_cond(states):
    cls={}
    for s in states:
        f=Path(FIA)/f"{s}_COND.csv"
        if not f.exists(): continue
        c=pd.read_csv(f,usecols=lambda x:x in ("PLT_CN","TRTCD1","TRTCD2","TRTCD3","DSTRBCD1","DSTRBCD2","DSTRBCD3"),low_memory=False)
        for cn,g2 in c.groupby("PLT_CN"):
            try: cn=int(cn)
            except: continue
            trt=g2[["TRTCD1","TRTCD2","TRTCD3"]].fillna(0).values; dst=g2[["DSTRBCD1","DSTRBCD2","DSTRBCD3"]].fillna(0).values
            cls[cn]="harvested" if (trt==10).any() else ("disturbed" if (dst>0).any() else "undisturbed")
    return cls
cls=load_cond(NE_STATES)
d=pd.read_csv(VD,low_memory=False)
print("validation_data rows:",len(d))
d["cls"]=d.PLT_CN_t2.astype("int64").map(lambda c: cls.get(int(c),"unk"))
print("class counts:",d.cls.value_counts().to_dict())
def biaspct(pred,obs):
    m=(~pred.isna())&(~obs.isna())&(obs>0);
    return float("nan") if m.sum()==0 else 100*(pred[m].sum()-obs[m].sum())/obs[m].sum()
for subset in ["undisturbed","harvested","ALL"]:
    s=d if subset=="ALL" else d[d.cls==subset]
    if len(s)==0: continue
    print("=== %s (n=%d) ==="%(subset,len(s)))
    for metric,obs,pc,pd_ in [("BA","BA_t2","BA_pred_calib","BA_pred_default"),
                              ("TPA","TPA_t2","TPA_pred_calib","TPA_pred_default"),
                              ("QMD","QMD_t2","QMD_pred_calib","QMD_pred_default"),
                              ("VOL_CFNET","VOL_CFNET_t2","VOL_CFNET_pred_calib","VOL_CFNET_pred_default")]:
        if obs in s and pc in s and pd_ in s:
            print("  %-9s default %+6.1f%%  fvs-conus %+6.1f%%"%(metric,biaspct(s[pd_],s[obs]),biaspct(s[pc],s[obs])))
print("DONE_FOURARM")
