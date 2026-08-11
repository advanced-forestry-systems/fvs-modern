#!/usr/bin/env python3
"""FIA-observed maximum Stand Density Index per variant (2026-06-17).
Reineke summation SDI per undisturbed FIA plot: SDI = sum_i TPA_i*(DBH_i/10)^1.605 (English, DBH inches,
TPA trees/acre). High percentiles per variant approximate the empirical self-thinning limit = revised max
SDI for SDIMAX. Reports English and metric (x2.471) percentiles. Env VARS,NSAMP,SEED."""
import os,sys,math,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P)
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
NSAMP=int(os.environ.get("NSAMP","2500")); SEEDR=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/maxsdi_fia.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"acd":["ME","NH","VT"],"sn":["AL","GA","SC","MS"],"ls":["MI","MN","WI"],
        "cs":["IL","IN","MO"],"ie":["ID","MT"],"kt":["MT","ID"],"ci":["ID"],"cr":["CO","WY"],"ut":["UT","NV"],
        "ca":["CA"],"nc":["CA","OR"],"ec":["OR","WA"],"wc":["OR","WA"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","").split(",") if os.environ.get("VARS") else list(VARMAP)
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
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
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    plots=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","MEASYEAR","STATECD"),low_memory=False) for s in states if (Path(FIA)/f"{s}_PLOT.csv").exists()]
    plot=pd.concat(plots,ignore_index=True); plot["CN"]=plot.CN.astype("int64")
    cls=load_cond(states)
    # use most recent measurement of each plot; sample plot CNs
    cns=plot.CN.drop_duplicates()
    if len(cns)>NSAMP: cns=cns.sample(n=NSAMP,random_state=SEEDR)
    tr=G.load_fia_trees(var,cns.tolist(),Path(FIA))
    sdis=[]
    for cn,grp in tr.groupby("PLT_CN"):
        if cls.get(int(cn))!="undisturbed": continue
        d=grp[(grp.STATUSCD==1)&(grp.DIA>0)&(grp.TPA_UNADJ>0)]
        if len(d)<3: continue
        sdi=(d.TPA_UNADJ*(d.DIA/10.0)**1.605).sum()   # English summation SDI (per acre)
        ba=(d.TPA_UNADJ*0.005454*d.DIA**2).sum()
        if sdi>0 and ba>20: sdis.append(sdi)            # require some stocking
    if len(sdis)<30: print(var,"too few plots (n=%d)"%len(sdis)); sys.stdout.flush(); continue
    s=np.array(sdis)
    p50,p90,p95,p99,pmax=np.percentile(s,[50,90,95,99]).tolist()+[s.max()]
    row={"variant":var,"n_plots":len(s),"SDI_p50":round(p50),"SDI_p90":round(p90),"SDI_p95":round(p95),"SDI_p99":round(p99),"SDI_max":round(pmax),
         "maxSDI_metric_p95":round(p95*2.471)}
    rows.append(row)
    print("%-4s n=%-5d  p50 %4d  p90 %4d  p95 %4d  p99 %4d  max %4d   (English SDI/ac)"%(var,len(s),p50,p90,p95,p99,pmax)); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_MAXSDI")
