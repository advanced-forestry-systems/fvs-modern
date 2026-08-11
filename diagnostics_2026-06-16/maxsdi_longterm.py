#!/usr/bin/env python3
"""Long-term (100-yr) max-SDI leverage test (2026-06-17).
Project dense undisturbed stands 10 cycles x 10yr under (a) default and (b) revised SDIMAX = FIA p95, and
record the asymptotic SDI/BA/TPH. Tests whether FVS default lets stands accumulate past the FIA-observed
self-thinning limit (long-term over-stocking) and whether the revised max SDI caps the trajectory.
Reads per-variant FIA max SDI from ms.csv. Env VARS,NSAMP,SEED,NCYC."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108
NSAMP=int(os.environ.get("NSAMP","250")); SEEDR=int(os.environ.get("SEED","5")); NCYC=int(os.environ.get("NCYC","10"))
MSV=os.environ.get("MSV",os.path.expanduser("~/overthin_work/ms.csv"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/maxsdi_long.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"acd":["ME","NH","VT"],"sn":["AL","GA","SC","MS"],"ls":["MI","MN","WI"],
        "kt":["MT","ID"],"ci":["ID"],"ca":["CA"],"nc":["CA","OR"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","ne,acd,sn,ls,kt,ca,pn").split(",")
msd={}
try:
    ms=pd.read_csv(MSV)
    for _,r in ms.iterrows(): msd[r.variant]=float(r["SDI_p95"])
except Exception as e: print("no ms.csv yet",e)
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
def sdimax_kw(val): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,val) for i in range(MAXSP))
def run(std,tdf,sid,kw):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=NCYC,cycle_length=10))
    r=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+var),os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def sdi_eng(tpa_ph,qmd_cm):  # convert metric back: tpa_ph trees/ha, qmd_cm -> english
    tpa_ac=tpa_ph/TPHc; qmd_in=qmd_cm/CMc
    return tpa_ac*(qmd_in/10.0)**1.605
rows=[]
for var in VARS:
    states=VARMAP.get(var); fiamax=msd.get(var)
    if not states or not fiamax: print(var,"skip (no states or no FIA max)"); sys.stdout.flush(); continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    plots=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states]
    plot=pd.concat(plots,ignore_index=True)
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    cls=load_cond(states)
    rem=rem.sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA))
    defF=[]; revF=[]; defSDI=[]; revSDI=[]; nd=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN)
        if cls.get(t1)!="undisturbed": continue
        al=tr1[(tr1.PLT_CN==t1)&(tr1.STATUSCD==1)&(tr1.DIA>0)]
        if len(al)<8: continue
        sdi0=(al.TPA_UNADJ*(al.DIA/10.0)**1.605).sum()
        if sdi0 < 0.45*fiamax: continue            # only reasonably-stocked stands so they approach the limit
        sid=str(t1)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
        except Exception: continue
        sd=run(std,tdf,sid,""); sr=run(std,tdf,sid,sdimax_kw(fiamax))
        if sd is None or len(sd)==0 or sr is None or len(sr)==0: continue
        ld=sd.iloc[-1]; lr=sr.iloc[-1]
        bd=g(ld,"BA")*M2HA; br=g(lr,"BA")*M2HA
        if bd<=0 or br<=0: continue
        nd+=1
        defF.append(bd); revF.append(br)
        defSDI.append(sdi_eng(g(ld,"Tpa")*TPHc,g(ld,"QMD")*CMc)); revSDI.append(sdi_eng(g(lr,"Tpa")*TPHc,g(lr,"QMD")*CMc))
    if nd<8: print("%-4s too few dense stands (n=%d)"%(var,nd)); sys.stdout.flush(); continue
    import numpy as _np
    row={"variant":var,"n":nd,"FIA_maxSDI_p95":round(fiamax),
         "def_100yr_SDI":round(_np.mean(defSDI)),"rev_100yr_SDI":round(_np.mean(revSDI)),
         "def_100yr_BA_m2ha":round(_np.mean(defF),1),"rev_100yr_BA_m2ha":round(_np.mean(revF),1),
         "BA_reduction_pct":round(100*(_np.mean(revF)-_np.mean(defF))/_np.mean(defF),1)}
    rows.append(row)
    print("%-4s n=%-3d FIAmaxSDI %4d | 100yr default SDI %4d BA %5.1f | revised SDI %4d BA %5.1f | dBA %+5.1f%%"%(
        var,nd,fiamax,row["def_100yr_SDI"],row["def_100yr_BA_m2ha"],row["rev_100yr_SDI"],row["rev_100yr_BA_m2ha"],row["BA_reduction_pct"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_LONG")
