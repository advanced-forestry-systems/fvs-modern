#!/usr/bin/env python3
"""Validate fvs2py in-process tree-attribute injection (Route A unblock test, 2026-06-18).

Proves the new get_tree_attr/set_tree_attr binding enables true in-engine growth override:
load one stand in process, step cycles stopping at restart code 5 (growth/mortality computed but not
applied), read per-tree dg, multiply it (stand-in for an fvs-conus prediction), write it back, resume.
Compare final basal area against an un-injected run. A change proves the write reaches the engine."""
import os, sys, tempfile, sqlite3, shutil
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [P+"/deployment/fvs2py", os.path.expanduser("~/overthin_work"), CONUS+"/python", P]:
    sys.path.insert(0, p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE
try:
    from fvs2py import FVS
except Exception:
    from fvs2py._base import FVS
LIB=os.environ.get("LIB", P+"/lib/FVSne.so"); VAR="ne"

# build one NE remeasurement stand
G.VARIANT_STATES[VAR]=(23,33,50)  # ME,NH,VT
plot=pd.read_csv(Path(FIA)/"ME_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False)
yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=40,random_state=5)
tr1=G.load_fia_trees(VAR,rem.PREV_PLT_CN.tolist(),Path(FIA))
sel=None
for _,r in rem.iterrows():
    a=tr1[tr1.PLT_CN==int(r.PREV_PLT_CN)]; al=a[(a.STATUSCD==1)&(a.DIA>0)]
    if len(al)>=10: sel=(int(r.PREV_PLT_CN),int(r.STATECD),al); break
assert sel, "no suitable stand"
t1,sc,al=sel; sid=str(t1)

def make_key(tag):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    std=build_fvs_standinit({"INVYR":2000,"STATECD":sc,"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
    tdf=build_fvs_treeinit(al,sid)
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    kp=os.path.join(tmp,"t.key")
    open(kp,"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=2,cycle_length=10))
    return tmp,kp

def final_ba(inject):
    tmp,kp=make_key("inj" if inject else "base")
    fvs=FVS(LIB); fvs.load_keyfile(kp)
    fvs.run(stop_point_code=7)            # read input
    n0=fvs.dims["ntrees"]
    dbh0=fvs.get_tree_attr("dbh") if n0>0 else np.array([])
    steps=0
    while fvs.itrncd==0:
        fvs.run(stop_point_code=5, stop_point_year=-1)   # stop after growth computed, before applied
        if fvs.restart_code==5:
            if inject:
                dg=fvs.get_tree_attr("dg")
                fvs.set_tree_attr("dg", dg*1.30)         # +30% diameter growth (stand-in prediction)
            steps+=1
        if fvs.restart_code==100 or fvs.itrncd!=0:
            break
    while fvs.itrncd==0:            # finalize so the summary is fully written
        fvs.run(stop_point_code=0, stop_point_year=0)
        if fvs.restart_code==100 or fvs.itrncd!=0:
            break
    summ=fvs.summary                # property, not a method
    shutil.rmtree(tmp,ignore_errors=True)
    bacol=next((c for c in (summ.columns if summ is not None else []) if str(c).upper()=="BA"),None)
    ba=float(summ[bacol].iloc[-1]) if (summ is not None and len(summ) and bacol) else float("nan")
    print("  [%s] ntrees=%d meanDBH=%.2f stops=%d cols=%s BA=%.2f"%(
        "inject" if inject else "base", n0, float(dbh0.mean()) if len(dbh0) else float("nan"),
        steps, list(summ.columns) if summ is not None else None, ba))
    return n0, float(dbh0.mean()) if len(dbh0) else float("nan"), steps, ba

n0,dbhm,st_b,ba_base=final_ba(False)
_,_,st_i,ba_inj=final_ba(True)
print("ntrees read in-process:", n0, " mean dbh:", round(dbhm,2))
print("stop-point-5 stops:", st_b)
print("final BA  base:", round(ba_base,2), " injected(+30%% dg):", round(ba_inj,2), " delta:", round(ba_inj-ba_base,2))
print("INJECTION_WORKS" if (ba_inj>ba_base+0.01) else "NO_EFFECT_CHECK")
print("DONE_TREEATTR_TEST")
