#!/usr/bin/env python3
"""Diagnostic: isolate where the in-process FVS run segfaults. Heavy step prints + faulthandler."""
import os, sys, tempfile, sqlite3, faulthandler
faulthandler.enable()
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [P+"/deployment/fvs2py", os.path.expanduser("~/overthin_work"), CONUS+"/python", P]:
    sys.path.insert(0, p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit
try:
    from fvs2py import FVS
except Exception:
    from fvs2py._base import FVS
LIB=os.environ.get("LIB", P+"/lib-test/FVSne.so"); VAR="ne"
def pr(*a): print(*a, flush=True)

INPROC_KEY = """STDIDENT
{stand_id}
DATABASE
DSNIN
{db_path}
STANDSQL
SELECT * FROM fvs_standinit WHERE stand_id = '{stand_id}'
ENDSQL
TREESQL
SELECT * FROM fvs_treeinit WHERE stand_id = '{stand_id}'
ENDSQL
END
TIMEINT            0        10
NUMCYCLE           2
PROCESS
STOP
"""
G.VARIANT_STATES[VAR]=(23,33,50)
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
t1,sc,al=sel; sid=str(t1); pr("stand",sid,"ntrees_input",len(al))
tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
std=build_fvs_standinit({"INVYR":2000,"STATECD":sc,"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
tdf=build_fvs_treeinit(al,sid)
c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
kp=os.path.join(tmp,"t.key"); open(kp,"w").write(INPROC_KEY.format(stand_id=sid,db_path=db))
pr("keyfile written; loading lib")
fvs=FVS(LIB); pr("lib ok; dims=",fvs.dims)
fvs.load_keyfile(kp); pr("keyfile loaded; itrncd=",fvs.itrncd)
pr("stepping straight to stop point 5 (let FVS fully initialize first); read/inject only when restart==5")
k=0
while fvs.itrncd==0 and k<6:
    fvs.run(stop_point_code=5, stop_point_year=-1); k+=1
    rc=fvs.restart_code; n=fvs.dims["ntrees"]
    pr("  step",k,"restart=",rc,"itrncd=",fvs.itrncd,"ntrees=",n)
    if rc==5 and n>0:
        dbh=fvs.get_tree_attr("dbh"); dg=fvs.get_tree_attr("dg")
        pr("    at stop5: meanDBH=",round(float(np.nanmean(dbh)),2)," meanDG=",round(float(np.nanmean(dg)),3))
        fvs.set_tree_attr("dg", dg*1.30); pr("    injected dg*1.30")
    if rc==100 or fvs.itrncd!=0: break
pr("summary cols:", list(fvs.summary.columns) if fvs.summary is not None else None)
if fvs.summary is not None and len(fvs.summary):
    bacol=next((c for c in fvs.summary.columns if str(c).upper()=="BA"),None)
    pr("final BA:", round(float(fvs.summary[bacol].iloc[-1]),2) if bacol else "n/a")
pr("DONE_DIAG")
