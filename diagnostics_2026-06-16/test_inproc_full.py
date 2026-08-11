#!/usr/bin/env python3
"""In-process FVS with the FULL keyword setup but DSNOUT on a separate file (2026-06-18).
The trimmed keyfile segfaults (incomplete init); the full keyfile runs further but errored on the output
DB when DSNIN and DSNOUT shared a file. Here DSNOUT is a separate scratch db; results read via the API."""
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
KEY = """STDIDENT
{sid}
DATABASE
DSNIN
{indb}
DSNOUT
{outdb}
STANDSQL
SELECT * FROM fvs_standinit WHERE stand_id = '%StandID%'
ENDSQL
TREESQL
SELECT * FROM fvs_treeinit WHERE stand_id = '%StandID%'
ENDSQL
END
DATABASE
SUMMARY            2
END
TREELIST           0
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
t1,sc,al=sel; sid=str(t1)
def go(inject):
    tmp=tempfile.mkdtemp(); indb=os.path.join(tmp,"in.db"); outdb=os.path.join(tmp,"out.db")
    std=build_fvs_standinit({"INVYR":2000,"STATECD":sc,"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
    tdf=build_fvs_treeinit(al,sid)
    c=sqlite3.connect(indb); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    kp=os.path.join(tmp,"t.key"); open(kp,"w").write(KEY.format(sid=sid,indb=indb,outdb=outdb))
    fvs=FVS(LIB); fvs.load_keyfile(kp)
    k=0
    while fvs.itrncd==0 and k<8:
        fvs.run(stop_point_code=5, stop_point_year=-1); k+=1
        if fvs.restart_code==5 and inject and fvs.dims["ntrees"]>0:
            dg=fvs.get_tree_attr("dg"); fvs.set_tree_attr("dg", dg*1.30)
        if fvs.restart_code==100 or fvs.itrncd!=0: break
    s=fvs.summary; bacol=next((c for c in (s.columns if s is not None else []) if str(c).upper()=="BA"),None)
    ba=float(s[bacol].iloc[-1]) if (s is not None and len(s) and bacol) else float("nan")
    return ba, (list(s.columns) if s is not None else None)
ba_b,cols=go(False); pr("base BA:",round(ba_b,2)," cols:",cols)
ba_i,_=go(True); pr("inject(+30%% dg) BA:",round(ba_i,2)," delta:",round(ba_i-ba_b,2))
pr("IN_ENGINE_INJECTION_WORKS" if (ba_i>ba_b+0.01) else "RAN_NO_INJECT_EFFECT")
pr("DONE_FULL")
