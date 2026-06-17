#!/usr/bin/env python3
"""Diagnose whether FVS ESTAB/NATURAL is parsed and adds trees on NE (2026-06-17).
One NE stand, 5 cycles x 10yr, strong NATURAL regen. Dump FVS .out tail + summary TPA by cycle + treelist
row counts by year, default vs +ESTAB. Reveals parse vs tally-threshold/timing issue."""
import os,sys,tempfile,sqlite3,shutil,subprocess,glob
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, FVS_LIB_DIR
VAR="ne"; FIPS={"ME":23}; G.VARIANT_STATES[VAR]=(23,)
INV={23:"ME"}; _o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
KT="""STDIDENT
{sid}
DATABASE
DSNIN
{db}
DSNOUT
{db}
STANDSQL
SELECT * FROM fvs_standinit WHERE stand_id = '%StandID%'
ENDSQL
TREESQL
SELECT * FROM fvs_treeinit WHERE stand_id = '%StandID%'
ENDSQL
END
DATABASE
SUMMARY            2
TREELIDB           2         2
END
TREELIST           0
TIMEINT            0        10
NUMCYCLE           5
{kw}
PROCESS
STOP
"""
plot=pd.read_csv(Path(FIA)/"ME_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","STATECD"),low_memory=False).dropna(subset=["PREV_PLT_CN"])
t1=int(plot.PREV_PLT_CN.astype("int64").iloc[0]); tr=G.load_fia_trees(VAR,[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]; i=0
while len(tr[(tr.STATUSCD==1)&(tr.DIA>0)])<8: i+=1; t1=int(plot.PREV_PLT_CN.astype("int64").iloc[i]); tr=G.load_fia_trees(VAR,[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]
al=tr[(tr.STATUSCD==1)&(tr.DIA>0)]; sid=str(t1)
std=build_fvs_standinit({"INVYR":2000,"STATECD":23,"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
# get a valid species num
def run(kw,tag):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    key=os.path.join(tmp,"t.key"); open(key,"w").write(KT.format(sid=sid,db=db,kw=kw or "** D"))
    r=subprocess.run([os.path.join(FVS_LIB_DIR,"FVS"+VAR),"--keywordfile="+key],cwd=tmp,capture_output=True,text=True,timeout=120)
    out=glob.glob(os.path.join(tmp,"*.out"))
    outtxt=open(out[0]).read() if out else ""
    c=sqlite3.connect(db)
    try: summ=pd.read_sql_query("SELECT Year,Tpa,BA FROM FVS_Summary2 ORDER BY Year",c)
    except Exception as e: summ=None; print("summ err",e)
    c.close()
    # detect establishment mention / errors
    estab_lines=[ln for ln in outtxt.splitlines() if ("ESTAB" in ln.upper() or "NATURAL" in ln.upper() or "REGEN" in ln.upper() or "ESTABLISHMENT" in ln.upper())]
    print("=== %s (rc=%d) ==="%(tag,r.returncode))
    if summ is not None: print(summ.to_string(index=False))
    print("estab/natural/regen lines in .out:",len(estab_lines))
    for ln in estab_lines[:6]: print("  |",ln.strip()[:90])
    err=[ln for ln in outtxt.splitlines() if "ERROR" in ln.upper() or "INVALID" in ln.upper()]
    for ln in err[:4]: print("  ERR|",ln.strip()[:90])
    shutil.rmtree(tmp,ignore_errors=True)
# probe species num
tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db"); c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
open(os.path.join(tmp,"t.key"),"w").write(KT.format(sid=sid,db=db,kw="** D")); subprocess.run([os.path.join(FVS_LIB_DIR,"FVS"+VAR),"--keywordfile="+os.path.join(tmp,"t.key")],cwd=tmp,capture_output=True,text=True)
c=sqlite3.connect(db);
try: ir=pd.read_sql_query("SELECT SpeciesNum,SpeciesFIA FROM FVS_InvReference",c); sp=int(ir.SpeciesNum.iloc[0])
except: sp=4
c.close(); shutil.rmtree(tmp,ignore_errors=True)
print("using species num",sp,"; input live trees:",len(al))
run("","DEFAULT")
kw="ESTAB             0\nNATURAL           1%10d     300.0     100.0      10.0      10.0\nEND"%sp
run(kw,"ESTAB_strong")
print("DONE_DIAG")
