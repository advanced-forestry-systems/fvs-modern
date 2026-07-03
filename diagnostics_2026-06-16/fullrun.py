import os,sys,sqlite3,tempfile
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [P,CONUS+"/python",P+"/deployment/fvs2py"]: sys.path.insert(0,p)
import pandas as pd
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE
from fvs2py import FVS
plot=pd.read_csv(Path(FIA)/"ME_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","STATECD"),low_memory=False).dropna(subset=["PREV_PLT_CN"])
t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[0]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]; i=0
while len(tr)<8: i+=1; t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[i]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]
sid=str(t1)
std=build_fvs_standinit({"INVYR":2000,"STATECD":23,"COUNTYCD":0},sid,"NE"); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
tdf=build_fvs_treeinit(tr,sid)
tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
kp=os.path.join(tmp,"ne.key"); open(kp,"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** DEFAULT",num_cycles=2,cycle_length=10))
os.chdir(tmp)
fvs=FVS(lib_path=P+"/lib/FVSne.so"); fvs.load_keyfile(kp)
print("stand",sid,"trees_in",len(tdf))
try:
    fvs.run()  # FULL, no stop point
    sm=fvs.summary
    print("full-run summary rows:", 0 if sm is None else len(sm))
    if sm is not None and len(sm):
        print("cols:", [c for c in sm.columns][:8])
        r0=sm.iloc[0]
        print("row0 Tpa/BA:", r0.get("Tpa"), r0.get("BA"))
except Exception as e:
    print("full run error:",type(e).__name__,str(e)[:120])
# check the output DB for FVS_Summary2 / treelist (did trees get processed?)
con=sqlite3.connect(db)
for t in ["FVS_Summary2","FVS_TreeList","FVS_Cases"]:
    try:
        n=pd.read_sql_query(f"SELECT COUNT(*) c FROM {t}",con).c.iloc[0]; print(t,"rows:",n)
    except Exception as e: print(t,"missing")
con.close()
# fvsAddTrees availability
import subprocess
syms=subprocess.run(["nm","-D",P+"/lib/FVSne.so"],capture_output=True,text=True).stdout
print("fvsAddTrees symbols:", [l.split()[-1] for l in syms.splitlines() if "addtree" in l.lower()][:4])
print("DONE")
