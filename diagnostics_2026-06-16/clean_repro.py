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
t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[0]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]
i=0
while len(tr)<8: i+=1; t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[i]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]
sid=str(t1)
std=build_fvs_standinit({"INVYR":2000,"STATECD":23,"COUNTYCD":0},sid,"NE"); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
tdf=build_fvs_treeinit(tr,sid)
tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
kp=os.path.join(tmp,"ne.key"); open(kp,"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** DEFAULT",num_cycles=2,cycle_length=10))
print("stand",sid,"trees_in",len(tdf),"cwd->tmp")
os.chdir(tmp)  # so fort.* scratch lands here
fvs=FVS(lib_path=P+"/lib/FVSne.so")
fvs.load_keyfile(kp)
print("keyfile set. running to stop_point 5 ...")
try:
    fvs.run(stop_point_code=5, stop_point_year=-1)
    print("dims:",fvs.dims,"exit:",getattr(fvs,'exit_code',None))
    print("NTREES=",fvs.dims.get("ntrees"))
except Exception as e:
    print("RUN ERROR:",type(e).__name__,e)
print("DONE")
