import os,sys,sqlite3,tempfile,ctypes as ct
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [P,CONUS+"/python"]: sys.path.insert(0,p)
import pandas as pd
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE
# build one NE stand
plot=pd.read_csv(Path(FIA)/"ME_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","STATECD"),low_memory=False).dropna(subset=["PREV_PLT_CN"])
t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[0]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]; i=0
while len(tr)<8: i+=1; t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[i]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]
sid=str(t1); std=build_fvs_standinit({"INVYR":2000,"STATECD":23,"COUNTYCD":0},sid,"NE"); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
tdf=build_fvs_treeinit(tr,sid)
tmp=os.path.expanduser("~/overthin_work/ctx"); os.makedirs(tmp,exist_ok=True); db=os.path.join(tmp,"FVS_Data.db")
c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
kp=os.path.join(tmp,"ne.key"); open(kp,"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** DEFAULT",num_cycles=2,cycle_length=10))
os.chdir(tmp)
print("stand",sid,"trees_in",len(tdf))
# write 4 filenames to a stdin file, redirect fd 0 (replicate subprocess)
sin=os.path.join(tmp,"stdin.txt"); open(sin,"w").write("ne.key\nnone.tre\nne.out\nne.trl\n")
fin=os.open(sin,os.O_RDONLY); os.dup2(fin,0)
# load lib raw, blank cmdline, then fvs_ loop
lib=ct.CDLL(P+"/lib/FVSne.so", mode=ct.RTLD_GLOBAL)
setcl=lib.fvssetcmdline_; setcl.restype=None
blank=b" "
setcl(ct.c_char_p(blank), ct.byref(ct.c_int(1)), ct.byref(ct.c_int(0)), ct.c_int(1))
fvs=lib.fvs_; fvs.restype=None
rc=ct.c_int(0)
for _ in range(200):
    fvs(ct.byref(rc))
    # check restart/return code
    if rc.value!=0: break
print("fvs_ loop done, rc=",rc.value)
con=sqlite3.connect(db)
for t in ["FVS_Summary2","FVS_TreeList"]:
    try: print(t,"rows:",pd.read_sql_query("SELECT COUNT(*) c FROM "+t,con).c.iloc[0])
    except Exception as e: print(t,"missing")
try:
    s=pd.read_sql_query("SELECT Year,Tpa,BA FROM FVS_Summary2",con)
    print(s.to_string(index=False))
except: pass
con.close()
print("CTYPES_DONE")
