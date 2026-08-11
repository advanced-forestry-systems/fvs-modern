import os,sys,sqlite3,tempfile,shutil
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR
import subprocess
MAXSP=108
def allsp(kw,val,extra=""): return "\n".join("%-16s%10d%10.4f%s"%(kw,i+1,val,extra) for i in range(MAXSP))
plot=pd.read_csv(Path(FIA)/"ME_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","STATECD"),low_memory=False).dropna(subset=["PREV_PLT_CN"])
t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[0]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]; i=0
while len(tr)<10: i+=1; t1=int(plot["PREV_PLT_CN"].astype("int64").iloc[i]); tr=G.load_fia_trees("NE",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]
sid=str(t1); std=build_fvs_standinit({"INVYR":2000,"STATECD":23,"COUNTYCD":0},sid,"NE"); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
tdf=build_fvs_treeinit(tr,sid)
arms={
 "default":"",
 "NOCALIB":"NOCALIB",
 "BAIMULT0.0":allsp("BAIMULT",0.0),
 "BAIMULT2.0":allsp("BAIMULT",2.0),
 "BAIMULT0.5+NOCALIB":"NOCALIB\n"+allsp("BAIMULT",0.5),
 "CRNMULT0.2":allsp("CRNMULT",0.2),
 "CRNMULT3.0":allsp("CRNMULT",3.0),
 "REGHMULT0.5":allsp("REGHMULT",0.5),
 "REGHMULT2.0":allsp("REGHMULT",2.0),
 "MORTMULT5.0":allsp("MORTMULT",5.0,"       0.0     999.0"),
 "MORTMULT0.0":allsp("MORTMULT",0.0,"       0.0     999.0"),
}
print("stand",sid,"trees_in",len(tdf),"(3 cycles x 10yr)")
print("%-20s %18s %14s" % ("arm","BA(yr0,1,2,3)","QMD_final"))
for lab,kw in arms.items():
    d=tempfile.mkdtemp(); db=os.path.join(d,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(d,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=3,cycle_length=10))
    subprocess.run([os.path.join(FVS_LIB_DIR,"FVSne")],input="t.key\nnone.tre\no.out\no.trl\n",capture_output=True,text=True,cwd=d,timeout=120)
    c=sqlite3.connect(db)
    try: s=pd.read_sql_query("SELECT Year,BA,Tpa,QMD FROM FVS_Summary2",c)
    except: s=None
    c.close(); shutil.rmtree(d,ignore_errors=True)
    if s is None or len(s)==0: print("%-20s  NO OUTPUT"%lab); continue
    ba=",".join("%.0f"%x for x in s.BA); print("%-20s %18s %14.2f  TPHf=%.0f"%(lab, ba, s.QMD.iloc[-1], s.Tpa.iloc[-1]))
