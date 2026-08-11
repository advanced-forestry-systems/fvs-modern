import os,sys,tempfile,sqlite3,shutil
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
G.VARIANT_STATES["ne"]=(23,); INV={23:"ME"}; _o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
plot=pd.read_csv(Path(FIA)/"ME_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","STATECD"),low_memory=False).dropna(subset=["PREV_PLT_CN"])
t1=int(plot.PREV_PLT_CN.astype("int64").iloc[0]); tr=G.load_fia_trees("ne",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]; i=0
while len(tr[(tr.STATUSCD==1)&(tr.DIA>0)])<10: i+=1; t1=int(plot.PREV_PLT_CN.astype("int64").iloc[i]); tr=G.load_fia_trees("ne",[t1],Path(FIA)); tr=tr[tr.PLT_CN==t1]
al=tr[(tr.STATUSCD==1)&(tr.DIA>0)].copy(); sid=str(t1)
print("load_fia_trees cols:",list(al.columns))
std=build_fvs_standinit({"INVYR":2000,"STATECD":23,"COUNTYCD":0},sid,"ne"); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
def run(tdf,tag):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=1,cycle_length=10))
    _run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVSne"),os.path.join(tmp,"t.key"),db,tmp)
    c=sqlite3.connect(db)
    cols=[d[1] for d in c.execute("PRAGMA table_info(FVS_TreeList)").fetchall()]
    tl=pd.read_sql_query("SELECT * FROM FVS_TreeList WHERE Year=(SELECT MIN(Year) FROM FVS_TreeList)",c)
    c.close(); shutil.rmtree(tmp,ignore_errors=True)
    print("=== %s FVS_TreeList cols ==="%tag,cols)
    keep=[x for x in ["TreeId","TreeIndex","SpeciesFIA","DBH","Ht","TPA"] if x in tl.columns]
    print(tl[keep].head(12).to_string(index=False) if keep else "no std cols")
    return tl
tdf=build_fvs_treeinit(al,sid)
print("treeinit ht (FIA-provided) sample:",list(tdf.ht.head(8)))
print("FIA ACTUALHT sample:",list(al.get("ACTUALHT", al.HT).head(8)) if "ACTUALHT" in al.columns else "no ACTUALHT in loader; HT=",list(al.HT.head(8)))
run(tdf,"WITH_HT")
tdf0=tdf.copy(); tdf0["ht"]=0.0
run(tdf0,"BLANK_HT_imputed")
print("DONE_PROBE")
