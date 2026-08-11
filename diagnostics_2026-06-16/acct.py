import os,sys,math,tempfile,sqlite3,shutil
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
STATES="CT,ME,MA,NH,NY,RI,VT".split(","); M2HA=0.2296;TPHc=2.4710538;CMc=2.54
_REV={"CT":9,"ME":23,"MA":25,"NH":33,"NY":36,"RI":44,"VT":50}
G.VARIANT_STATES["ne"]=tuple(_REV.values()); _o=G._state_abbrev; _A={v:k for k,v in _REV.items()}; G._state_abbrev=lambda c:_A.get(c) or _o(c)
fr=[pd.read_csv(Path(FIA)/f"{ab}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for ab in STATES if (Path(FIA)/f"{ab}_PLOT.csv").exists()]
plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=120,random_state=5)
tr1=G.load_fia_trees("ne",rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees("ne",rem.CN.tolist(),Path(FIA))
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    return (0,0) if len(d)==0 else ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc)
def run(std,tdf,sid,yrs):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=1,cycle_length=yrs))
    r=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVSne"),os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
R={"ot1":[],"fv0":[],"ot2":[],"fvN":[],"ft0":[],"otp1":[]}
n=0
for _,r in rem.iterrows():
    t1=int(r.PREV_PLT_CN); a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==int(r.CN)]
    if len(a)<5 or len(b)<3: continue
    oBA1,oTP1=metr(a); oBA2,oTP2=metr(b)
    if oBA1<=0 or oBA2<=0: continue
    sid=str(t1); yrs=int(r.interval); std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,"ne"); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(a,sid)
    s=run(std,tdf,sid,yrs)
    if s is None or len(s)==0: continue
    f0=s.iloc[0]; fN=s.iloc[-1]; fba0=g(f0,"BA")*M2HA; fbaN=g(fN,"BA")*M2HA; ftp0=g(f0,"Tpa")*TPHc
    if fba0<=0: continue
    R["ot1"].append(oBA1); R["fv0"].append(fba0); R["ot2"].append(oBA2); R["fvN"].append(fbaN); R["ft0"].append(ftp0); R["otp1"].append(oTP1); n+=1
def mr(p,o):
    m=[(x,y) for x,y in zip(p,o) if y>0]; return 100*sum(x-y for x,y in m)/sum(y for _,y in m)
print("n",n)
print("FVS year-0 BA vs OBSERVED t1 BA (same input trees):  %+.1f%%"%mr(R["fv0"],R["ot1"]))
print("FVS year-0 TPH vs OBSERVED t1 TPH                  :  %+.1f%%"%mr(R["ft0"],R["otp1"]))
print("FVS final BA vs OBSERVED t2 BA                     :  %+.1f%%"%mr(R["fvN"],R["ot2"]))
print("OBSERVED t2/t1 BA growth                           :  %+.1f%%"%mr(R["ot2"],R["ot1"]))
print("FVS final/year0 BA growth                          :  %+.1f%%"%mr(R["fvN"],R["fv0"]))
print("mean obs t1 BA %.2f  FVS y0 BA %.2f  obs t2 %.2f  FVS yN %.2f"%(np.mean(R["ot1"]),np.mean(R["fv0"]),np.mean(R["ot2"]),np.mean(R["fvN"])))
print("DONE_ACCT")
