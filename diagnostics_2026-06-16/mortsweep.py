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
    return (0,0,0) if len(d)==0 else ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA,d.TPA_UNADJ.sum()*TPHc,math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
def run(std,tdf,sid,kw,yrs):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=1,cycle_length=yrs))
    r=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVSne"),os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def stt(p,o):
    m=[(a,b) for a,b in zip(p,o) if a==a and b==b and b>0]; k=len(m)
    return (float("nan"),0) if k==0 else (100*sum(a-b for a,b in m)/sum(b for _,b in m),k)
MAXSP=108
def allsp(kw,val,extra=""): return "\n".join("%-16s%10d%10.4f%s"%(kw,i+1,val,extra) for i in range(MAXSP))
DBH="       0.0     999.0"
ARMS=[("default",""),("MORT1.5",allsp("MORTMULT",1.5,DBH)),("MORT2.0",allsp("MORTMULT",2.0,DBH)),("MORT2.5",allsp("MORTMULT",2.5,DBH)),("MORT3.0",allsp("MORTMULT",3.0,DBH))]
A={lab:{x:[] for x in ["BA","QMD","TPH","oBA","oQMD","oTPH"]} for lab,_ in ARMS}; n=0
for _,r in rem.iterrows():
    t1=int(r.PREV_PLT_CN); a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==int(r.CN)]
    if len(a)<5 or len(b)<3: continue
    oBA,oTPH,oQMD=metr(b)
    if oBA<=0: continue
    sid=str(t1); yrs=int(r.interval); std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,"ne"); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(a,sid); n+=1
    for lab,kw in ARMS:
        s=run(std,tdf,sid,kw,yrs)
        if s is None or len(s)==0: continue
        l=s.iloc[-1]; ba=g(l,"BA")*M2HA
        if ba<=0: continue
        A[lab]["BA"].append(ba); A[lab]["QMD"].append(g(l,"QMD")*CMc); A[lab]["TPH"].append(g(l,"Tpa")*TPHc); A[lab]["oBA"].append(oBA); A[lab]["oQMD"].append(oQMD); A[lab]["oTPH"].append(oTPH)
print("n",n)
print("%-9s %8s %8s %8s %5s"%("arm","BA%","QMD%","TPH%","n"))
for lab,_ in ARMS:
    bb=stt(A[lab]["BA"],A[lab]["oBA"]); qq=stt(A[lab]["QMD"],A[lab]["oQMD"]); tt=stt(A[lab]["TPH"],A[lab]["oTPH"]); print("%-9s %8.1f %8.1f %8.1f %5d"%(lab,bb[0],qq[0],tt[0],bb[1]))
print("DONE_MORT")
