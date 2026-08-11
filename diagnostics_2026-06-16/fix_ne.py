#!/usr/bin/env python3
"""Ordered NE component fix test (2026-06-16). Signed, region-level adjustments applied to the FVS
engine in dependency order, each measured against observed FIA remeasurement:
  default
  +allometry   : CRNMULT (crown) + REGHMULT (height) from the NE refit (foundation; small stand effect)
  +growth      : signed BAIMULT (diameter-growth slowdown) to attack QMD/BA over-prediction
  +mortality   : signed MORTMULT (<1 reduces mortality, keeps small trees) to attack TPH/QMD
  +combined    : allometry + growth + mortality + localized max SDI
Reports BA/TPH/QMD bias vs observed. Env NSAMP, SEED."""
import os,sys,math,tempfile,sqlite3,json
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from config.config_loader import FvsConfigLoader
from perseus_100yr_projection import build_fvs_treeinit, build_fvs_standinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
VAR="ne"; STATES="CT,ME,MA,NH,NY,RI,VT".split(","); MAXSP=108
NSAMP=int(os.environ.get("NSAMP","120")); SEED=int(os.environ.get("SEED","5"))
M2HA=0.2296; TPHc=2.4710538; CMc=2.54
_REV={"CT":9,"ME":23,"MA":25,"NH":33,"NY":36,"RI":44,"VT":50}
G.VARIANT_STATES[VAR]=tuple(_REV.values()); _o=G._state_abbrev; _A={v:k for k,v in _REV.items()}; G._state_abbrev=lambda c:_A.get(c) or _o(c)
ld=FvsConfigLoader(VAR,version="calibrated",config_dir=P+"/config")
cm=ld.config.get("calibration_multipliers",{}); cr=cm.get("cr_multiplier",[]); hd=cm.get("htdbh_multiplier",[])
sdi=ld._find_sdi_param(ld.config.get("categories",{})) or []
# keyword builders (verified field layouts)
def crn(): return "\n".join("CRNMULT         %10d%10.4f"%(i+1,float(cr[i])) for i in range(len(cr)) if isinstance(cr[i],(int,float)))
def regh(): return "\n".join("REGHMULT        %10d%10.4f"%(i+1,float(hd[i])) for i in range(len(hd)) if isinstance(hd[i],(int,float)))
def bai(s): return "\n".join("BAIMULT         %10d%10.4f"%(i+1,s) for i in range(1,MAXSP+1))
def mort(s): return "\n".join("MORTMULT        %10d%10.4f       0.0     999.0"%(i+1,s) for i in range(1,MAXSP+1))
def maxsdi(): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,float(v)*1.4) for i,v in enumerate(sdi) if isinstance(v,(int,float)) and v>0)
ARMS={"default":"", "+allometry":crn()+"\n"+regh(), "+growth(bai.7)":bai(0.7),
      "+mortality(.6)":mort(0.6), "+combined":crn()+"\n"+regh()+"\n"+bai(0.8)+"\n"+mort(0.6)+"\n"+maxsdi()}
def runkw(std,tdf,sid,kw,yrs):
    with tempfile.TemporaryDirectory() as tmp:
        db=os.path.join(tmp,"FVS_Data.db"); con=sqlite3.connect(db)
        std.to_sql("fvs_standinit",con,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",con,if_exists="replace",index=False); con.close()
        open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=1,cycle_length=yrs))
        return _run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVSne"),os.path.join(tmp,"t.key"),db,tmp)
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc, math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def stt(p,o):
    m=[(a,b) for a,b in zip(p,o) if a==a and b==b and b>0]; k=len(m)
    return (float("nan"),0) if k==0 else (100*sum(a-b for a,b in m)/sum(b for _,b in m), k)
def main():
    fr=[pd.read_csv(Path(FIA)/f"{ab}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for ab in STATES if (Path(FIA)/f"{ab}_PLOT.csv").exists()]
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEED)
    tr1=G.load_fia_trees(VAR,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(VAR,rem.CN.tolist(),Path(FIA))
    A={a:{x:[] for x in ["BA","TPH","QMD","oBA","oTPH","oQMD"]} for a in ARMS}; n=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); tr=tr1[tr1.PLT_CN==t1]; tro=tr2[tr2.PLT_CN==int(r.CN)]
        if len(tr)<5 or len(tro)<3: continue
        oBA,oTPH,oQMD=metr(tro)
        if oBA<=0: continue
        sid=str(t1); yrs=int(r.interval)
        std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
        tdf=build_fvs_treeinit(tr,sid); n+=1
        for a,kw in ARMS.items():
            try:
                s=runkw(std,tdf,sid,kw,yrs).get("summary")
                if s is None or len(s)==0: continue
                l=s.iloc[-1]; ba=g(l,"BA")*M2HA
                if ba<=0: continue
                A[a]["BA"].append(ba); A[a]["TPH"].append(g(l,"Tpa")*TPHc); A[a]["QMD"].append(g(l,"QMD")*CMc)
                A[a]["oBA"].append(oBA); A[a]["oTPH"].append(oTPH); A[a]["oQMD"].append(oQMD)
            except Exception: pass
    print(f"NE ordered component fix, n_run={n}")
    print(f"{'arm':<16}{'BA bias%':>10}{'TPH bias%':>11}{'QMD bias%':>11}{'n':>6}")
    for a in ARMS:
        b=stt(A[a]["BA"],A[a]["oBA"]); t=stt(A[a]["TPH"],A[a]["oTPH"]); q=stt(A[a]["QMD"],A[a]["oQMD"])
        print(f"{a:<16}{b[0]:>10.1f}{t[0]:>11.1f}{q[0]:>11.1f}{t[1]:>6d}")
    print("DONE_FIXNE")
if __name__=="__main__": main()
