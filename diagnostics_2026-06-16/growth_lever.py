#!/usr/bin/env python3
"""Growth-lever test (2026-06-16): does a systematic diameter-growth slowdown fix the
BA over-prediction? Emit a global BAIMULT (basal-area-increment multiplier) at scales
< 1 for all species on default FVS, measure BA/QMD/TPH bias vs observed. If BA bias
drops toward 0 as the scale drops, the BA over-prediction is growth-driven and a signed
growth recalibration is the lever. Env: VAR, STATES, MAXSP, NSAMP, SEED."""
import os, sys, math, tempfile, sqlite3
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
WORK=os.path.expanduser("~/overthin_work")
for k,v in {"FIA_DATA_DIR":FIA,"FVS_PROJECT_ROOT":P,"FVS_LIB_DIR":P+"/lib","FVS_CONFIG_DIR":P+"/config"}.items(): os.environ[k]=v
for p in [WORK,CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_treeinit, build_fvs_standinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
VAR=os.environ.get("VAR","ne"); STATES=os.environ.get("STATES","CT,ME,MA,NH,NY,RI,VT").split(",")
MAXSP=int(os.environ.get("MAXSP","108")); NSAMP=int(os.environ.get("NSAMP","120")); SEED=int(os.environ.get("SEED","5"))
M2HA=0.2296; TPHc=2.4710538; CMc=2.54
_REV={"AL":1,"CA":6,"CO":8,"CT":9,"FL":12,"GA":13,"ID":16,"IL":17,"IN":18,"IA":19,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NM":35,"NY":36,"OR":41,"RI":44,"SC":45,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
G.VARIANT_STATES[VAR]=tuple(_REV[x] for x in STATES if x in _REV)
_o=G._state_abbrev; _A={v:k for k,v in _REV.items()}; G._state_abbrev=lambda c:_A.get(c) or _o(c)
SCALES=[0.55,0.65,0.75,0.85,1.0]
def baimult(s): return "\n".join("BAIMULT         %10d%10.4f"%(i,s) for i in range(1,MAXSP+1))
def run_kw(std,tdf,sid,kw,yrs):
    with tempfile.TemporaryDirectory() as tmp:
        db=os.path.join(tmp,"FVS_Data.db"); con=sqlite3.connect(db)
        std.to_sql("fvs_standinit",con,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",con,if_exists="replace",index=False); con.close()
        open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** DEFAULT",num_cycles=1,cycle_length=yrs))
        return _run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+VAR),os.path.join(tmp,"t.key"),db,tmp)
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc, math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def stt(p,o):
    m=[(a,b) for a,b in zip(p,o) if a==a and b==b]; k=len(m)
    if k==0: return (float("nan"),float("nan"),0)
    e=[a-b for a,b in m]; mo=sum(b for _,b in m)/k
    return (100*math.sqrt(sum(x*x for x in e)/k)/mo, 100*sum(e)/k/mo, k)
def main():
    print(f"VAR={VAR} MAXSP={MAXSP} NSAMP={NSAMP}")
    frames=[pd.read_csv(Path(FIA)/f"{ab}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for ab in STATES if (Path(FIA)/f"{ab}_PLOT.csv").exists()]
    plot=pd.concat(frames,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEED)
    tr1=G.load_fia_trees(VAR,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(VAR,rem.CN.tolist(),Path(FIA))
    arms=["default"]+["bai_%.2f"%s for s in SCALES if s<1.0]
    kwm={"default":""};
    for s in SCALES:
        if s<1.0: kwm["bai_%.2f"%s]=baimult(s)
    A={a:{x:[] for x in ["BA","TPH","QMD","oBA","oTPH","oQMD"]} for a in arms}; n=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); tr=tr1[tr1.PLT_CN==t1]; tro=tr2[tr2.PLT_CN==int(r.CN)]
        if len(tr)<5 or len(tro)<3: continue
        oBA,oTPH,oQMD=metr(tro)
        if oBA<=0: continue
        sid=str(t1); yrs=int(r.interval)
        std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
        tdf=build_fvs_treeinit(tr,sid); n+=1
        for a in arms:
            try:
                res=run_kw(std,tdf,sid,kwm[a],yrs); s=res.get("summary")
                if s is None or len(s)==0: continue
                l=s.iloc[-1]; ba=g(l,"BA")*M2HA
                if ba<=0: continue
                A[a]["BA"].append(ba); A[a]["TPH"].append(g(l,"Tpa")*TPHc); A[a]["QMD"].append(g(l,"QMD")*CMc)
                A[a]["oBA"].append(oBA); A[a]["oTPH"].append(oTPH); A[a]["oQMD"].append(oQMD)
            except: pass
    print(f"n_run={n}")
    print(f"{'arm':<10}{'BA bias%':>9}{'BA RMSE%':>9}{'TPH bias%':>10}{'QMD bias%':>10}{'n':>5}")
    for a in arms:
        ba=stt(A[a]["BA"],A[a]["oBA"]); t=stt(A[a]["TPH"],A[a]["oTPH"]); q=stt(A[a]["QMD"],A[a]["oQMD"])
        print(f"{a:<10}{ba[1]:>9.1f}{ba[0]:>9.1f}{t[1]:>10.1f}{q[1]:>10.1f}{t[2]:>5d}")
    print("DONE_GROWTH")
if __name__=="__main__": main()
