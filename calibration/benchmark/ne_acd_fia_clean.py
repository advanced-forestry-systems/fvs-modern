import os, sys, math
from pathlib import Path
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA=Path("/fs/scratch/PUOM0008/crsfaaron/FIA")
os.environ["FIA_DATA_DIR"]=str(FIA); os.environ["FVS_PROJECT_ROOT"]=P
os.environ["FVS_LIB_DIR"]=P+"/lib"; os.environ["FVS_CONFIG_DIR"]=P+"/config"
for p in [P, P+"/calibration/python", P+"/calibration", P+"/deployment/fvs2py"]: sys.path.insert(0,p)
import pandas as pd, numpy as np
import fia_stand_generator as G
from perseus_100yr_projection import run_fvs_projection, build_fvs_treeinit, build_fvs_standinit
NE_ABBR=["CT","ME","MA","NH","NY","RI","VT"]
# 1) build remeasurement linkage from PLOT tables
plotcols=["CN","PREV_PLT_CN","MEASYEAR","STATECD","SICOND"]
frames=[]
for ab in NE_ABBR:
    f=FIA/f"{ab}_PLOT.csv"
    if not f.exists(): continue
    df=pd.read_csv(f, usecols=lambda c: c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"), low_memory=False)
    frames.append(df)
plot=pd.concat(frames, ignore_index=True)
yr=dict(zip(plot["CN"].astype("int64"), plot["MEASYEAR"]))
rem=plot.dropna(subset=["PREV_PLT_CN"]).copy()
rem["PREV_PLT_CN"]=rem["PREV_PLT_CN"].astype("int64"); rem["CN"]=rem["CN"].astype("int64")
rem["interval"]=rem.apply(lambda r: r["MEASYEAR"]-yr.get(r["PREV_PLT_CN"], np.nan), axis=1)
rem=rem[(rem["interval"]>=5)&(rem["interval"]<=15)]
rem=rem.sample(n=min(35,len(rem)), random_state=11)
print("remeasured NE plots sampled:", len(rem))
t1_cns=rem["PREV_PLT_CN"].tolist(); t2_cns=rem["CN"].tolist()
tr1=G.load_fia_trees("NE", t1_cns, FIA); tr2=G.load_fia_trees("NE", t2_cns, FIA)
print("t1 trees:", len(tr1), " t2 trees:", len(tr2))
def metr(df):  # FIA trees -> metric stand metrics
    d=df[(df["DIA"]>0)&(df["TPA_UNADJ"]>0)].copy()
    if d.empty: return (0,0,0)
    ba=(d["TPA_UNADJ"]*0.005454*d["DIA"]**2).sum()*0.2296
    tph=d["TPA_UNADJ"].sum()*2.4710538
    qmd=math.sqrt((d["DIA"]**2*d["TPA_UNADJ"]).sum()/d["TPA_UNADJ"].sum())*2.54
    return (ba,tph,qmd)
M2HA=0.2296; TPHc=2.4710538; CMc=2.54
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
runs=[("NE",None,"real_NE"),("ACD",None,"real_ACD"),("NE","calibrated","unified")]
A={lab:{k:[] for k in ["BA","TPH","QMD","oBA","oTPH","oQMD","BA0","TPH0","oBA1","oTPH1"]} for _,_,lab in runs}
nrun=0
for _,r in rem.iterrows():
    t1=int(r["PREV_PLT_CN"]); tr=tr1[tr1["PLT_CN"]==t1].copy(); t2=int(r["CN"]); tro=tr2[tr2["PLT_CN"]==t2]
    if len(tr)<5 or len(tro)<3: continue
    oBA1,oTPH1,_=metr(tr); oBA2,oTPH2,oQMD2=metr(tro)
    if oBA2<=0: continue
    cond=pd.Series({"STDAGE":40,"ASPECT":0,"SLOPE":15,"ELEV":1000,"SICOND":60,
                    "STATECD":int(r["STATECD"]),"COUNTYCD":0,"FORTYPCD":0})
    sid=str(t1); yrs=int(r["interval"])
    std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r["STATECD"]),"COUNTYCD":0}, sid, "NE"); tdf=build_fvs_treeinit(tr, sid); nrun+=1
    for variant,cfg,label in runs:
        try:
            res=run_fvs_projection(std,tdf,sid,variant,config_version=cfg,num_cycles=1,cycle_length=yrs)
            s=res.get("summary")
            if s is None or len(s)==0: continue
            f0=s.iloc[0]; l=s.iloc[-1]; ba=g(l,"BA")*M2HA
            if ba<=0: continue
            a=A[label]
            a["BA"].append(ba); a["TPH"].append(g(l,"Tpa")*TPHc); a["QMD"].append(g(l,"QMD")*CMc)
            a["oBA"].append(oBA2); a["oTPH"].append(oTPH2); a["oQMD"].append(oQMD2)
            a["BA0"].append(g(f0,"BA")*M2HA); a["TPH0"].append(g(f0,"Tpa")*TPHc); a["oBA1"].append(oBA1); a["oTPH1"].append(oTPH1)
        except Exception as e: pass
def stt(pred,obs):
    n=len(pred)
    if n==0: return (0,0,0)
    err=[a-b for a,b in zip(pred,obs)]; mo=sum(obs)/n
    return (100*math.sqrt(sum(e*e for e in err)/n)/mo, 100*(sum(err)/n)/mo, n)
print(f"\nplots run: {nrun}")
print("=== YEAR-0 input check (FVS0 vs observed t1, FIA-derived) ===")
for _,_,lab in runs:
    a=A[lab]; b=stt(a["BA0"],a["oBA1"]); t=stt(a["TPH0"],a["oTPH1"])
    if b[2]: print(f"{lab:9s} | BA0 %RMSE {b[0]:6.1f} bias {b[1]:6.1f} | TPH0 %RMSE {t[0]:6.1f} bias {t[1]:6.1f} (n={b[2]})")
print("\n=== PROJECTED vs OBSERVED t2 (CLEAN VERDICT) ===")
print(f"{'model':9s} | {'BA %RMSE':>8s} {'BA bias':>8s} | {'TPH %RMSE':>9s} {'TPH bias':>8s} | {'QMD %RMSE':>9s} {'QMD bias':>8s}")
for _,_,lab in runs:
    a=A[lab]; ba=stt(a["BA"],a["oBA"]); tph=stt(a["TPH"],a["oTPH"]); qmd=stt(a["QMD"],a["oQMD"])
    print(f"{lab:9s} | {ba[0]:7.1f}% {ba[1]:7.1f}% | {tph[0]:8.1f}% {tph[1]:7.1f}% | {qmd[0]:8.1f}% {qmd[1]:7.1f}%  (n={ba[2]})")
print("DONE_CLEAN")
