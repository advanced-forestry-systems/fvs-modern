import os, sys, math
P="/users/PUOM0008/crsfaaron/fvs-modern"
os.environ["FIA_DATA_DIR"]="/fs/scratch/PUOM0008/crsfaaron/FIA"
os.environ["FVS_PROJECT_ROOT"]=P; os.environ["FVS_LIB_DIR"]=P+"/lib"; os.environ["FVS_CONFIG_DIR"]=P+"/config"
for p in [P, P+"/calibration/python", P+"/calibration", P+"/deployment/fvs2py"]: sys.path.insert(0,p)
import pandas as pd
from fia_stand_generator import build_stand_init, build_tree_init
from perseus_100yr_projection import run_fvs_projection
B="/users/PUOM0008/crsfaaron/fvs-conus/output/conus/stand_level/"
trees=pd.read_csv(B+"ne_bench_trees.csv"); plots=pd.read_csv(B+"ne_bench_plots.csv")
print("plots:",len(plots)," trees:",len(trees))
M2HA=0.2296; TPHc=2.4710538; CMc=2.54
runs=[("NE",None,"real_NE"),("ACD",None,"real_ACD"),("NE","calibrated","unified")]
res={lab:{"BA":[], "TPH":[], "QMD":[], "obsBA":[], "obsTPH":[], "obsQMD":[]} for _,_,lab in runs}
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
for _,prow in plots.iterrows():
    pk=str(prow["plot_key"]); yrs=int(prow["YEARS"]) if prow["YEARS"]>=1 else 10
    cond=pd.Series({"STDAGE":prow.get("STDAGE",40),"ASPECT":prow.get("ASPECT",0),"SLOPE":prow.get("SLOPE",15),
                    "ELEV":prow.get("ELEV",1000),"SICOND":prow.get("SICOND",60),"STATECD":prow.get("STATECD",0),
                    "COUNTYCD":prow.get("COUNTYCD",0),"FORTYPCD":prow.get("FORTYPCD",0)})
    tr=trees[trees["plot_key"]==prow["plot_key"]].copy()
    if tr.empty: continue
    std=build_stand_init(cond, pk, "NE"); std["num_plots"]=4; tdf=build_tree_init(tr, pk)
    for variant,cfg,label in runs:
        try:
            r=run_fvs_projection(std,tdf,pk,variant,config_version=cfg,num_cycles=1,cycle_length=yrs)
            s=r.get("summary")
            if s is None or len(s)==0: continue
            f0=s.iloc[0]; l=s.iloc[-1]; ba=g(l,"BA")*M2HA
            res[label].setdefault("BA0",[]).append(g(f0,"BA")*M2HA); res[label].setdefault("TPH0",[]).append(g(f0,"Tpa")*TPHc); res[label].setdefault("oBA1",[]).append(float(prow["BA1"])); res[label].setdefault("oTPH1",[]).append(float(prow["TPH1"]))
            if ba<=0: continue
            res[label]["BA"].append(ba); res[label]["TPH"].append(g(l,"Tpa")*TPHc); res[label]["QMD"].append(g(l,"QMD")*CMc)
            res[label]["obsBA"].append(float(prow["BA2_obs"])); res[label]["obsTPH"].append(float(prow["TPH2_obs"])); res[label]["obsQMD"].append(float(prow["QMD2_obs"]))
        except Exception as e: pass
def stats(pred,obs):
    n=len(pred);
    if n==0: return (0,0,0)
    err=[p-o for p,o in zip(pred,obs)]; mo=sum(obs)/n
    rmse=math.sqrt(sum(e*e for e in err)/n); bias=sum(err)/n
    return (100*rmse/mo, 100*bias/mo, n)
print("\n=== NE/ACD vs unified vs OBSERVED (pctRMSE / pctBias, n) ===")
print(f"{'model':9s} | {'BA pctRMSE':>10s} {'BA bias':>8s} | {'TPH pctRMSE':>11s} {'TPH bias':>8s} | {'QMD pctRMSE':>11s} {'QMD bias':>8s}")
for _,_,lab in runs:
    r=res[lab]
    ba=stats(r["BA"],r["obsBA"]); tph=stats(r["TPH"],r["obsTPH"]); qmd=stats(r["QMD"],r["obsQMD"])
    print(f"{lab:9s} | {ba[0]:9.1f}% {ba[1]:7.1f}% | {tph[0]:10.1f}% {tph[1]:7.1f}% | {qmd[0]:10.1f}% {qmd[1]:7.1f}%  (n={ba[2]})")
print("\n=== INITIAL (FVS year 0) vs observed t1 [input-loading check] ===")
for _,_,lab in runs:
    r=res[lab]
    if r.get("BA0"):
        b=stats(r["BA0"],r["oBA1"]); t=stats(r["TPH0"],r["oTPH1"])
        print(f"{lab:9s} | BA0 %RMSE {b[0]:6.1f} bias {b[1]:6.1f} | TPH0 %RMSE {t[0]:6.1f} bias {t[1]:6.1f}")
        # show a couple raw rows
import numpy as np
r=res["real_NE"]
print("sample year0 BA/TPH vs obs t1 (first 5):")
for i in range(min(5,len(r.get("BA0",[])))):
    b0=r["BA0"][i]; t0=r["TPH0"][i]; ob=r["oBA1"][i]; ot=r["oTPH1"][i]
    print("  FVS0 BA=%5.1f TPH=%6.0f | obs t1 BA=%5.1f TPH=%6.0f" % (b0,t0,ob,ot))
print("DONE_OBS_BENCH")
