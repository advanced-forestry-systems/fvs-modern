import os,sys,math,json
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"
os.environ["FIA_DATA_DIR"]=FIA; os.environ["FVS_PROJECT_ROOT"]=P; os.environ["FVS_LIB_DIR"]=P+"/lib"; os.environ["FVS_CONFIG_DIR"]=P+"/config"
for p in [P,P+"/calibration/python",P+"/calibration",P+"/deployment/fvs2py"]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import run_fvs_projection, build_fvs_treeinit, build_fvs_standinit
VAR="pn"; STATES=["OR","WA"]; MAXSP=int(os.environ.get("MAXSP","63"))
NSAMP=int(os.environ.get("NSAMP","40")); SEED=int(os.environ.get("SEED","7"))
M2HA=0.2296; TPHc=2.471; CMc=2.54
frames=[]
for ab in STATES:
    f=Path(FIA)/f"{ab}_PLOT.csv"
    if f.exists():
        d=pd.read_csv(f,usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD","UNITCD","COUNTYCD","PLOT"),low_memory=False)
        frames.append(d)
plot=pd.concat(frames,ignore_index=True)
plot["plot_key"]=plot["STATECD"].astype(int).astype(str)+"-"+plot["UNITCD"].astype(int).astype(str)+"-"+plot["COUNTYCD"].astype(int).astype(str)+"-"+plot["PLOT"].astype(int).astype(str)
cn2key=dict(zip(plot["CN"].astype("int64"),plot["plot_key"])); yr=dict(zip(plot["CN"].astype("int64"),plot["MEASYEAR"]))
b=pd.read_csv(P.replace("fvs-modern","fvs-conus")+"/data/brms_SDImax.csv"); b.columns=[c.strip().strip(chr(34)) for c in b.columns]
key2sdi=dict(zip(b["ID"].astype(str),b["SDImax.mean"]))
rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem["PREV_PLT_CN"].astype("int64")
rem["interval"]=rem.apply(lambda r:r["MEASYEAR"]-yr.get(r["PREV_PLT_CN"],np.nan),axis=1)
rem=rem[(rem["interval"]>=5)&(rem["interval"]<=15)]
rem=rem[rem["plot_key"].map(lambda k:k in key2sdi)]
rem=rem.sample(n=min(NSAMP,len(rem)),random_state=SEED)
t1=rem["PREV_PLT_CN"].tolist(); tr1=G.load_fia_trees(VAR,t1,Path(FIA)); tr2=G.load_fia_trees(VAR,rem["CN"].astype("int64").tolist(),Path(FIA))
def metr(df):
    d=df[(df["DIA"]>0)&(df["TPA_UNADJ"]>0)]
    if len(d)==0: return (0,0,0)
    return ((d["TPA_UNADJ"]*0.005454*d["DIA"]**2).sum()*M2HA, d["TPA_UNADJ"].sum()*TPHc, math.sqrt((d["DIA"]**2*d["TPA_UNADJ"]).sum()/d["TPA_UNADJ"].sum())*CMc)
def sdi1(df):
    d=df[(df["DIA"]>0)&(df["TPA_UNADJ"]>0)]
    if len(d)==0: return 0.0
    # metric summation SDI, 25 cm reference, trees/ha
    return float((d["TPA_UNADJ"]*TPHc*(d["DIA"]*CMc/25.0)**1.605).sum())
def sdimax_kw(val): return "\n".join("SDIMAX  %10d%10.1f"%(i,val) for i in range(1,MAXSP+1))
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def runproj(std,tdf,sid,ek,ncyc,clen):
    return run_fvs_projection(std,tdf,sid,VAR,config_version=None,num_cycles=ncyc,cycle_length=clen,extra_keywords=ek)
modes=[("default",""),("localized","brms")]
shortrows=[]; traj=[]
ndef_sdi=[]; nbrms=[]
nrun=0
for _,r in rem.iterrows():
    t1cn=int(r["PREV_PLT_CN"]); tr=tr1[tr1["PLT_CN"]==t1cn].copy(); tro=tr2[tr2["PLT_CN"]==int(r["CN"])]
    if len(tr)<5 or len(tro)<3: continue
    pk=cn2key.get(t1cn) or cn2key.get(int(r["CN"])); brms=key2sdi.get(pk)
    if brms is None or not np.isfinite(brms): continue
    oBA,oTPH,oQMD=metr(tro)
    if oBA<=0: continue
    nbrms.append(brms)
    sid=str(t1cn); yrs=int(r["interval"]); nrun+=1
    std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r["STATECD"])},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
    tdf=build_fvs_treeinit(tr,sid)
    # short-interval validation
    SDI1=sdi1(tr)
    row={"sid":sid,"interval":yrs,"brms":brms,"SDI1":SDI1,"RD_brms":(SDI1/brms if brms>0 else 0),"oBA":oBA,"oTPH":oTPH,"oQMD":oQMD,"oselfthin":(oTPH<sdi1(tr))}
    for lab,m in modes:
        ek=sdimax_kw(brms) if m=="brms" else ""
        try:
            res=runproj(std,tdf,sid,ek,1,yrs); s=res.get("summary")
            if s is not None and len(s)>0:
                l=s.iloc[-1]; ba=g(l,"BA")*M2HA
                if ba>0:
                    row[lab+"_BA"]=ba; row[lab+"_TPH"]=g(l,"Tpa")*TPHc; row[lab+"_QMD"]=g(l,"QMD")*CMc
        except Exception as e: pass
    shortrows.append(row)
    # long-horizon 100yr trajectory
    for lab,m in modes:
        ek=sdimax_kw(brms) if m=="brms" else ""
        try:
            res=runproj(std,tdf,sid,ek,10,10); s=res.get("summary")
            if s is None or len(s)==0: continue
            for _,cy in s.iterrows():
                ba=g(cy,"BA")*M2HA
                traj.append({"sid":sid,"mode":lab,"year":int(g(cy,"Year")),"age":g(cy,"Age"),"BA":ba,"TPH":g(cy,"Tpa")*TPHc,"QMD":g(cy,"QMD")*CMc,"brms":brms})
        except Exception as e: pass
import pandas as pd
S=pd.DataFrame(shortrows); T=pd.DataFrame(traj)
out=P.replace("fvs-modern","fvs-conus")+"/output/pn_sdimax_demo"
os.makedirs(out,exist_ok=True); S.to_csv(out+"/short.csv",index=False); T.to_csv(out+"/traj.csv",index=False)
print("n run:",nrun,"short rows:",len(S),"traj rows:",len(T))
print("mean brms SDImax (tph):",round(np.mean(nbrms),0) if nbrms else None)
def stt(p,o):
    e=[a-b for a,b in zip(p,o) if (a==a and b==b)]; oo=[b for a,b in zip(p,o) if (a==a and b==b)]
    n=len(e)
    if n==0: return (float("nan"),float("nan"),0)
    mo=sum(oo)/n; import math; return (100*math.sqrt(sum(x*x for x in e)/n)/mo,100*(sum(e)/n)/mo,n)
if len(S):
    for mtr in ["BA","TPH","QMD"]:
        print("---",mtr,"vs observed t2 ---")
        for lab,_ in modes:
            c=lab+"_"+mtr
            if c in S:
                sub=S.dropna(subset=[c])
                rm,bi,n=stt(sub[c].tolist(),sub["o"+mtr].tolist())
                print(f"  {lab:10s} RMSE {rm:6.1f}%  bias {bi:6.1f}%  n={n}")
print("DONE_PN_DEMO")
