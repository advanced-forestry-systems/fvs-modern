import os,sys,math
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"
os.environ["FIA_DATA_DIR"]=FIA; os.environ["FVS_PROJECT_ROOT"]=P; os.environ["FVS_LIB_DIR"]=P+"/lib"; os.environ["FVS_CONFIG_DIR"]=P+"/config"
for p in [P,P+"/calibration/python",P+"/calibration",P+"/deployment/fvs2py"]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
G.VARIANT_STATES.setdefault("cr",(8,56,49,35)); G.VARIANT_STATES.setdefault("pn",(41,53)); G.VARIANT_STATES.setdefault("ls",(27,55,26)); G.VARIANT_STATES.setdefault("sn",(1,12,13,28,45,47)); G.VARIANT_STATES.setdefault("ne",(9,23,25,33,36,44,50))
_ABBR={1:"AL",6:"CA",8:"CO",9:"CT",12:"FL",13:"GA",16:"ID",17:"IL",18:"IN",19:"IA",23:"ME",25:"MA",26:"MI",27:"MN",28:"MS",29:"MO",30:"MT",32:"NV",33:"NH",35:"NM",36:"NY",41:"OR",44:"RI",45:"SC",47:"TN",49:"UT",50:"VT",53:"WA",55:"WI",56:"WY"}; _o=G._state_abbrev
G._state_abbrev=lambda c:_ABBR.get(c) or _o(c)
from perseus_100yr_projection import run_fvs_projection, build_fvs_treeinit, build_fvs_standinit
import os as _os
VAR=_os.environ.get("VAR","cr"); STATES=_os.environ.get("STATES","CO,WY,UT,NM").split(","); MAXSP=int(_os.environ.get("MAXSP","38")); M2HA=0.2296;TPHc=2.471;CMc=2.54
SCALES=[0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
frames=[]
for ab in STATES:
    f=Path(FIA)/f"{ab}_PLOT.csv"
    if f.exists(): frames.append(pd.read_csv(f,usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD","UNITCD","COUNTYCD","PLOT"),low_memory=False))
plot=pd.concat(frames,ignore_index=True)
plot["plot_key"]=plot.STATECD.astype(int).astype(str)+"-"+plot.UNITCD.astype(int).astype(str)+"-"+plot.COUNTYCD.astype(int).astype(str)+"-"+plot.PLOT.astype(int).astype(str)
cn2key=dict(zip(plot.CN.astype("int64"),plot.plot_key)); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
b=pd.read_csv(P.replace("fvs-modern","fvs-conus")+"/data/brms_SDImax.csv"); b.columns=[c.strip().strip(chr(34)) for c in b.columns]
key2sdi=dict(zip(b.ID.astype(str),b["SDImax.mean"]))
rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64")
rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
rem=rem[(rem.interval>=5)&(rem.interval<=15)]; rem=rem[rem.plot_key.map(lambda k:k in key2sdi)]
rem=rem.sample(n=min(110,len(rem)),random_state=5)
tr1=G.load_fia_trees(VAR,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(VAR,rem.CN.astype("int64").tolist(),Path(FIA))
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc)
def kw(v): return "\n".join("SDIMAX  %10d%10.1f"%(i,v/TPHc) for i in range(1,MAXSP+1))
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
cols=["default"]+["brms_x%.1f"%s for s in SCALES]
acc={c:{"p":[],"o":[]} for c in cols}; n=0
for _,r in rem.iterrows():
    t1cn=int(r.PREV_PLT_CN); tr=tr1[tr1.PLT_CN==t1cn]; tro=tr2[tr2.PLT_CN==int(r.CN)]
    if len(tr)<5 or len(tro)<3: continue
    brms=key2sdi.get(cn2key.get(t1cn))
    if brms is None or not np.isfinite(brms): continue
    oBA,oTPH=metr(tro)
    if oBA<=0: continue
    sid=str(t1cn); yrs=int(r.interval); n+=1
    std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD)},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
    tdf=build_fvs_treeinit(tr,sid)
    settings=[("default","")]+[("brms_x%.1f"%s, kw(brms*s)) for s in SCALES]
    for lab,ek in settings:
        try:
            res=run_fvs_projection(std,tdf,sid,VAR,config_version=None,num_cycles=1,cycle_length=yrs,extra_keywords=ek)
            s=res.get("summary")
            if s is None or len(s)==0: continue
            l=s.iloc[-1]; tph=g(l,"Tpa")*TPHc
            if tph>0: acc[lab]["p"].append(tph); acc[lab]["o"].append(oTPH)
        except: pass
def stt(p,o):
    m=[(a,b) for a,b in zip(p,o)]; k=len(m); e=[a-b for a,b in m]; mo=sum(b for _,b in m)/k
    return 100*math.sqrt(sum(x*x for x in e)/k)/mo, 100*sum(e)/k/mo, k
print("CR scale diagnostic, n=",n,"mean brms=",round(np.mean([key2sdi[cn2key[int(r.PREV_PLT_CN)]] for _,r in rem.iterrows() if cn2key.get(int(r.PREV_PLT_CN)) in key2sdi]),0))
print("%-12s %9s %7s n" % ("setting","TPH RMSE%","bias%"))
for c in cols:
    a=acc[c]
    if a["p"]: rm,bi,k=stt(a["p"],a["o"]); print(f"{c:12s} {rm:9.1f} {bi:7.1f} {k}")
print("DONE_CR_SCALE")
