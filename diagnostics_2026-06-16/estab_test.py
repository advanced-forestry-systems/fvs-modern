#!/usr/bin/env python3
"""ESTAB natural-regeneration proof-of-concept on NE undisturbed plots (2026-06-17).
Inject FVS natural regeneration (ESTAB/NATURAL) of the stand's dominant species at the observed ingrowth
rate and test whether TPH recovers toward observed (without wrecking BA). Arms: default; +ESTAB at ~half
and ~full observed ingrowth density. Species number from FVS_InvReference. Env NSAMP,SEED."""
import os,sys,math,tempfile,sqlite3,shutil
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
TPHc=2.4710538; M2HA=0.2296; TPAc=1/2.4710538
NSAMP=int(os.environ.get("NSAMP","90")); SEED=int(os.environ.get("SEED","5")); VAR="ne"
STATES=["CT","ME","MA","NH","NY","VT"]; FIPS={"CT":9,"ME":23,"MA":25,"NH":33,"NY":36,"VT":50}
INV={v:k for k,v in FIPS.items()}; G.VARIANT_STATES[VAR]=tuple(FIPS.values())
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
def load_cond(states):
    cls={}
    for s in states:
        f=Path(FIA)/f"{s}_COND.csv"
        if not f.exists(): continue
        c=pd.read_csv(f,usecols=lambda x:x in ("PLT_CN","TRTCD1","TRTCD2","TRTCD3","DSTRBCD1","DSTRBCD2","DSTRBCD3"),low_memory=False)
        for cn,g2 in c.groupby("PLT_CN"):
            try: cn=int(cn)
            except: continue
            trt=g2[["TRTCD1","TRTCD2","TRTCD3"]].fillna(0).values; dst=g2[["DSTRBCD1","DSTRBCD2","DSTRBCD3"]].fillna(0).values
            cls[cn]="harvested" if (trt==10).any() else ("disturbed" if (dst>0).any() else "undisturbed")
    return cls
SPCD2IDX={}
def runkw(std,tdf,sid,kw,yrs):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=1,cycle_length=yrs))
    res=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+VAR),os.path.join(tmp,"t.key"),db,tmp)
    if not SPCD2IDX:
        try:
            c=sqlite3.connect(db); ir=pd.read_sql_query("SELECT SpeciesNum,SpeciesFIA FROM FVS_InvReference",c); c.close()
            for _,row in ir.iterrows():
                try: SPCD2IDX[int(row.SpeciesFIA)]=int(row.SpeciesNum)
                except: pass
        except: pass
    shutil.rmtree(tmp,ignore_errors=True); return res.get("summary")
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc)
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def estab(spnum,tpa,yrs):
    # NATURAL: cycle1, species, density(tpa), %surv, age, height(ft); add saplings that grow into tally
    return "ESTAB             0\nNATURAL           1%10d%10.1f     100.0       5.0       4.0\nEND"%(spnum,tpa)
def bias(f,o):
    m=[(x,y) for x,y in zip(f,o) if y>0]; return float("nan") if not m else 100*sum(x-y for x,y in m)/sum(y for _,y in m)
plots=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in STATES]
plot=pd.concat(plots,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=NSAMP,random_state=SEED)
cls=load_cond(STATES)
tr1=G.load_fia_trees(VAR,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(VAR,rem.CN.tolist(),Path(FIA))
# prime species map with one run
ARMS=["default","ESTAB_half","ESTAB_full"]
A={k:{"BA":{"f":[],"o":[]},"TPH":{"f":[],"o":[]}} for k in ARMS}; nu=0
for _,r in rem.iterrows():
    t1=int(r.PREV_PLT_CN); cn=int(r.CN)
    if cls.get(cn)!="undisturbed": continue
    a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
    if len(al)<5 or len(bl)<3: continue
    oBA,oTPH=metr(bl)
    if oBA<=0: continue
    sid=str(t1); yrs=int(r.interval)
    try:
        std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
    except Exception: continue
    # dominant species (FIA SPCD) -> FVS num
    dom=int(al.groupby("SPCD").TPA_UNADJ.sum().idxmax()) if len(al) else 0
    # observed ingrowth tpa for this plot (new trees), used to size NATURAL
    t1set=set(zip(a.SUBP,a.TREE)); ing=bl[~bl.set_index(["SUBP","TREE"]).index.isin(t1set)]
    ing_tpa=ing.TPA_UNADJ.sum()  # trees/acre
    res={}
    sdef=runkw(std,tdf,sid,"",yrs)
    spnum=SPCD2IDX.get(dom,1)
    if sdef is None or len(sdef)==0: continue
    res["default"]=sdef
    res["ESTAB_half"]=runkw(std,tdf,sid,estab(spnum,max(ing_tpa*0.5,5),yrs),yrs)
    res["ESTAB_full"]=runkw(std,tdf,sid,estab(spnum,max(ing_tpa,10),yrs),yrs)
    if any(res[k] is None or len(res[k])==0 for k in ARMS): continue
    nu+=1
    for k in ARMS:
        l=res[k].iloc[-1]; A[k]["BA"]["f"].append(g(l,"BA")*M2HA); A[k]["BA"]["o"].append(oBA); A[k]["TPH"]["f"].append(g(l,"Tpa")*TPHc); A[k]["TPH"]["o"].append(oTPH)
print("NE undisturbed ESTAB proof-of-concept, n=%d"%nu)
for k in ARMS:
    print("%-12s BA %+6.1f  TPH %+6.1f"%(k,bias(A[k]["BA"]["f"],A[k]["BA"]["o"]),bias(A[k]["TPH"]["f"],A[k]["TPH"]["o"])));
print("DONE_ESTAB")
