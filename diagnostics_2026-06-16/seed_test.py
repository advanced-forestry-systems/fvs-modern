#!/usr/bin/env python3
"""Ingrowth-injection demonstration on undisturbed plots (2026-06-17).
Since FVS ESTAB cannot add background ingrowth undisturbed, inject the recruitment cohort directly: append
small recruit trees to the initial treelist at the observed per-plot ingrowth rate, let FVS grow them.
Arms: default; +SEED (recruits at observed ingrowth TPA); +SEED+BAIMULT0.85. Report BA/TPH/QMD bias.
Env VARS,NSAMP,SEED_FRAC,SEED."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
TPHc=2.4710538; M2HA=0.2296; TPAac=1/2.4710538; CMc=2.54; MAXSP=108
NSAMP=int(os.environ.get("NSAMP","120")); SEEDR=int(os.environ.get("SEED","5")); FRAC=float(os.environ.get("SEED_FRAC","0.7"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/seed_test.csv"))
FIPS={"AL":1,"CA":6,"CT":9,"GA":13,"ID":16,"ME":23,"MA":25,"MS":28,"MT":30,"NH":33,"NY":36,"OR":41,"SC":45,"VT":50,"WA":53}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"sn":["AL","GA","SC","MS"],"kt":["MT","ID"]}
VARS=os.environ.get("VARS","ne,sn,kt").split(",")
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
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
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc, math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
def run(std,tdf,sid,kw,yrs,var):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=1,cycle_length=yrs))
    r=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+var),os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def bias(f,o):
    m=[(x,y) for x,y in zip(f,o) if y>0]; return float("nan") if not m else 100*sum(x-y for x,y in m)/sum(y for _,y in m)
def seed_rows(tdf, recruits_tpa, spp):
    # append small recruit trees (1.0 in, 13 ft, cr 40) at recruits_tpa, split across up to 3 records
    if recruits_tpa<=0: return tdf
    base_id=int(tdf.tree_id.max())+1 if len(tdf) else 1; rows=[]
    per=recruits_tpa/ max(len(spp),1)
    for j,sp in enumerate(spp):
        rows.append({"stand_id":tdf.stand_id.iloc[0],"plot_id":1,"tree_id":base_id+j,"tree_count":per,"species":int(sp),"diameter":1.0,"ht":13.0,"crratio":40})
    return pd.concat([tdf,pd.DataFrame(rows)],ignore_index=True)
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    plots=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states]
    plot=pd.concat(plots,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    cls=load_cond(states)
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    ARMS=["default","SEED","SEED+BAI0.85"]
    A={k:{m:{"f":[],"o":[]} for m in ("BA","TPH","QMD")} for k in ARMS}; nu=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed": continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)<5 or len(bl)<3: continue
        oBA,oTPH,oQMD=metr(bl)
        if oBA<=0: continue
        t1set=set(zip(a.SUBP,a.TREE)); ing=bl[~bl.set_index(["SUBP","TREE"]).index.isin(t1set)]
        ing_tpa=ing.TPA_UNADJ.sum()*FRAC                  # trees/acre to seed (fraction of observed ingrowth)
        spp=list(ing.SPCD.value_counts().index[:3]) if len(ing) else list(al.SPCD.value_counts().index[:1])
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
        except Exception: continue
        tdf_seed=seed_rows(tdf,ing_tpa,spp)
        res={}
        res["default"]=run(std,tdf,sid,"",yrs,var)
        res["SEED"]=run(std,tdf_seed,sid,"",yrs,var)
        res["SEED+BAI0.85"]=run(std,tdf_seed,sid,allsp("BAIMULT",0.85),yrs,var)
        if any(res[k] is None or len(res[k])==0 for k in ARMS): continue
        nu+=1
        for k in ARMS:
            l=res[k].iloc[-1]
            A[k]["BA"]["f"].append(g(l,"BA")*M2HA); A[k]["BA"]["o"].append(oBA)
            A[k]["TPH"]["f"].append(g(l,"Tpa")*TPHc); A[k]["TPH"]["o"].append(oTPH)
            A[k]["QMD"]["f"].append(g(l,"QMD")*CMc); A[k]["QMD"]["o"].append(oQMD)
    for k in ARMS:
        rows.append({"variant":var,"arm":k,"n":nu,"BA%":round(bias(A[k]["BA"]["f"],A[k]["BA"]["o"]),1),"TPH%":round(bias(A[k]["TPH"]["f"],A[k]["TPH"]["o"]),1),"QMD%":round(bias(A[k]["QMD"]["f"],A[k]["QMD"]["o"]),1)})
        print("%-4s %-13s n=%-3d BA %+6.1f TPH %+6.1f QMD %+6.1f"%(var,k,nu,rows[-1]["BA%"],rows[-1]["TPH%"],rows[-1]["QMD%"])); sys.stdout.flush()
    print("---"); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_SEED")
