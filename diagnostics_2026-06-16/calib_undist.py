#!/usr/bin/env python3
"""Targeted BAIMULT calibration on COND-undisturbed plots (2026-06-17).
For over-predicting variants, sweep all-species BAIMULT on the truly-undisturbed FIA subset and report
BA/TPH/QMD bias per arm. Demonstrates the diameter-growth lever drives BA & QMD bias toward zero (TPH
is orthogonal -> needs recruitment). Env VARS,NSAMP,SEED,OUT."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108
NSAMP=int(os.environ.get("NSAMP","120")); SEED=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/calib_undist.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"sn":["AL","GA","SC","MS"],"kt":["MT","ID"],"nc":["CA","OR"],"ls":["MI","MN","WI"],"cs":["IL","IN","MO"]}
VARS=os.environ.get("VARS","ne,sn,kt,nc").split(",")
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
ARMS=[("default",""),("BAIMULT0.90",allsp("BAIMULT",0.90)),("BAIMULT0.80",allsp("BAIMULT",0.80)),("BAIMULT0.70",allsp("BAIMULT",0.70))]
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
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states if (Path(FIA)/f"{s}_PLOT.csv").exists()]
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEED)
    cls=load_cond(states)
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    A={lab:{m:{"f":[],"o":[]} for m in ("BA","TPH","QMD")} for lab,_ in ARMS}; nu=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed": continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]
        if len(a)<5 or len(b)<3: continue
        oBA,oTPH,oQMD=metr(b)
        if oBA<=0: continue
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(a,sid)
        except Exception: continue
        ok=True; tmpres={}
        for lab,kw in ARMS:
            s=run(std,tdf,sid,kw,yrs,var)
            if s is None or len(s)==0: ok=False; break
            l=s.iloc[-1]; ba=g(l,"BA")*M2HA
            if ba<=0: ok=False; break
            tmpres[lab]=(ba,g(l,"Tpa")*TPHc,g(l,"QMD")*CMc)
        if not ok: continue
        nu+=1
        for lab,_ in ARMS:
            ba,tph,qmd=tmpres[lab]
            A[lab]["BA"]["f"].append(ba); A[lab]["BA"]["o"].append(oBA)
            A[lab]["TPH"]["f"].append(tph); A[lab]["TPH"]["o"].append(oTPH)
            A[lab]["QMD"]["f"].append(qmd); A[lab]["QMD"]["o"].append(oQMD)
    for lab,_ in ARMS:
        bb=bias(A[lab]["BA"]["f"],A[lab]["BA"]["o"]); tt=bias(A[lab]["TPH"]["f"],A[lab]["TPH"]["o"]); qq=bias(A[lab]["QMD"]["f"],A[lab]["QMD"]["o"])
        rows.append({"variant":var,"arm":lab,"n_undist":nu,"BA%":round(bb,1),"TPH%":round(tt,1),"QMD%":round(qq,1)})
        print("%-4s %-12s n=%-3d BA %+6.1f TPH %+6.1f QMD %+6.1f"%(var,lab,nu,bb,tt,qq)); sys.stdout.flush()
    print("---"); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_CALIBU")
