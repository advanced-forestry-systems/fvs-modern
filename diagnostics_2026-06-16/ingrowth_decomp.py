#!/usr/bin/env python3
"""Decompose the undisturbed TPH gap into ingrowth vs mortality (2026-06-17).
On COND-undisturbed plots per variant: observed ingrowth = TPA of t2 live trees whose (SUBP,TREE) is new
(not present at t1); observed mortality = TPA of t1 live trees not alive at t2. FVS net change = final-initial
TPH (default run, ~no ingrowth). Shows the TPH under-prediction is missing recruitment and quantifies the
observed ingrowth rate (%/decade of initial TPH) to calibrate ESTAB. Env VARS,NSAMP,SEED,OUT."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
TPHc=2.4710538
NSAMP=int(os.environ.get("NSAMP","130")); SEED=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/ingrowth_decomp.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"acd":["ME","NH","VT"],"sn":["AL","GA","SC","MS"],"ls":["MI","MN","WI"],
        "kt":["MT","ID"],"ci":["ID"],"nc":["CA","OR"],"ec":["OR","WA"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","ne,acd,sn,ls,kt,ci,nc,ec,pn").split(",")
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
def run(std,tdf,sid,yrs,var):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=1,cycle_length=yrs))
    r=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+var),os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
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
    oI=[]; oM=[]; oT1=[]; oT2=[]; fNet=[]; fT0=[]; yrs_=[]
    nu=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed": continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]
        if len(a)<5 or len(b)<3: continue
        al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)==0 or len(bl)==0: continue
        t1set=set(zip(a.SUBP,a.TREE))                 # all t1 trees (any status)
        ingrow=bl[~bl.set_index(["SUBP","TREE"]).index.isin(t1set)]    # t2 live not present at t1
        surv_keys=set(zip(bl.SUBP,bl.TREE))
        mort=al[~al.set_index(["SUBP","TREE"]).index.isin(surv_keys)]   # t1 live not live at t2
        ph=lambda df: df.TPA_UNADJ.sum()*TPHc
        oI.append(ph(ingrow)); oM.append(ph(mort)); oT1.append(ph(al)); oT2.append(ph(bl)); yrs_.append(int(r.interval))
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
        except Exception: oI.pop();oM.pop();oT1.pop();oT2.pop();yrs_.pop(); continue
        s=run(std,tdf,sid,yrs,var)
        if s is None or len(s)==0: oI.pop();oM.pop();oT1.pop();oT2.pop();yrs_.pop(); continue
        fT0.append(g(s.iloc[0],"Tpa")*TPHc); fNet.append((g(s.iloc[-1],"Tpa")-g(s.iloc[0],"Tpa"))*TPHc); nu+=1
    if nu<8: print(var,"too few undist (n=%d)"%nu); sys.stdout.flush(); continue
    import numpy as _np
    mI=_np.mean(oI); mM=_np.mean(oM); mT1=_np.mean(oT1); mfNet=_np.mean(fNet); myr=_np.mean(yrs_)
    ingr_rate=100*mI/mT1/myr*10   # %/decade of initial TPH
    obs_net=_np.mean(_np.array(oT2)-_np.array(oT1))
    row={"variant":var,"n":nu,"obs_T1_TPH":round(mT1),"obs_ingrowth_TPH":round(mI),"obs_mort_TPH":round(mM),
         "obs_net_TPH":round(obs_net),"fvs_net_TPH":round(mfNet),"ingrowth_rate_pct_per_decade":round(ingr_rate,1)}
    rows.append(row)
    print("%-4s n=%-3d  obsT1 %4d  obs_ingrowth %4d  obs_mort %4d  obs_net %+4d  FVS_net %+4d  ingrowth %.1f%%/decade"%(
        var,nu,row["obs_T1_TPH"],row["obs_ingrowth_TPH"],row["obs_mort_TPH"],row["obs_net_TPH"],row["fvs_net_TPH"],ingr_rate)); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_DECOMP")
