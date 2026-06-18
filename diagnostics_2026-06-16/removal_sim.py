#!/usr/bin/env python3
"""Removal-simulation converse test (2026-06-17). The decisive proof of the disturbance-artifact headline.
On COND-harvested plots, the default FVS run never cuts -> apparent +40-55% BA over-prediction. Here we
SIMULATE the recorded harvest: remove from the t1 FVS input the trees FIA records as cut (STATUSCD=3 at
t2, or live-at-t1 absent-at-t2), project the residual stand to t2, and compare to observed t2. If the
harvested-plot bias then collapses toward the undisturbed level, harvest (not growth) caused the pooled
over-prediction. Reports harvested-plot BA bias: default (no removal) vs removal-simulated, plus the
undisturbed reference. Env VARS,NSAMP,SEED."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54
NSAMP=int(os.environ.get("NSAMP","180")); SEEDR=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/removal_sim.csv"))
FIPS={"AL":1,"GA":13,"ID":16,"ME":23,"MS":28,"MT":30,"NH":33,"SC":45,"VT":50,"OR":41,"WA":53,"CT":9,"MA":25,"NY":36}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"sn":["AL","GA","SC","MS"],"kt":["MT","ID"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","ne,sn,kt,pn").split(",")
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
def ba(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    return 0 if len(d)==0 else (d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA
def run(std,tdf,sid,yrs,var):
    if len(tdf)==0: return None
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=1,cycle_length=yrs))
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
    fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states]
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    cls=load_cond(states)
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    H={"def":{"f":[],"o":[]},"rem":{"f":[],"o":[]}}; U={"f":[],"o":[]}; nh=0; nu=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN); k=cls.get(cn)
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)<5 or len(bl)<3: continue
        oBA=ba(bl)
        if oBA<=0: continue
        sid=str(t1); yrs=int(r.interval)
        std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
        if k=="undisturbed":
            s=run(std,build_fvs_treeinit(al,sid),sid,yrs,var)
            if s is not None and len(s) and g(s.iloc[-1],"BA")*M2HA>0: U["f"].append(g(s.iloc[-1],"BA")*M2HA); U["o"].append(oBA); nu+=1
            continue
        if k!="harvested": continue
        # cut trees: live at t1, at t2 either STATUSCD==3 (cut) or absent (removed). Use cut+absent as removed.
        b_all=b.set_index(["SUBP","TREE"]);
        def fate(row):
            key=(row.SUBP,row.TREE)
            if key in b_all.index:
                st=b_all.loc[key,"STATUSCD"]; st=st.iloc[0] if hasattr(st,"iloc") else st
                return "cut" if st==3 else "kept"
            return "gone"
        al2=al.copy(); al2["fate"]=al2.apply(fate,axis=1)
        residual=al2[al2.fate=="kept"]              # trees NOT cut/removed = post-harvest residual stand
        if len(residual)<3 or len(residual)==len(al2): continue   # need an actual removal
        sdef=run(std,build_fvs_treeinit(al,sid),sid,yrs,var)        # default: grow everything (no cut)
        srem=run(std,build_fvs_treeinit(residual,sid),sid,yrs,var)  # removal-sim: grow residual only
        if sdef is None or len(sdef)==0 or srem is None or len(srem)==0: continue
        bd=g(sdef.iloc[-1],"BA")*M2HA; br=g(srem.iloc[-1],"BA")*M2HA
        if bd<=0 or br<=0: continue
        nh+=1
        H["def"]["f"].append(bd); H["def"]["o"].append(oBA); H["rem"]["f"].append(br); H["rem"]["o"].append(oBA)
    if nh<8: print(var,"too few harvested (n=%d)"%nh); sys.stdout.flush(); continue
    row={"variant":var,"n_harv":nh,"n_undist":nu,
         "harv_BA_default":round(bias(H["def"]["f"],H["def"]["o"]),1),
         "harv_BA_removalsim":round(bias(H["rem"]["f"],H["rem"]["o"]),1),
         "undist_BA_ref":round(bias(U["f"],U["o"]),1) if nu>=8 else None}
    rows.append(row)
    print("%-4s harv n=%d | BA default(no cut) %+6.1f -> removal-sim %+6.1f | undisturbed ref %s"%(
        var,nh,row["harv_BA_default"],row["harv_BA_removalsim"],str(row["undist_BA_ref"]))); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_REMOVAL")
