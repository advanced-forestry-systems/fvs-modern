#!/usr/bin/env python3
"""Fully-calibrated FVS vs default across all variants (2026-06-17).
Calibrated stand = (1) SDIMAX set to the plot's brms site-specific max SDI (metric/2.471 -> English; fall
back to variant median where unmatched); (2) ingrowth injected at the per-variant observed rate
(deployable: rate x interval x initial TPA, no future knowledge); (3) signed BAIMULT for size. Benchmark
on COND-undisturbed FIA remeasurement plots; report BA/TPH/QMD bias default vs calibrated per variant.
Env VARS,NSAMP,SEED,FRAC,BAIMULT."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108; SDIc=2.471
NSAMP=int(os.environ.get("NSAMP","130")); SEEDR=int(os.environ.get("SEED","5")); FRAC=float(os.environ.get("FRAC","0.7")); BAI=float(os.environ.get("BAIMULT","0.90"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/calib_full.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"acd":["ME","NH","VT"],"sn":["AL","GA","SC","MS"],"ls":["MI","MN","WI"],
        "cs":["IL","IN","MO"],"ie":["ID","MT"],"kt":["MT","ID"],"ci":["ID"],"cr":["CO","WY"],"ut":["UT","NV"],
        "ca":["CA"],"nc":["CA","OR"],"ec":["OR","WA"],"wc":["OR","WA"],"pn":["OR","WA"]}
# per-variant observed ingrowth rate (%/decade) from dec.csv; median fallback 30
INGR={"ne":30.5,"acd":36.1,"sn":60.7,"ls":23.5,"kt":69.8,"ci":33.1,"nc":16.9,"ec":21.7,"pn":21.7}
VARS=os.environ.get("VARS","").split(",") if os.environ.get("VARS") else list(VARMAP)
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["key"]=brms.STATECD.astype(str)+"-"+brms.UNITCD.astype(str)+"-"+brms.COUNTYCD.astype(str)+"-"+brms.PLOT.astype(str)
BRMS={k:v/SDIc for k,v in zip(brms.key,brms["SDImax.median"])}   # English SDImax per plot
def load_meta(states):  # CN -> (key, undisturbed?)
    key={}; cls={}
    for s in states:
        pf=Path(FIA)/f"{s}_PLOT.csv"
        if pf.exists():
            pp=pd.read_csv(pf,usecols=lambda c:c in ("CN","STATECD","UNITCD","COUNTYCD","PLOT"),low_memory=False)
            for _,r in pp.iterrows():
                try: key[int(r.CN)]="%d-%d-%d-%d"%(r.STATECD,r.UNITCD,r.COUNTYCD,r.PLOT)
                except: pass
        cf=Path(FIA)/f"{s}_COND.csv"
        if cf.exists():
            c=pd.read_csv(cf,usecols=lambda x:x in ("PLT_CN","TRTCD1","TRTCD2","TRTCD3","DSTRBCD1","DSTRBCD2","DSTRBCD3"),low_memory=False)
            for cn,g2 in c.groupby("PLT_CN"):
                try: cn=int(cn)
                except: continue
                trt=g2[["TRTCD1","TRTCD2","TRTCD3"]].fillna(0).values; dst=g2[["DSTRBCD1","DSTRBCD2","DSTRBCD3"]].fillna(0).values
                cls[cn]="harvested" if (trt==10).any() else ("disturbed" if (dst>0).any() else "undisturbed")
    return key,cls
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc, math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
def sdimax_kw(val): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,val) for i in range(MAXSP))
def baimult_kw(val): return "\n".join("%-16s%10d%10.4f"%("BAIMULT",i+1,val) for i in range(MAXSP))
def seed_rows(tdf,recruits_tpa,spp):
    if recruits_tpa<=0 or not spp: return tdf
    base=int(tdf.tree_id.max())+1 if len(tdf) else 1; per=recruits_tpa/len(spp); rows=[]
    for j,sp in enumerate(spp): rows.append({"stand_id":tdf.stand_id.iloc[0],"plot_id":1,"tree_id":base+j,"tree_count":per,"species":int(sp),"diameter":1.0,"ht":13.0,"crratio":40})
    return pd.concat([tdf,pd.DataFrame(rows)],ignore_index=True)
def run(std,tdf,sid,kw,yrs,var):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=1,cycle_length=yrs))
    lib=os.path.join(FVS_LIB_DIR,"FVS"+var)
    if not os.path.exists(lib): shutil.rmtree(tmp,ignore_errors=True); return None
    r=_run_via_subprocess(lib,os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def bias(f,o):
    m=[(x,y) for x,y in zip(f,o) if y>0]; return float("nan") if not m else 100*sum(x-y for x,y in m)/sum(y for _,y in m)
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states); rate=INGR.get(var,30.0)/100.0
    fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states if (Path(FIA)/f"{s}_PLOT.csv").exists()]
    if not fr: continue
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    key,cls=load_meta(states)
    # per-variant brms fallback = median of matched plots' English SDImax
    fb=[BRMS[key[c]] for c in rem.CN.astype("int64") if c in key and key[c] in BRMS]
    fb=float(np.median(fb)) if fb else None
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    A={k:{m:{"f":[],"o":[]} for m in ("BA","TPH","QMD")} for k in ("default","calibrated")}; nu=0; nsdi=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed": continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)<5 or len(bl)<3: continue
        oBA,oTPH,oQMD=metr(bl)
        if oBA<=0: continue
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
        except Exception: continue
        # calibrated keyword block
        sdival=BRMS.get(key.get(cn,""),fb)
        if key.get(cn,"") in BRMS: nsdi+=1
        kwc=""
        if sdival and sdival>0: kwc+=sdimax_kw(sdival)+"\n"
        kwc+=baimult_kw(BAI)+"\n"
        init_tpa=float(al.TPA_UNADJ.sum())  # trees/acre
        recruits=rate*(yrs/10.0)*init_tpa*FRAC
        spp=list(al.SPCD.value_counts().index[:3])
        tdf_c=seed_rows(tdf,recruits,spp)
        sd=run(std,tdf,sid,"",yrs,var); sc=run(std,tdf_c,sid,kwc,yrs,var)
        if sd is None or len(sd)==0 or sc is None or len(sc)==0: continue
        nu+=1
        for k,s in [("default",sd),("calibrated",sc)]:
            l=s.iloc[-1]
            A[k]["BA"]["f"].append(g(l,"BA")*M2HA); A[k]["BA"]["o"].append(oBA)
            A[k]["TPH"]["f"].append(g(l,"Tpa")*TPHc); A[k]["TPH"]["o"].append(oTPH)
            A[k]["QMD"]["f"].append(g(l,"QMD")*CMc); A[k]["QMD"]["o"].append(oQMD)
    if nu<8: print(var,"too few (n=%d)"%nu); sys.stdout.flush(); continue
    row={"variant":var,"n":nu,"sdimax_matched":nsdi,
         "def_BA%":round(bias(A["default"]["BA"]["f"],A["default"]["BA"]["o"]),1),"cal_BA%":round(bias(A["calibrated"]["BA"]["f"],A["calibrated"]["BA"]["o"]),1),
         "def_TPH%":round(bias(A["default"]["TPH"]["f"],A["default"]["TPH"]["o"]),1),"cal_TPH%":round(bias(A["calibrated"]["TPH"]["f"],A["calibrated"]["TPH"]["o"]),1),
         "def_QMD%":round(bias(A["default"]["QMD"]["f"],A["default"]["QMD"]["o"]),1),"cal_QMD%":round(bias(A["calibrated"]["QMD"]["f"],A["calibrated"]["QMD"]["o"]),1)}
    rows.append(row)
    print("%-4s n=%-3d | BA %+6.1f->%+6.1f | TPH %+6.1f->%+6.1f | QMD %+6.1f->%+6.1f"%(var,nu,row["def_BA%"],row["cal_BA%"],row["def_TPH%"],row["cal_TPH%"],row["def_QMD%"],row["cal_QMD%"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_FULL")
