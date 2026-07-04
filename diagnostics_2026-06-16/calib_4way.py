#!/usr/bin/env python3
"""FULL STRESS TEST: default vs fully-calibrated FVS, all variants, all metrics (2026-06-17).
Calibrated = brms site-specific SDImax + SIGN-AWARE ingrowth (only where variant under-predicts TPH) +
BAIMULT 0.90. Metrics: BA, TPH, QMD, merch volume (FVS MCuFt vs FIA VOLCFNET). COND-undisturbed FIA
remeasurement. Env VARS,NSAMP,SEED,OUT."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108; SDIc=2.471; VOLc=0.0699055
NSAMP=int(os.environ.get("NSAMP","130")); SEEDR=int(os.environ.get("SEED","5")); FRAC=0.7; BAI=0.90
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/calib_final.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"acd":["ME","NH","VT"],"sn":["AL","GA","SC","MS"],"ls":["MI","MN","WI"],
        "cs":["IL","IN","MO"],"ie":["ID","MT"],"kt":["MT","ID"],"ci":["ID"],"cr":["CO","WY"],"ut":["UT","NV"],
        "ca":["CA"],"nc":["CA","OR"],"ec":["OR","WA"],"wc":["OR","WA"],"pn":["OR","WA"]}
INGR={"ne":30.5,"acd":36.1,"sn":60.7,"ls":23.5,"kt":69.8,"ci":33.1,"nc":16.9,"ec":21.7,"pn":21.7,"ie":25.0,"ut":30.0,"wc":21.7}
INJECT={"ne","acd","sn","ie","kt","ci","ut","nc","ec","wc","pn"}  # sign-aware: variants where default TPH is under-predicted
VARS=os.environ.get("VARS","").split(",") if os.environ.get("VARS") else list(VARMAP)
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["key"]=brms.STATECD.astype(str)+"-"+brms.UNITCD.astype(str)+"-"+brms.COUNTYCD.astype(str)+"-"+brms.PLOT.astype(str)
BRMS={k:v/SDIc for k,v in zip(brms.key,brms["SDImax.median"])}
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
def sdimax_kw(val): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,val) for i in range(MAXSP))
def seed_rows(tdf,rt,spp):
    if rt<=0 or not spp: return tdf
    base=int(tdf.tree_id.max())+1; per=rt/len(spp); rows=[{"stand_id":tdf.stand_id.iloc[0],"plot_id":1,"tree_id":base+j,"tree_count":per,"species":int(sp),"diameter":1.0,"ht":13.0,"crratio":40} for j,sp in enumerate(spp)]
    return pd.concat([tdf,pd.DataFrame(rows)],ignore_index=True)
def load_meta(states):
    key={}; cls={}; obs={}
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
        tf=Path(FIA)/f"{s}_TREE.csv"
        if tf.exists():
            t=pd.read_csv(tf,usecols=lambda c:c in ("PLT_CN","STATUSCD","DIA","TPA_UNADJ","VOLCFNET"),low_memory=False)
            t=t[(t.STATUSCD==1)&(t.DIA>0)&(t.TPA_UNADJ>0)]; t["v"]=t.VOLCFNET.fillna(0)*t.TPA_UNADJ
            for cn,gg in t.groupby("PLT_CN"):
                try: obs[int(cn)]=gg.v.sum()*VOLc
                except: pass
    return key,cls,obs
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc, math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
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
rows=[]; cond_rows=[]   # cond_rows: per-condition obs/default/calibrated, keyed by CN for the 4-way join
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states); rate=INGR.get(var,30.0)/100.0; inject=var in INJECT
    fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states if (Path(FIA)/f"{s}_PLOT.csv").exists()]
    if not fr: continue
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    key,cls,obs=load_meta(states)
    fb=[BRMS[key[c]] for c in rem.CN.astype("int64") if c in key and key[c] in BRMS]; fb=float(np.median(fb)) if fb else None
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    A={k:{m:{"f":[],"o":[]} for m in ("BA","TPH","QMD","VOL")} for k in ("default","calibrated")}; nu=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed": continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)<5 or len(bl)<3 or cn not in obs: continue
        oBA,oTPH,oQMD=metr(bl); oVOL=obs[cn]
        if oBA<=0 or oVOL<=0: continue
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
        except Exception: continue
        sdival=BRMS.get(key.get(cn,""),fb); kwc=""
        if sdival and sdival>0: kwc+=sdimax_kw(sdival)+"\n"
        kwc+=allsp("BAIMULT",BAI)+"\n"
        tdf_c=tdf
        if inject: tdf_c=seed_rows(tdf,rate*(yrs/10.0)*float(al.TPA_UNADJ.sum())*FRAC,list(al.SPCD.value_counts().index[:3]))
        sd=run(std,tdf,sid,"",yrs,var); sc=run(std,tdf_c,sid,kwc,yrs,var)
        if sd is None or len(sd)==0 or sc is None or len(sc)==0: continue
        nu+=1
        for k,s in [("default",sd),("calibrated",sc)]:
            l=s.iloc[-1]
            A[k]["BA"]["f"].append(g(l,"BA")*M2HA); A[k]["BA"]["o"].append(oBA)
            A[k]["TPH"]["f"].append(g(l,"Tpa")*TPHc); A[k]["TPH"]["o"].append(oTPH)
            A[k]["QMD"]["f"].append(g(l,"QMD")*CMc); A[k]["QMD"]["o"].append(oQMD)
            A[k]["VOL"]["f"].append(g(l,"MCuFt")*VOLc); A[k]["VOL"]["o"].append(oVOL)
        ld=sd.iloc[-1]; lc=sc.iloc[-1]
        cond_rows.append({"cn":cn,"t1":t1,"variant":var,
            "oBA":oBA,"dBA":g(ld,"BA")*M2HA,"cBA":g(lc,"BA")*M2HA,
            "oTPH":oTPH,"dTPH":g(ld,"Tpa")*TPHc,"cTPH":g(lc,"Tpa")*TPHc,
            "oQMD":oQMD,"dQMD":g(ld,"QMD")*CMc,"cQMD":g(lc,"QMD")*CMc,
            "oVOL":oVOL,"dVOL":g(ld,"MCuFt")*VOLc,"cVOL":g(lc,"MCuFt")*VOLc})
    if nu<8: print(var,"too few (n=%d)"%nu); sys.stdout.flush(); continue
    row={"variant":var,"n":nu,"inject":int(inject)}
    for m in ("BA","TPH","QMD","VOL"):
        row["def_%s"%m]=round(bias(A["default"][m]["f"],A["default"][m]["o"]),1)
        row["cal_%s"%m]=round(bias(A["calibrated"][m]["f"],A["calibrated"][m]["o"]),1)
    rows.append(row)
    print("%-4s n=%-3d inj=%d | BA %+5.1f>%+5.1f TPH %+6.1f>%+6.1f QMD %+5.1f>%+5.1f VOL %+5.1f>%+5.1f"%(
        var,nu,inject,row["def_BA"],row["cal_BA"],row["def_TPH"],row["cal_TPH"],row["def_QMD"],row["cal_QMD"],row["def_VOL"],row["cal_VOL"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
if cond_rows:
    cout=OUT.replace(".csv","_percond.csv")
    with open(cout,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(cond_rows[0].keys())); w.writeheader(); w.writerows(cond_rows)
    print("wrote",cout,"(%d conditions)"%len(cond_rows))
# summary
import numpy as _np
def med(k):
    v=[r[k] for r in rows if r[k]==r[k]]; return _np.median(v) if v else float("nan")
print("MEDIAN |bias| default vs calibrated:")
for m in ("BA","TPH","QMD","VOL"):
    d=_np.median([abs(r["def_%s"%m]) for r in rows]); c=_np.median([abs(r["cal_%s"%m]) for r in rows]); print("  %-3s  %.1f -> %.1f"%(m,d,c))
print("DONE_FINAL")
