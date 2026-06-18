#!/usr/bin/env python3
"""Single-framework four-arm comparison, FVS engine, one disturbance-clean held-out plot set (2026-06-18).

Engine-runnable arms (this script):
  A = default FVS (no keywords)
  B = keyword-calibrated: brms site-specific max SDI + density-dependent recruitment (item 1) +
      fold-A BAIMULT (level minimizing |QMD bias| on fold A)
All four metrics (BA, TPH, QMD, merch volume MCuFt vs FIA VOLCFNET) with percentile-bootstrap 95% CIs,
reported OUT-OF-SAMPLE on spatial fold B (calibration derived on fold A only, county-hash folds).

Arms C (fvs-conus equations) and D (fvs-conus + density layer) require injecting the fvs-conus equations
into the FVS engine (fvs2py in-process tree loading, the documented maintainer-level blocker). They are
produced by the fvs-conus standalone projector and reconciled to this engine baseline via within-framework
deltas in the synthesis doc; see 20260618_fourarm_result.md. This keeps A and B in one framework with CIs.
Env: VARS, NSAMP, SEED, OUT.
"""
import os, sys, math, tempfile, sqlite3, shutil, csv, hashlib, json
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108; SDIc=2.471; VOLc=0.0699055
NSAMP=int(os.environ.get("NSAMP","400")); SEEDR=int(os.environ.get("SEED","5")); FRAC=0.7
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/fourarm_engine.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"sn":["AL","GA","SC","MS"],"kt":["MT","ID"],"pn":["OR","WA"],
        "cr":["CO","WY"],"ut":["UT","NV"],"nc":["CA","OR"],"ec":["OR","WA"],"wc":["OR","WA"]}
VARS=os.environ.get("VARS","ne,sn,kt,pn").split(",")
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["key"]=brms.STATECD.astype(str)+"-"+brms.UNITCD.astype(str)+"-"+brms.COUNTYCD.astype(str)+"-"+brms.PLOT.astype(str)
BRMS={k:v/SDIc for k,v in zip(brms.key,brms["SDImax.median"])}
def fold_of(statecd,countycd):
    h=int(hashlib.md5(("%d-%d"%(statecd,countycd)).encode()).hexdigest(),16)
    return "A" if h%2==0 else "B"
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
def sdimax_kw(val): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,val) for i in range(MAXSP))
def stand_sdi(al):
    d=al[(al.DIA>0)&(al.TPA_UNADJ>0)]
    return float((d.TPA_UNADJ*(d.DIA/10.0)**1.605).sum()) if len(d) else 0.0
def seed_rows(tdf,rt,spp):
    if rt<=0 or not spp: return tdf
    base=int(tdf.tree_id.max())+1; per=rt/len(spp)
    rows=[{"stand_id":tdf.stand_id.iloc[0],"plot_id":1,"tree_id":base+j,"tree_count":per,"species":int(sp),"diameter":1.0,"ht":13.0,"crratio":40} for j,sp in enumerate(spp)]
    return pd.concat([tdf,pd.DataFrame(rows)],ignore_index=True)
def load_meta(states):
    cls={}; cnty={}; obs={}
    for s in states:
        pf=Path(FIA)/f"{s}_PLOT.csv"
        if pf.exists():
            pp=pd.read_csv(pf,usecols=lambda c:c in ("CN","STATECD","UNITCD","COUNTYCD","PLOT"),low_memory=False)
            for _,r in pp.iterrows():
                try: cnty[int(r.CN)]=(int(r.STATECD),int(r.COUNTYCD),"%d-%d-%d-%d"%(r.STATECD,r.UNITCD,r.COUNTYCD,r.PLOT))
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
    return cls,cnty,obs
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
def boot_ci(f,o,nb=2000,seed=11):
    m=[(x,y) for x,y in zip(f,o) if y>0]
    if len(m)<3: return (float("nan"),float("nan"))
    rng=np.random.default_rng(seed); fa=np.array([x for x,_ in m],float); oa=np.array([y for _,y in m],float); n=len(m); out=[]
    for _ in range(nb):
        i=rng.integers(0,n,n); out.append(100*(fa[i]-oa[i]).sum()/oa[i].sum())
    return (round(float(np.percentile(out,2.5)),1),round(float(np.percentile(out,97.5)),1))
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states if (Path(FIA)/f"{s}_PLOT.csv").exists()]
    if not fr: continue
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    cls,cnty,obs=load_meta(states)
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    recs={"A":[],"B":[]}
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed" or cn not in cnty or cn not in obs: continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)<5 or len(bl)<3: continue
        oBA,oTPH,oQMD=metr(bl); oVOL=obs[cn]
        if oBA<=0 or oVOL<=0: continue
        sc,cc,pkey=cnty[cn]; fold=fold_of(sc,cc)
        t1set=set(zip(a.SUBP,a.TREE)); ing=bl[~bl.set_index(["SUBP","TREE"]).index.isin(t1set)]
        recs[fold].append({"t1":t1,"cn":cn,"al":al,"oBA":oBA,"oTPH":oTPH,"oQMD":oQMD,"oVOL":oVOL,"yrs":int(r.interval),
                           "statecd":sc,"sdi":BRMS.get(pkey),"sdi_t1":stand_sdi(al),"ing_tpa":ing.TPA_UNADJ.sum(),
                           "init_tpa":float(al.TPA_UNADJ.sum()),"spp":list(al.SPCD.value_counts().index[:3])})
    if len(recs["A"])<10 or len(recs["B"])<10:
        print(var,"fold too small A=%d B=%d"%(len(recs["A"]),len(recs["B"]))); sys.stdout.flush(); continue
    # density-dependent recruitment fit on fold A (item 1)
    def headroom(rec):
        return max(0.0,1.0-rec["sdi_t1"]/rec["sdi"]) if rec.get("sdi") and rec["sdi"]>0 and rec["sdi_t1"]>0 else 0.0
    _num=sum(rec["ing_tpa"] for rec in recs["A"] if headroom(rec)>0)
    _den=sum(headroom(rec)*rec["yrs"] for rec in recs["A"] if headroom(rec)>0)
    Rmax=_num/_den if _den>0 else 0.0
    def recruit_tpa(rec): return Rmax*headroom(rec)*rec["yrs"]
    def eval_fold(fold,baimult,inject):
        A={m:{"f":[],"o":[]} for m in ("BA","TPH","QMD","VOL")}
        for rec in recs[fold]:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":rec["statecd"],"COUNTYCD":0},str(rec["t1"]),var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
            tdf=build_fvs_treeinit(rec["al"],str(rec["t1"]))
            if baimult is None and not inject:
                s=run(std,tdf,str(rec["t1"]),"",rec["yrs"],var)
            else:
                kw=""
                if rec["sdi"] and rec["sdi"]>0: kw+=sdimax_kw(rec["sdi"])+"\n"
                if baimult is not None: kw+=allsp("BAIMULT",baimult)+"\n"
                tdf2=seed_rows(tdf,recruit_tpa(rec)*FRAC,rec["spp"]) if inject else tdf
                s=run(std,tdf2,str(rec["t1"]),kw,rec["yrs"],var)
            if s is None or len(s)==0: continue
            l=s.iloc[-1]; ba=g(l,"BA")*M2HA
            if ba<=0: continue
            A["BA"]["f"].append(ba); A["BA"]["o"].append(rec["oBA"])
            A["TPH"]["f"].append(g(l,"Tpa")*TPHc); A["TPH"]["o"].append(rec["oTPH"])
            A["QMD"]["f"].append(g(l,"QMD")*CMc); A["QMD"]["o"].append(rec["oQMD"])
            A["VOL"]["f"].append(g(l,"MCuFt")*VOLc); A["VOL"]["o"].append(rec["oVOL"])
        return A
    defA=eval_fold("A",None,False)
    injectA = bias(defA["TPH"]["f"],defA["TPH"]["o"]) < -2
    best_bai=None; best=1e9
    for bm in [1.0,0.9,0.8,0.7]:
        Ac=eval_fold("A",bm,injectA); q=abs(bias(Ac["QMD"]["f"],Ac["QMD"]["o"]))
        if q<best: best=q; best_bai=bm
    defB=eval_fold("B",None,False); calB=eval_fold("B",best_bai,injectA)   # OOS: arm A and arm B on held-out fold
    row={"variant":var,"nA":len(recs["A"]),"nB":len(recs["B"]),"Rmax_tpa_ac_yr":round(Rmax,3),"baimult":best_bai,"inject":int(injectA)}
    for m in ("BA","TPH","QMD","VOL"):
        row["A_%s"%m]=round(bias(defB[m]["f"],defB[m]["o"]),1)
        row["B_%s"%m]=round(bias(calB[m]["f"],calB[m]["o"]),1)
        row["A_%s_ci"%m]="[%+.1f,%+.1f]"%boot_ci(defB[m]["f"],defB[m]["o"])
        row["B_%s_ci"%m]="[%+.1f,%+.1f]"%boot_ci(calB[m]["f"],calB[m]["o"])
    rows.append(row)
    print("%-4s nB=%-3d BAIMULT=%.2f inj=%d | OOS A->B  BA %+.1f>%+.1f  TPH %+.1f>%+.1f  QMD %+.1f>%+.1f  VOL %+.1f>%+.1f"%(
        var,row["nB"],best_bai,injectA,row["A_BA"],row["B_BA"],row["A_TPH"],row["B_TPH"],row["A_QMD"],row["B_QMD"],row["A_VOL"],row["B_VOL"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
    print("MEDIAN |bias| arm A (default) -> arm B (calibrated), OOS:")
    for m in ("BA","TPH","QMD","VOL"):
        da=np.median([abs(r["A_%s"%m]) for r in rows]); cb=np.median([abs(r["B_%s"%m]) for r in rows]); print("  %-3s  %.1f -> %.1f"%(m,da,cb))
print("DONE_FOURARM")
