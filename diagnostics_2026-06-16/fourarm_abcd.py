#!/usr/bin/env python3
"""Four-arm A/B/C/D comparison in ONE framework (FVS engine), all variants, disturbance-clean OOS (2026-06-18).

  A = default FVS
  B = fvs-modern keyword layer: brms maxSDI + density-dependent recruitment + fold-A global BAIMULT
  C = fvs-conus growth signal in-engine: per-species BAIMULT = observed DG / default-FVS DG ({var}_baimult_calib.json,
      built by calib_ne.py; the fvs-conus DG equations are calibrated to that same observed growth, so the
      per-species multiplier emulates the fvs-conus growth equations as engine keywords)
  D = combined: per-species BAIMULT (C) + brms maxSDI + density-dependent recruitment (B's density levers)

All four arms run in the FVS engine, so the projector-vs-engine sign problem does not arise. All four
metrics (BA/TPH/QMD/merch volume) with percentile-bootstrap 95% CIs, reported out-of-sample on spatial
fold B (calibration derived on fold A only). Arm C is an emulation of the fvs-conus growth equations via
their per-species DG signal; it captures the diameter-growth level effect, not the full trait-driven /
height-diameter / annualized structure, which still awaits the fvs2py in-engine injection. Env VARS,NSAMP,SEED,OUT.
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
NSAMP=int(os.environ.get("NSAMP","600")); SEEDR=int(os.environ.get("SEED","5")); FRAC=0.7
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/fourarm_abcd.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"ND":38,"NV":32,"NH":33,"NY":36,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56,"AK":2}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"acd":["ME","NH","VT"],"sn":["AL","GA","SC","MS"],"ls":["MI","MN","WI"],
        "cs":["IL","IN","MO"],"ie":["ID","MT","WA"],"kt":["MT","ID"],"ci":["ID"],"em":["MT","ND","SD"],"bm":["OR","WA"],
        "cr":["CO","WY"],"tt":["WY","ID"],"ut":["UT","NV"],"ca":["CA"],"ws":["CA"],"nc":["CA","OR"],"so":["CA","OR"],
        "ec":["OR","WA"],"wc":["OR","WA"],"oc":["OR","WA"],"op":["WA"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","").split(",") if os.environ.get("VARS") else list(VARMAP)
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["key"]=brms.STATECD.astype(str)+"-"+brms.UNITCD.astype(str)+"-"+brms.COUNTYCD.astype(str)+"-"+brms.PLOT.astype(str)
BRMS={k:v/SDIc for k,v in zip(brms.key,brms["SDImax.median"])}
def fold_of(sc,cc):
    return "A" if int(hashlib.md5(("%d-%d"%(sc,cc)).encode()).hexdigest(),16)%2==0 else "B"
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
def sdimax_kw(val): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,val) for i in range(MAXSP))
def persp_kw(d): return "\n".join("%-16s%10d%10.4f"%("BAIMULT",int(i),float(v)) for i,v in d.items()) if d else ""
def stand_sdi(al):
    x=al[(al.DIA>0)&(al.TPA_UNADJ>0)]; return float((x.TPA_UNADJ*(x.DIA/10.0)**1.605).sum()) if len(x) else 0.0
def seed_rows(tdf,rt,spp):
    if rt<=0 or not spp: return tdf
    base=int(tdf.tree_id.max())+1; per=rt/len(spp)
    r=[{"stand_id":tdf.stand_id.iloc[0],"plot_id":1,"tree_id":base+j,"tree_count":per,"species":int(sp),"diameter":1.0,"ht":13.0,"crratio":40} for j,sp in enumerate(spp)]
    return pd.concat([tdf,pd.DataFrame(r)],ignore_index=True)
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
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states if s in FIPS)
    bj=os.path.expanduser("~/overthin_work/%s_baimult_calib.json"%var)
    PERSP=persp_kw(json.load(open(bj))) if os.path.exists(bj) else ""
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
    def headroom(rec): return max(0.0,1.0-rec["sdi_t1"]/rec["sdi"]) if rec.get("sdi") and rec["sdi"]>0 and rec["sdi_t1"]>0 else 0.0
    _num=sum(rec["ing_tpa"] for rec in recs["A"] if headroom(rec)>0); _den=sum(headroom(rec)*rec["yrs"] for rec in recs["A"] if headroom(rec)>0)
    Rmax=_num/_den if _den>0 else 0.0
    def recruit_tpa(rec): return Rmax*headroom(rec)*rec["yrs"]
    def eval_fold(fold,sdimax_on,baimult_level,persp_on,inject):
        A={m:{"f":[],"o":[]} for m in ("BA","TPH","QMD","VOL")}
        for rec in recs[fold]:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":rec["statecd"],"COUNTYCD":0},str(rec["t1"]),var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
            tdf=build_fvs_treeinit(rec["al"],str(rec["t1"]))
            kw=""
            if sdimax_on and rec["sdi"] and rec["sdi"]>0: kw+=sdimax_kw(rec["sdi"])+"\n"
            if persp_on and PERSP: kw+=PERSP+"\n"
            elif baimult_level is not None: kw+=allsp("BAIMULT",baimult_level)+"\n"
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
    defA=eval_fold("A",False,None,False,False)
    injectA = bias(defA["TPH"]["f"],defA["TPH"]["o"]) < -2
    best_bai=None; best=1e9
    for bm in [1.0,0.9,0.8,0.7]:
        Ac=eval_fold("A",True,bm,False,injectA); q=abs(bias(Ac["QMD"]["f"],Ac["QMD"]["o"]))
        if q<best: best=q; best_bai=bm
    arms={"A":eval_fold("B",False,None,False,False),
          "B":eval_fold("B",True,best_bai,False,injectA),
          "C":eval_fold("B",False,None,True,False),
          "D":eval_fold("B",True,None,True,injectA)}
    row={"variant":var,"nB":len(recs["B"]),"Rmax":round(Rmax,3),"baimult":best_bai,"inject":int(injectA),"n_persp_spp":len(json.load(open(bj))) if os.path.exists(bj) else 0}
    for m in ("BA","TPH","QMD","VOL"):
        for arm in ("A","B","C","D"):
            row["%s_%s"%(arm,m)]=round(bias(arms[arm][m]["f"],arms[arm][m]["o"]),1)
        row["B_%s_ci"%m]="[%+.1f,%+.1f]"%boot_ci(arms["B"][m]["f"],arms["B"][m]["o"])
        row["D_%s_ci"%m]="[%+.1f,%+.1f]"%boot_ci(arms["D"][m]["f"],arms["D"][m]["o"])
    rows.append(row)
    print("%-4s nB=%-3d persp=%-3d | OOS QMD A%+.0f B%+.0f C%+.0f D%+.0f | VOL A%+.0f B%+.0f C%+.0f D%+.0f | BA A%+.0f D%+.0f | TPH A%+.0f D%+.0f"%(
        var,row["nB"],row["n_persp_spp"],row["A_QMD"],row["B_QMD"],row["C_QMD"],row["D_QMD"],row["A_VOL"],row["B_VOL"],row["C_VOL"],row["D_VOL"],row["A_BA"],row["D_BA"],row["A_TPH"],row["D_TPH"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
    print("MEDIAN |bias| across variants, arms A/B/C/D:")
    for m in ("BA","TPH","QMD","VOL"):
        vals={arm:np.median([abs(r["%s_%s"%(arm,m)]) for r in rows if r["%s_%s"%(arm,m)]==r["%s_%s"%(arm,m)]]) for arm in ("A","B","C","D")}
        print("  %-3s  A %.1f  B %.1f  C %.1f  D %.1f"%(m,vals["A"],vals["B"],vals["C"],vals["D"]))
print("DONE_ABCD")
