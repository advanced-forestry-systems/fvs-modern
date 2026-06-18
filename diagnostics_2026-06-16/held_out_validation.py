#!/usr/bin/env python3
"""Held-out spatial-fold validation of the fvs-modern keyword calibration (2026-06-17).
Addresses the red-team's #1 issue (in-sample leakage). Split COND-undisturbed plots by COUNTY into two
spatial folds (hash of county -> fold A/B). Derive the calibration on fold A only: ingrowth rate =
fold-A observed ingrowth %/decade; BAIMULT = the level (1.0/0.9/0.8/0.7) minimizing |QMD bias| on fold A.
brms SDImax is plot-level (not fold-derived). Apply the fold-A-derived calibration to held-out fold B.
Report default vs calibrated bias on BOTH folds: in-sample (A) vs out-of-sample (B). Env VARS,NSAMP,SEED."""
import os,sys,math,tempfile,sqlite3,shutil,csv,hashlib
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108; SDIc=2.471; TPAac=1/2.4710538
NSAMP=int(os.environ.get("NSAMP","220")); SEEDR=int(os.environ.get("SEED","5")); FRAC=0.7
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/held_out.csv"))
FIPS={"AL":1,"CA":6,"CT":9,"GA":13,"ID":16,"ME":23,"MA":25,"MS":28,"MT":30,"NH":33,"NY":36,"OR":41,"SC":45,"VT":50,"WA":53,"MI":26,"MN":27,"WI":55}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"sn":["AL","GA","SC","MS"],"kt":["MT","ID"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","ne,sn,kt,pn").split(",")
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["key"]=brms.STATECD.astype(str)+"-"+brms.UNITCD.astype(str)+"-"+brms.COUNTYCD.astype(str)+"-"+brms.PLOT.astype(str)
BRMS={k:v/SDIc for k,v in zip(brms.key,brms["SDImax.median"])}
def fold_of(statecd,countycd):  # spatial block: hash of state-county -> A/B
    h=int(hashlib.md5(("%d-%d"%(statecd,countycd)).encode()).hexdigest(),16)
    return "A" if h%2==0 else "B"
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
def sdimax_kw(val): return "\n".join("%-10s%10d%10.1f"%("SDIMAX",i+1,val) for i in range(MAXSP))
def stand_sdi(al):  # Reineke summation SDI, English per acre: sum TPA*(DIA/10)^1.605
    d=al[(al.DIA>0)&(al.TPA_UNADJ>0)]
    return float((d.TPA_UNADJ*(d.DIA/10.0)**1.605).sum()) if len(d) else 0.0
def seed_rows(tdf,rt,spp):
    if rt<=0 or not spp: return tdf
    base=int(tdf.tree_id.max())+1; per=rt/len(spp); rows=[{"stand_id":tdf.stand_id.iloc[0],"plot_id":1,"tree_id":base+j,"tree_count":per,"species":int(sp),"diameter":1.0,"ht":13.0,"crratio":40} for j,sp in enumerate(spp)]
    return pd.concat([tdf,pd.DataFrame(rows)],ignore_index=True)
def load_meta(states):
    cls={}; cnty={}
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
    return cls,cnty
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
def boot_ci(f,o,nb=2000,seed=11):  # percentile bootstrap 95% CI on the aggregate % bias, resampling plots
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
    fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states]
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    cls,cnty=load_meta(states)
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    # build per-plot records with fold + observed ingrowth
    recs={"A":[],"B":[]}
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed" or cn not in cnty: continue
        a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==cn]; al=a[(a.STATUSCD==1)&(a.DIA>0)]; bl=b[(b.STATUSCD==1)&(b.DIA>0)]
        if len(al)<5 or len(bl)<3: continue
        oBA,oTPH,oQMD=metr(bl)
        if oBA<=0: continue
        sc,cc,pkey=cnty[cn]; fold=fold_of(sc,cc)
        t1set=set(zip(a.SUBP,a.TREE)); ing=bl[~bl.set_index(["SUBP","TREE"]).index.isin(t1set)]
        recs[fold].append({"t1":t1,"cn":cn,"al":al,"oBA":oBA,"oTPH":oTPH,"oQMD":oQMD,"yrs":int(r.interval),
                           "statecd":sc,"sdi":BRMS.get(pkey),"sdi_t1":stand_sdi(al),"ing_tpa":ing.TPA_UNADJ.sum(),"init_tpa":float(al.TPA_UNADJ.sum()),
                           "spp":list(al.SPCD.value_counts().index[:3]),"yr10frac":int(r.interval)/10.0})
    if len(recs["A"])<10 or len(recs["B"])<10: print(var,"fold too small A=%d B=%d"%(len(recs["A"]),len(recs["B"]))); sys.stdout.flush(); continue
    # derive calibration on fold A: ingrowth rate, and BAIMULT minimizing |QMD bias| on A
    rateA=100*np.mean([r["ing_tpa"] for r in recs["A"]])/max(np.mean([r["init_tpa"] for r in recs["A"]]),1e-9)/np.mean([r["yrs"] for r in recs["A"]])*10/100.0
    # DENSITY-DEPENDENT RECRUITMENT (replaces fixed rate*interval*initialTPA). headroom shuts recruitment
    # off as the stand approaches its brms site-specific max SDI; R_max is fit on fold A so that
    # recruits = R_max * headroom * interval reproduces fold-A observed ingrowth, then applied unchanged to B.
    def headroom(rec):
        return max(0.0,1.0-rec["sdi_t1"]/rec["sdi"]) if rec.get("sdi") and rec["sdi"]>0 and rec["sdi_t1"]>0 else 0.0
    _num=sum(rec["ing_tpa"] for rec in recs["A"] if headroom(rec)>0)
    _den=sum(headroom(rec)*rec["yrs"] for rec in recs["A"] if headroom(rec)>0)
    Rmax=_num/_den if _den>0 else 0.0   # recruits per acre per year at zero density (headroom=1)
    def recruit_tpa(rec): return Rmax*headroom(rec)*rec["yrs"]
    # QMD over-predicted? decide injection sign from fold A default TPH
    def run_arm(rec,baimult,inject):
        std=build_fvs_standinit({"INVYR":2000,"STATECD":rec["statecd"],"COUNTYCD":0},str(rec["t1"]),var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
        tdf=build_fvs_treeinit(rec["al"],str(rec["t1"]))
        kw=""
        if rec["sdi"] and rec["sdi"]>0: kw+=sdimax_kw(rec["sdi"])+"\n"
        if baimult is not None: kw+=allsp("BAIMULT",baimult)+"\n"
        tdf_c=seed_rows(tdf,recruit_tpa(rec)*FRAC,rec["spp"]) if inject else tdf
        s=run(std,tdf_c,str(rec["t1"]),kw,rec["yrs"],var); return s
    # default fold A TPH sign
    defA={m:{"f":[],"o":[]} for m in ("BA","TPH","QMD")}
    for rec in recs["A"]:
        s=run(build_fvs_standinit({"INVYR":2000,"STATECD":rec["statecd"],"COUNTYCD":0},str(rec["t1"]),var)|{},None,None,0,0) if False else None
    # simpler: compute default on both folds + sweep BAIMULT on A
    def eval_fold(fold,baimult,inject):
        A={m:{"f":[],"o":[]} for m in ("BA","TPH","QMD")}
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
            A["BA"]["f"].append(ba); A["BA"]["o"].append(rec["oBA"]); A["TPH"]["f"].append(g(l,"Tpa")*TPHc); A["TPH"]["o"].append(rec["oTPH"]); A["QMD"]["f"].append(g(l,"QMD")*CMc); A["QMD"]["o"].append(rec["oQMD"])
        return A
    defA=eval_fold("A",None,False)
    tphA_def=bias(defA["TPH"]["f"],defA["TPH"]["o"]); injectA = tphA_def < -2
    best_bai=None; best=1e9
    for bm in [1.0,0.9,0.8,0.7]:
        Ac=eval_fold("A",bm,injectA); q=abs(bias(Ac["QMD"]["f"],Ac["QMD"]["o"]))
        if q<best: best=q; best_bai=bm
    # apply fold-A calibration (best_bai, rateA, injectA) to held-out fold B
    defB=eval_fold("B",None,False); calB=eval_fold("B",best_bai,injectA); calA=eval_fold("A",best_bai,injectA)
    row={"variant":var,"nA":len(recs["A"]),"nB":len(recs["B"]),"ingrowth_rateA_pct_dec":round(rateA*100,1),"Rmax_tpa_ac_yr":round(Rmax,4),"baimultA":best_bai,"injectA":int(injectA),
         "INSAMPLE_A_BA":"%+.1f>%+.1f"%(bias(defA["BA"]["f"],defA["BA"]["o"]),bias(calA["BA"]["f"],calA["BA"]["o"])),
         "INSAMPLE_A_QMD":"%+.1f>%+.1f"%(bias(defA["QMD"]["f"],defA["QMD"]["o"]),bias(calA["QMD"]["f"],calA["QMD"]["o"])),
         "INSAMPLE_A_TPH":"%+.1f>%+.1f"%(bias(defA["TPH"]["f"],defA["TPH"]["o"]),bias(calA["TPH"]["f"],calA["TPH"]["o"])),
         "OOS_B_BA":"%+.1f>%+.1f"%(bias(defB["BA"]["f"],defB["BA"]["o"]),bias(calB["BA"]["f"],calB["BA"]["o"])),
         "OOS_B_QMD":"%+.1f>%+.1f"%(bias(defB["QMD"]["f"],defB["QMD"]["o"]),bias(calB["QMD"]["f"],calB["QMD"]["o"])),
         "OOS_B_TPH":"%+.1f>%+.1f"%(bias(defB["TPH"]["f"],defB["TPH"]["o"]),bias(calB["TPH"]["f"],calB["TPH"]["o"]))}
    # bootstrap 95% CIs on the CALIBRATED bias (the number that must transfer), both folds, all metrics
    row.update({
        "INSAMPLE_A_BA_ci95":"[%+.1f,%+.1f]"%boot_ci(calA["BA"]["f"],calA["BA"]["o"]),
        "INSAMPLE_A_QMD_ci95":"[%+.1f,%+.1f]"%boot_ci(calA["QMD"]["f"],calA["QMD"]["o"]),
        "INSAMPLE_A_TPH_ci95":"[%+.1f,%+.1f]"%boot_ci(calA["TPH"]["f"],calA["TPH"]["o"]),
        "OOS_B_BA_ci95":"[%+.1f,%+.1f]"%boot_ci(calB["BA"]["f"],calB["BA"]["o"]),
        "OOS_B_QMD_ci95":"[%+.1f,%+.1f]"%boot_ci(calB["QMD"]["f"],calB["QMD"]["o"]),
        "OOS_B_TPH_ci95":"[%+.1f,%+.1f]"%boot_ci(calB["TPH"]["f"],calB["TPH"]["o"])})
    rows.append(row)
    print("%-4s nA=%d nB=%d Rmax=%.3f tpa/ac/yr BAIMULT=%.2f inj=%d | OOS(B) QMD %s %s BA %s %s TPH %s %s"%(
        var,row["nA"],row["nB"],Rmax,best_bai,injectA,row["OOS_B_QMD"],row["OOS_B_QMD_ci95"],row["OOS_B_BA"],row["OOS_B_BA_ci95"],row["OOS_B_TPH"],row["OOS_B_TPH_ci95"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_HELDOUT")
