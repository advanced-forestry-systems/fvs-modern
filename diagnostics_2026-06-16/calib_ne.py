#!/usr/bin/env python3
"""Per-species diameter-growth calibration of an FVS variant from FIA (2026-06-16).
BAIMULT_species = mean observed DG / mean FVS-predicted DG, observed from matched FIA remeasurement,
predicted from FVS_TreeList.DG, species index from FVS_InvReference (SpeciesNum<->SpeciesFIA). Apply
per-species BAIMULT, validate stand BA/QMD/TPH bias default vs calibrated. Env NSAMP,SEED,VAR,STATES,MAXSP,OUT."""
import os,sys,math,tempfile,sqlite3,json,shutil
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
VAR=os.environ.get("VAR","ne"); STATES=os.environ.get("STATES","CT,ME,MA,NH,NY,RI,VT").split(","); MAXSP=int(os.environ.get("MAXSP","108"))
NSAMP=int(os.environ.get("NSAMP","150")); SEED=int(os.environ.get("SEED","5")); OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/%s_baimult_calib.json"%VAR.lower()))
M2HA=0.2296; TPHc=2.4710538; CMc=2.54
_REV={"AL":1,"CA":6,"CO":8,"CT":9,"FL":12,"GA":13,"ID":16,"IL":17,"IN":18,"IA":19,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NM":35,"NY":36,"OR":41,"RI":44,"SC":45,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
G.VARIANT_STATES[VAR]=tuple(_REV[x] for x in STATES if x in _REV); _o=G._state_abbrev; _A={v:k for k,v in _REV.items()}; G._state_abbrev=lambda c:_A.get(c) or _o(c)
SPCD2IDX={}
def runkw(std,tdf,sid,kw,yrs,want_tl=False):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords=kw or "** D",num_cycles=1,cycle_length=yrs))
    res=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+VAR),os.path.join(tmp,"t.key"),db,tmp)
    tl=None
    c=sqlite3.connect(db)
    try:
        ir=pd.read_sql_query("SELECT SpeciesNum,SpeciesFIA FROM FVS_InvReference",c)
        for _,row in ir.iterrows():
            try: SPCD2IDX[int(row.SpeciesFIA)]=int(row.SpeciesNum)
            except: pass
    except: pass
    if want_tl:
        try: tl=pd.read_sql_query("SELECT Year,SpeciesFIA,DG FROM FVS_TreeList",c)
        except: tl=None
    c.close(); shutil.rmtree(tmp,ignore_errors=True)
    return res.get("summary"), tl
def metr(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    if len(d)==0: return (0,0,0)
    return ((d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA, d.TPA_UNADJ.sum()*TPHc, math.sqrt((d.DIA**2*d.TPA_UNADJ).sum()/d.TPA_UNADJ.sum())*CMc)
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def stt(p,o):
    m=[(a,b) for a,b in zip(p,o) if a==a and b==b and b>0]; k=len(m)
    return (float("nan"),0) if k==0 else (100*sum(a-b for a,b in m)/sum(b for _,b in m),k)
def main():
    fr=[pd.read_csv(Path(FIA)/f"{ab}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for ab in STATES if (Path(FIA)/f"{ab}_PLOT.csv").exists()]
    plot=pd.concat(fr,ignore_index=True); yrm=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yrm.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEED)
    tr1=G.load_fia_trees(VAR,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(VAR,rem.CN.tolist(),Path(FIA))
    obs={}; pred={}; pairs=[]
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==int(r.CN)]
        if len(a)<5 or len(b)<3: continue
        yrs=int(r.interval)
        m=pd.merge(a[(a.STATUSCD==1)&(a.DIA>0)][["SUBP","TREE","SPCD","DIA"]], b[(b.STATUSCD==1)&(b.DIA>0)][["SUBP","TREE","DIA"]], on=["SUBP","TREE"], suffixes=("1","2"))
        m=m[m.DIA2>=m.DIA1]
        for spcd,grp in m.groupby("SPCD"): obs.setdefault(int(spcd),[]).extend(list(grp.DIA2-grp.DIA1))
        sid=str(t1); std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,VAR); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
        tdf=build_fvs_treeinit(a,sid); pairs.append((std,tdf,sid,yrs,metr(b)))
        s,tl=runkw(std,tdf,sid,"",yrs,want_tl=True)
        if tl is not None and len(tl):
            t0=tl[tl.Year==tl.Year.max()]
            for spc,grp in t0.groupby("SpeciesFIA"):
                try: sp=int(spc)
                except: continue
                pred.setdefault(sp,[]).extend([x for x in grp.DG if x==x and x>0])
    BAI={}
    print("SPCD2IDX size:",len(SPCD2IDX))
    print("spcd  n   obsDG  predDG  BAIMULT")
    for spcd in sorted(set(obs)&set(pred)&set(SPCD2IDX)):
        if len(obs[spcd])<10 or len(pred[spcd])<5: continue
        o=np.mean(obs[spcd]); p=np.mean(pred[spcd])
        if p<=0: continue
        b=min(max(o/p,0.2),3.0); BAI[SPCD2IDX[spcd]]=b
        if len(obs[spcd])>=40: print("%5d %4d  %5.2f  %6.2f  %6.3f"%(spcd,len(obs[spcd]),o,p,b))
    kw="\n".join("BAIMULT         %10d%10.4f"%(idx,BAI[idx]) for idx in sorted(BAI))
    json.dump({str(k):v for k,v in BAI.items()}, open(OUT,"w"))
    print("calibrated species:",len(BAI))
    A={a:{x:[] for x in ["BA","QMD","TPH","oBA","oTPH","oQMD"]} for a in ["default","calibrated"]}
    for std,tdf,sid,yrs,(oBA,oTPH,oQMD) in pairs:
        if oBA<=0: continue
        for lab,k in [("default",""),("calibrated",kw)]:
            s,_=runkw(std,tdf,sid,k,yrs)
            if s is None or len(s)==0: continue
            l=s.iloc[-1]; ba=g(l,"BA")*M2HA
            if ba<=0: continue
            A[lab]["BA"].append(ba); A[lab]["QMD"].append(g(l,"QMD")*CMc); A[lab]["TPH"].append(g(l,"Tpa")*TPHc)
            A[lab]["oBA"].append(oBA); A[lab]["oQMD"].append(oQMD); A[lab]["oTPH"].append(oTPH)
    print(f"\n{'arm':<12}{'BA bias%':>10}{'QMD bias%':>11}{'TPH bias%':>11}{'n':>6}")
    for lab in ["default","calibrated"]:
        b=stt(A[lab]["BA"],A[lab]["oBA"]); q=stt(A[lab]["QMD"],A[lab]["oQMD"]); t=stt(A[lab]["TPH"],A[lab]["oTPH"])
        print(f"{lab:<12}{b[0]:>10.1f}{q[0]:>11.1f}{t[0]:>11.1f}{b[1]:>6d}")
    print("DONE_CALIB")
if __name__=="__main__": main()
