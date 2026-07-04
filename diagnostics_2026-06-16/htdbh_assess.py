#!/usr/bin/env python3
"""Tree-level HT-DBH validation + recalibration factors (2026-06-17).
Blank input heights so FVS imputes via its HT-DBH curve; compare imputed height to FIA MEASURED height
(ACTUALHT) by DBH class and species per variant. TopHt can be unbiased while the curve is biased across
sizes. Outputs DBH-class bias and per-species recalibration ratio (measured/imputed). Env VARS,NSAMP,SEED."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
FTm=0.3048; CMc=2.54
NSAMP=int(os.environ.get("NSAMP","200")); SEEDR=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/htdbh.csv"))
OUTSP=os.environ.get("OUTSP",os.path.expanduser("~/overthin_work/htdbh_species.csv"))
FIPS={"AL":1,"CA":6,"CT":9,"GA":13,"ID":16,"ME":23,"MA":25,"MS":28,"MT":30,"NH":33,"NY":36,"OR":41,"SC":45,"VT":50,"WA":53}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["ME","NH","VT"],"sn":["AL","GA","SC"],"kt":["MT","ID"],"ie":["ID"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","ne,sn,kt,ie,pn").split(",")
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
def load_actualht(states):  # PLT_CN -> df(SUBP,TREE,DIA,ACTUALHT,SPCD) measured only
    fr=[]
    for s in states:
        f=Path(FIA)/f"{s}_TREE.csv"
        if not f.exists(): continue
        d=pd.read_csv(f,usecols=lambda c:c in ("PLT_CN","SUBP","TREE","STATUSCD","DIA","ACTUALHT","SPCD"),low_memory=False)
        d=d[(d.STATUSCD==1)&(d.DIA>0)&(d.ACTUALHT>0)]
        fr.append(d)
    return pd.concat(fr,ignore_index=True) if fr else None
def run_impute(std,tdf,sid,var):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=1,cycle_length=10))
    _run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+var),os.path.join(tmp,"t.key"),db,tmp)
    c=sqlite3.connect(db)
    try: tl=pd.read_sql_query("SELECT TreeId,SpeciesFIA,DBH,Ht FROM FVS_TreeList WHERE Year=(SELECT MIN(Year) FROM FVS_TreeList)",c)
    except: tl=None
    c.close(); shutil.rmtree(tmp,ignore_errors=True); return tl
DBHBINS=[(1,3),(3,5),(5,9),(9,13),(13,19),(19,40)]
agg={b:{"imp":[],"obs":[]} for b in DBHBINS}; spagg={}
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    plot=pd.concat([pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","STATECD"),low_memory=False) for s in states],ignore_index=True)
    cns=plot.CN.astype("int64").drop_duplicates().sample(n=min(NSAMP,len(plot)),random_state=SEEDR)
    aht=load_actualht(states)
    if aht is None: continue
    ahtg=dict(tuple(aht.groupby("PLT_CN")))
    tr=G.load_fia_trees(var,cns.tolist(),Path(FIA))
    vbin={b:{"imp":[],"obs":[]} for b in DBHBINS}; nplt=0
    for cn,grp in tr.groupby("PLT_CN"):
        al=grp[(grp.STATUSCD==1)&(grp.DIA>0)].copy()
        meas=ahtg.get(int(cn))
        if len(al)<5 or meas is None or len(meas)<3: continue
        sid=str(int(cn))
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(al.STATECD.iloc[0]),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
            tdf=build_fvs_treeinit(al,sid); tdf["ht"]=0.0   # blank -> impute
        except Exception: continue
        tl=run_impute(std,tdf,sid,var)
        if tl is None or len(tl)==0: continue
        nplt+=1
        # match measured to imputed by SPCD + nearest DBH within plot
        m=meas.copy(); m["DIAr"]=m.DIA.round(1)
        tl["DIAr"]=tl.DBH.round(1)
        for _,mr in m.iterrows():
            cand=tl[(tl.SpeciesFIA.astype(str)==str(int(mr.SPCD)).zfill(3))&(abs(tl.DIAr-mr.DIAr)<=0.2)]
            if len(cand)==0: cand=tl[abs(tl.DIAr-mr.DIAr)<=0.2]
            if len(cand)==0: continue
            imp=float(cand.Ht.iloc[0])*FTm; obs=float(mr.ACTUALHT)*FTm; dbh=float(mr.DIA)
            for b in DBHBINS:
                if b[0]<=dbh<b[1]: vbin[b]["imp"].append(imp); vbin[b]["obs"].append(obs); agg[b]["imp"].append(imp); agg[b]["obs"].append(obs); break
            sp=int(mr.SPCD); spagg.setdefault((var,sp),{"imp":[],"obs":[]}); spagg[(var,sp)]["imp"].append(imp); spagg[(var,sp)]["obs"].append(obs)
    for b in DBHBINS:
        o=vbin[b]["obs"]; i=vbin[b]["imp"]
        if len(o)<10: continue
        bias=100*(sum(i)-sum(o))/sum(o)
        rows.append({"variant":var,"dbh_class":"%d-%d in"%b,"n":len(o),"obs_ht_m":round(np.mean(o),1),"imp_ht_m":round(np.mean(i),1),"ht_bias%":round(bias,1)})
        print("%-4s %-8s n=%-4d obsHt %4.1f impHt %4.1f bias %+5.1f%%"%(var,"%d-%d"%b,len(o),np.mean(o),np.mean(i),bias)); sys.stdout.flush()
    print("--- %s done (n plots %d) ---"%(var,nplt)); sys.stdout.flush()
print("=== ALL-VARIANT pooled by DBH class ===")
for b in DBHBINS:
    o=agg[b]["obs"]; i=agg[b]["imp"]
    if len(o)<10: continue
    print("DBH %-7s n=%-5d obsHt %4.1f impHt %4.1f bias %+5.1f%%"%("%d-%d"%b,len(o),np.mean(o),np.mean(i),100*(sum(i)-sum(o))/sum(o)))
# per-species recalibration ratio (obs/imp)
sprows=[]
for (var,sp),d in spagg.items():
    if len(d["obs"])<20: continue
    ratio=sum(d["obs"])/sum(d["imp"]) if sum(d["imp"])>0 else float("nan")
    sprows.append({"variant":var,"SPCD":sp,"n":len(d["obs"]),"recal_ratio_obs_over_imp":round(ratio,3)})
import csv as _csv
if rows:
    with open(OUT,"w",newline="") as fh: w=_csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
if sprows:
    with open(OUTSP,"w",newline="") as fh: w=_csv.DictWriter(fh,fieldnames=list(sprows[0].keys())); w.writeheader(); w.writerows(sprows)
    print("species ratios written:",len(sprows))
print("DONE_HTDBH")
