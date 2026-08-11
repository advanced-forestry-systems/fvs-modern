#!/usr/bin/env python3
"""Volume + height (HT-DBH) assessment and REGHMULT lever test (2026-06-17).
On COND-undisturbed plots per variant: observed t2 volume (sum VOLCFNET*TPA -> m3/ha), biomass
(sum DRYBIO_AG*TPA -> Mg/ha), top height (TPA-wtd mean ACTUALHT of dominant trees). FVS volume = TCuFt,
height = TopHt. Arms: default, REGHMULT 0.90, 0.80 (height-growth lever). Reports volume/height bias.
Env VARS,NSAMP,SEED."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; MAXSP=108
VOLc=0.0699055   # cuft/acre -> m3/ha
BIOc=0.00112085  # lb/acre -> Mg/ha (lb->kg 0.4536, /acre->/ha 2.4711, /1000) = 0.4536*2.4711/1000
FTc=0.3048
NSAMP=int(os.environ.get("NSAMP","130")); SEEDR=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/assess_vol.csv"))
FIPS={"AL":1,"CA":6,"CT":9,"GA":13,"ID":16,"ME":23,"MA":25,"MS":28,"MT":30,"NH":33,"NY":36,"OR":41,"SC":45,"VT":50,"WA":53,"MI":26,"MN":27,"WI":55}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"sn":["AL","GA","SC","MS"],"kt":["MT","ID"],"ie":["ID","MT"],"pn":["OR","WA"],"ls":["MI","MN","WI"]}
VARS=os.environ.get("VARS","ne,sn,kt,pn").split(",")
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
def allsp(kw,val): return "\n".join("%-16s%10d%10.4f"%(kw,i+1,val) for i in range(MAXSP))
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
def load_obs(states):  # PLT_CN -> dict of observed vol/bio/topht/ba
    frames=[]
    for s in states:
        f=Path(FIA)/f"{s}_TREE.csv"
        if not f.exists(): continue
        d=pd.read_csv(f,usecols=lambda c:c in ("PLT_CN","STATUSCD","DIA","TPA_UNADJ","VOLCFNET","DRYBIO_AG","ACTUALHT"),low_memory=False)
        frames.append(d)
    t=pd.concat(frames,ignore_index=True); t=t[(t.STATUSCD==1)&(t.DIA>0)&(t.TPA_UNADJ>0)]
    return t
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
    fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states]
    plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
    rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
    rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
    rem=rem[(rem.interval>=5)&(rem.interval<=15)].sample(n=min(NSAMP,len(rem)),random_state=SEEDR)
    cls=load_cond(states); obs=load_obs(states); obsg=dict(tuple(obs.groupby("PLT_CN")))
    tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA))
    ARMS=["default","REGH0.90","REGH0.80"]
    A={k:{m:{"f":[],"o":[]} for m in ("VOL","HT")} for k in ARMS}; nu=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); cn=int(r.CN)
        if cls.get(cn)!="undisturbed": continue
        al=tr1[(tr1.PLT_CN==t1)&(tr1.STATUSCD==1)&(tr1.DIA>0)]
        b=obsg.get(cn)
        if len(al)<5 or b is None or len(b)<3: continue
        # observed t2 volume (m3/ha), top height (m): dominant = largest 40 tpa
        oVOL=(b.VOLCFNET.fillna(0)*b.TPA_UNADJ).sum()*VOLc
        bb=b.sort_values("DIA",ascending=False).copy(); bb["cum"]=bb.TPA_UNADJ.cumsum()
        dom=bb[bb.cum<=40] if (bb.cum<=40).any() else bb.head(3)
        oHT=(dom.ACTUALHT.fillna(0)*dom.TPA_UNADJ).sum()/max(dom.TPA_UNADJ.sum(),1e-9)*FTc
        if oVOL<=0 or oHT<=0: continue
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(al,sid)
        except Exception: continue
        res={}
        res["default"]=run(std,tdf,sid,"",yrs,var)
        res["REGH0.90"]=run(std,tdf,sid,allsp("REGHMULT",0.90),yrs,var)
        res["REGH0.80"]=run(std,tdf,sid,allsp("REGHMULT",0.80),yrs,var)
        if any(res[k] is None or len(res[k])==0 for k in ARMS): continue
        nu+=1
        for k in ARMS:
            l=res[k].iloc[-1]
            A[k]["VOL"]["f"].append(g(l,"TCuFt")*VOLc); A[k]["VOL"]["o"].append(oVOL)
            A[k]["HT"]["f"].append(g(l,"TopHt")*FTc); A[k]["HT"]["o"].append(oHT)
    if nu<8: print(var,"too few (n=%d)"%nu); sys.stdout.flush(); continue
    row={"variant":var,"n":nu}
    for k in ARMS: row[k+"_VOL%"]=round(bias(A[k]["VOL"]["f"],A[k]["VOL"]["o"]),1); row[k+"_HT%"]=round(bias(A[k]["HT"]["f"],A[k]["HT"]["o"]),1)
    rows.append(row)
    print("%-4s n=%-3d VOL def %+6.1f R90 %+6.1f R80 %+6.1f | HT def %+6.1f R90 %+6.1f R80 %+6.1f"%(
        var,nu,row["default_VOL%"],row["REGH0.90_VOL%"],row["REGH0.80_VOL%"],row["default_HT%"],row["REGH0.90_HT%"],row["REGH0.80_HT%"])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_VOL")
