#!/usr/bin/env python3
"""Disturbance-stratified FIA benchmark of default FVS across variants (2026-06-17).
For each variant, sample FIA remeasurement pairs, run the default engine 1 cycle, and split the BA bias
by what actually happened in the real plot: grew (undisturbed), flat/declined (natural mortality), or
harvested (>15% of t1 live stems cut/removed). Shows how much of the pooled "over-prediction" is a
disturbance/benchmark artifact vs a real undisturbed-projection bias. Env NSAMP, SEED, OUT, VARS."""
import os,sys,math,tempfile,sqlite3,shutil,json,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
M2HA=0.2296; TPHc=2.4710538
NSAMP=int(os.environ.get("NSAMP","130")); SEED=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/multivar_disturb.csv"))
FIPS={"AL":1,"AZ":4,"AR":5,"CA":6,"CO":8,"CT":9,"DE":10,"FL":12,"GA":13,"ID":16,"IL":17,"IN":18,"IA":19,"KS":20,"KY":21,"LA":22,"ME":23,"MD":24,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NE":31,"NV":32,"NH":33,"NJ":34,"NM":35,"NY":36,"NC":37,"ND":38,"OH":39,"OK":40,"OR":41,"PA":42,"RI":44,"SC":45,"SD":46,"TN":47,"TX":48,"UT":49,"VT":50,"VA":51,"WA":53,"WV":54,"WI":55,"WY":56}
INV={v:k for k,v in FIPS.items()}
VARMAP={
 # Eastern / central / southern (state lists trimmed for I/O; enough plots for a bias estimate)
 "ne":["CT","ME","MA","NH","NY","VT"],
 "acd":["ME","NH","VT"],
 "sn":["AL","GA","SC","MS"],
 "ls":["MI","MN","WI"],
 "cs":["IL","IN","MO"],
 # Intermountain / Rockies
 "ie":["ID","MT","WA"],
 "kt":["MT","ID"],
 "ci":["ID"],
 "em":["MT","ND","SD"],
 "bm":["OR","WA"],
 "cr":["CO","WY"],
 "tt":["WY","ID"],
 "ut":["UT","NV"],
 # Pacific / California
 "ca":["CA"],
 "ws":["CA"],
 "nc":["CA","OR"],
 "so":["CA","OR"],
 "ec":["OR","WA"],
 "wc":["OR","WA"],
 "oc":["OR","WA"],
 "on":["OR","WA"],
 "op":["WA"],
 "pn":["OR","WA"],
 # No FIA coverage (will auto-skip): Alaska, British Columbia
 "ak":["AK"],
 "bc":["AK"],
}
VARS=os.environ.get("VARS","").split(",") if os.environ.get("VARS") else list(VARMAP)
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
def ba(df):
    d=df[(df.DIA>0)&(df.TPA_UNADJ>0)]
    return 0 if len(d)==0 else (d.TPA_UNADJ*0.005454*d.DIA**2).sum()*M2HA
def run(std,tdf,sid,yrs,var):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=1,cycle_length=yrs))
    lib=os.path.join(FVS_LIB_DIR,"FVS"+var)
    if not os.path.exists(lib): shutil.rmtree(tmp,ignore_errors=True); return None
    r=_run_via_subprocess(lib,os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def bias(f,o):
    m=[(x,y) for x,y in zip(f,o) if y>0]; return (float("nan"),0) if not m else (100*sum(x-y for x,y in m)/sum(y for _,y in m),len(m))
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    try:
        fr=[pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","PREV_PLT_CN","MEASYEAR","STATECD"),low_memory=False) for s in states if (Path(FIA)/f"{s}_PLOT.csv").exists()]
        plot=pd.concat(fr,ignore_index=True); yr=dict(zip(plot.CN.astype("int64"),plot.MEASYEAR))
        rem=plot.dropna(subset=["PREV_PLT_CN"]).copy(); rem["PREV_PLT_CN"]=rem.PREV_PLT_CN.astype("int64"); rem["CN"]=rem.CN.astype("int64")
        rem["interval"]=rem.apply(lambda r:r.MEASYEAR-yr.get(r.PREV_PLT_CN,np.nan),axis=1)
        rem=rem[(rem.interval>=5)&(rem.interval<=15)]
        if len(rem)==0: print(var,"no remeasurement"); sys.stdout.flush(); continue
        rem=rem.sample(n=min(NSAMP,len(rem)),random_state=SEED)
        tr1=G.load_fia_trees(var,rem.PREV_PLT_CN.tolist(),Path(FIA)); tr2=G.load_fia_trees(var,rem.CN.tolist(),Path(FIA))
    except Exception as e:
        print(var,"LOAD ERR",repr(e)[:120]); sys.stdout.flush(); continue
    grp={"grew":{"f":[],"o":[]},"flat_lost":{"f":[],"o":[]},"harvested":{"f":[],"o":[]}}
    nrec=0
    for _,r in rem.iterrows():
        t1=int(r.PREV_PLT_CN); a=tr1[tr1.PLT_CN==t1]; b=tr2[tr2.PLT_CN==int(r.CN)]
        if len(a)<5 or len(b)<3: continue
        oBA1=ba(a); oBA2=ba(b)
        if oBA1<=0 or oBA2<=0: continue
        al=a[(a.STATUSCD==1)&(a.DIA>0)]; m=pd.merge(al[["SUBP","TREE"]],b[["SUBP","TREE","STATUSCD"]],on=["SUBP","TREE"],how="left")
        cut=int((m.STATUSCD==3).sum()); gone=int(m.STATUSCD.isna().sum()); harv=(cut+gone)/max(len(al),1)>0.15
        sid=str(t1); yrs=int(r.interval)
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(r.STATECD),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0; tdf=build_fvs_treeinit(a,sid)
        except Exception: continue
        s=run(std,tdf,sid,yrs,var)
        if s is None or len(s)==0: continue
        fbaN=g(s.iloc[-1],"BA")*M2HA
        if fbaN<=0: continue
        nrec+=1
        k="harvested" if harv else ("grew" if oBA2>oBA1*1.02 else "flat_lost")
        grp[k]["f"].append(fbaN); grp[k]["o"].append(oBA2)
    allf=sum((grp[k]["f"] for k in grp),[]); allo=sum((grp[k]["o"] for k in grp),[])
    pooled=bias(allf,allo); grew=bias(grp["grew"]["f"],grp["grew"]["o"]); flat=bias(grp["flat_lost"]["f"],grp["flat_lost"]["o"]); harvb=bias(grp["harvested"]["f"],grp["harvested"]["o"])
    rows.append({"variant":var,"n":nrec,"pooled_BA%":round(pooled[0],1),"undist_BA%":round(grew[0],1),"undist_n":grew[1],"flatlost_BA%":round(flat[0],1),"flatlost_n":flat[1],"harvest_BA%":round(harvb[0],1),"harvest_n":harvb[1]})
    print("%-4s n=%-4d pooled %+6.1f  undist %+6.1f(n=%d)  flat %+6.1f(n=%d)  harv %+6.1f(n=%d)"%(var,nrec,pooled[0],grew[0],grew[1],flat[0],flat[1],harvb[0],harvb[1])); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_MULTI")
