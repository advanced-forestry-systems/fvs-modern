#!/usr/bin/env python3
"""Volume sensitivity to HT-DBH: FVS-imputed heights vs FIA-measured heights (2026-06-17).
Same stands, two runs: (a) blank heights -> FVS imputes via HT-DBH curve; (b) FIA measured heights.
Volume (MCuFt) difference isolates the volume impact of the HT-DBH bias. Env VARS,NSAMP,SEED."""
import os,sys,math,tempfile,sqlite3,shutil,csv
P="/users/PUOM0008/crsfaaron/fvs-modern"; FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; CONUS="/users/PUOM0008/crsfaaron/fvs-conus"
os.environ.update(FIA_DATA_DIR=FIA,FVS_PROJECT_ROOT=P,FVS_LIB_DIR=P+"/lib",FVS_CONFIG_DIR=P+"/config")
for p in [os.path.expanduser("~/overthin_work"),CONUS+"/python",P]: sys.path.insert(0,p)
import pandas as pd, numpy as np
from pathlib import Path
import fia_stand_generator as G
from perseus_100yr_projection import build_fvs_standinit, build_fvs_treeinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess
VOLc=0.0699055
NSAMP=int(os.environ.get("NSAMP","150")); SEEDR=int(os.environ.get("SEED","5"))
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/volht.csv"))
FIPS={"AL":1,"GA":13,"ID":16,"ME":23,"MS":28,"MT":30,"NH":33,"VT":50,"OR":41,"WA":53,"SC":45}
INV={v:k for k,v in FIPS.items()}
VARMAP={"ne":["ME","NH","VT"],"sn":["AL","GA","SC"],"kt":["MT","ID"],"pn":["OR","WA"]}
VARS=os.environ.get("VARS","ne,sn,kt,pn").split(",")
_o=G._state_abbrev; G._state_abbrev=lambda c:INV.get(c) or _o(c)
def run(std,tdf,sid,var):
    tmp=tempfile.mkdtemp(); db=os.path.join(tmp,"FVS_Data.db")
    c=sqlite3.connect(db); std.to_sql("fvs_standinit",c,if_exists="replace",index=False); tdf.to_sql("fvs_treeinit",c,if_exists="replace",index=False); c.close()
    open(os.path.join(tmp,"t.key"),"w").write(KEYFILE_TEMPLATE.format(stand_id=sid,db_path=db,calibration_keywords="** D",num_cycles=1,cycle_length=10))
    r=_run_via_subprocess(os.path.join(FVS_LIB_DIR,"FVS"+var),os.path.join(tmp,"t.key"),db,tmp); shutil.rmtree(tmp,ignore_errors=True); return r.get("summary")
def g(r,k):
    try: return float(r.get(k,0) or 0)
    except: return 0.0
rows=[]
for var in VARS:
    states=VARMAP.get(var)
    if not states: continue
    G.VARIANT_STATES[var]=tuple(FIPS[s] for s in states)
    plot=pd.concat([pd.read_csv(Path(FIA)/f"{s}_PLOT.csv",usecols=lambda c:c in ("CN","STATECD"),low_memory=False) for s in states],ignore_index=True)
    cns=plot.CN.astype("int64").drop_duplicates().sample(n=min(NSAMP,len(plot)),random_state=SEEDR)
    tr=G.load_fia_trees(var,cns.tolist(),Path(FIA))
    imp=[]; mea=[]; nu=0
    for cn,grp in tr.groupby("PLT_CN"):
        al=grp[(grp.STATUSCD==1)&(grp.DIA>0)&(grp.HT>0)].copy()
        if len(al)<8: continue
        sid=str(int(cn))
        try:
            std=build_fvs_standinit({"INVYR":2000,"STATECD":int(al.STATECD.iloc[0]),"COUNTYCD":0},sid,var); std["inv_plot_size"]=1.0; std["brk_dbh"]=99.0
            tdf_m=build_fvs_treeinit(al,sid)              # measured heights
            tdf_i=tdf_m.copy(); tdf_i["ht"]=0.0           # blank -> FVS imputes
        except Exception: continue
        sm=run(std,tdf_m,sid,var); si=run(std,tdf_i,sid,var)
        if sm is None or len(sm)==0 or si is None or len(si)==0: continue
        vm=g(sm.iloc[0],"MCuFt")*VOLc; vi=g(si.iloc[0],"MCuFt")*VOLc   # year-0 volume
        if vm<=0: continue
        mea.append(vm); imp.append(vi); nu+=1
    if nu<8: print(var,"too few (n=%d)"%nu); sys.stdout.flush(); continue
    mm=np.mean(mea); mi=np.mean(imp); d=100*(mi-mm)/mm
    rows.append({"variant":var,"n":nu,"vol_measured_ht":round(mm,1),"vol_imputed_ht":round(mi,1),"imputed_vs_measured%":round(d,1)})
    print("%-4s n=%-3d  vol(measured Ht) %5.1f  vol(FVS-imputed Ht) %5.1f m3/ha  ->  imputed is %+5.1f%%"%(var,nu,mm,mi,d)); sys.stdout.flush()
if rows:
    with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote",OUT)
print("DONE_VOLHT")
