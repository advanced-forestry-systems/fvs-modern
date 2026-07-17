import os,sys
import numpy as np, pandas as pd
PR=os.environ.get("FVS_PROJECT_ROOT",os.path.expanduser("~/fvs-modern"))
sys.path.insert(0,PR); sys.path.insert(0,os.path.join(PR,"calibration","python"))
sys.path.insert(0,"/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import perseus_100yr_projection as P
import run_gompit_projection as G
import run_conus_task_fvstreeinit as RC
FIPS=G.FIPS
SI="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant/standinit_NE.csv"
TI="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit"
si=pd.read_csv(SI,low_memory=False); si["STAND_CN"]=si["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
picked=None
for sf,grp in si.groupby("STATE"):
    try: state=FIPS[int(float(sf))]
    except: continue
    f=os.path.join(TI,f"{state}_FVS_TREEINIT_PLOT.csv")
    if not os.path.exists(f): continue
    tt=pd.read_csv(f,low_memory=False); tt["STAND_CN"]=tt["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
    by={k:v for k,v in tt.groupby("STAND_CN")}
    for _,stand in grp.iterrows():
        cn=stand["STAND_CN"]; r=by.get(cn)
        if r is None or r.empty: continue
        # require it to be a fully-projecting stand: inv_year>=2000
        iy=int(float(stand.get("INV_YEAR") or 2010))
        if iy<2000: continue
        picked=(state,cn,stand,r,iy); break
    if picked: break
state,cn,stand,rows,iy=picked
print("PICKED",state,cn,"iy",iy,"ntrees",len(rows))
sid=f"S{cn}"
pd_={"INVYR":iy,"LAT":stand.get("LATITUDE"),"LON":stand.get("LONGITUDE"),"ELEV":stand.get("ELEVFT") or 500,
     "SLOPE":stand.get("SLOPE") or 10,"ASPECT":stand.get("ASPECT") or 180,"STDAGE":stand.get("AGE") or 50}
gdir=G.build_gompit_config_dir(P.CONFIG_DIR,"ne","/tmp/gompit_config")
P.CONFIG_DIR=gdir
sdf=P.build_fvs_standinit(pd_,sid,"ne"); cur=RC.treeinit_for_stand(rows,sid)
print("init treeinit rows",len(cur))
for c in range(3):
    fr=P.run_fvs_projection(sdf,cur,sid,"ne",config_version="calibrated",num_cycles=1,cycle_length=5)
    yrs=sorted(fr["treelists"].keys())
    print(f"cyc{c}: yrs={yrs} exit={fr['exit_code']}")
    if not yrs: break
    g=fr["treelists"][yrs[-1]]
    print("   grown rows",len(g),"sumTPA",round(g['TPA'].sum(),3),"meanDBH",round(g['DBH'].mean(),2))
    nxt=G.treelist_to_treeinit(g,sid)
    print("   rebuilt rows",len(nxt), "cols",list(nxt.columns)[:4])
    cur=nxt
    if cur.empty: print("EMPTY break"); break
