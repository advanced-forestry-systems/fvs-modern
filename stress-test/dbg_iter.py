import os,sys
import numpy as np, pandas as pd
PR=os.environ.get("FVS_PROJECT_ROOT",os.path.expanduser("~/fvs-modern"))
sys.path.insert(0,PR); sys.path.insert(0,os.path.join(PR,"calibration","python"))
sys.path.insert(0,"/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import perseus_100yr_projection as P
import run_gompit_projection as G
import run_conus_task_fvstreeinit as RC
from greg_mortality import GregMortality
greg=GregMortality("/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv")
gdir=G.build_gompit_config_dir(P.CONFIG_DIR,"ne","/tmp/gompit_config")
# load one good stand
si=pd.read_csv("/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant/standinit_NE.csv",low_memory=False)
si["STAND_CN"]=si["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
tt=pd.read_csv("/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit/WI_FVS_TREEINIT_PLOT.csv",low_memory=False)
tt["STAND_CN"]=tt["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
cn="55919559010538"; rows=tt[tt.STAND_CN==cn]
stand=si[si.STAND_CN==cn].iloc[0]; iy=int(float(stand.get("INV_YEAR") or 2010))
pd_={"INVYR":iy,"LAT":stand.get("LATITUDE"),"LON":stand.get("LONGITUDE"),"ELEV":stand.get("ELEVFT") or 500,
     "SLOPE":stand.get("SLOPE") or 10,"ASPECT":stand.get("ASPECT") or 180,"STDAGE":stand.get("AGE") or 50}
sid=f"S{cn}"
P.CONFIG_DIR=gdir
sdf=P.build_fvs_standinit(pd_,sid,"ne"); cur=RC.treeinit_for_stand(rows,sid)
print("init treeinit rows",len(cur),"cols",list(cur.columns))
for c in range(3):
    fr=P.run_fvs_projection(sdf,cur,sid,"ne",config_version="calibrated",num_cycles=1,cycle_length=5)
    yrs=sorted(fr["treelists"].keys())
    print(f"cycle {c}: yrs={yrs} exit={fr['exit_code']} ntl0={len(fr['treelists'][yrs[0]]) if yrs else 0}")
    if not yrs: break
    grown=fr["treelists"][yrs[-1]]
    print("   grown rows",len(grown),"sumTPA",round(grown['TPA'].sum(),2))
    nxt=G.treelist_to_treeinit(grown,sid)
    print("   rebuilt treeinit rows",len(nxt))
    cur=nxt
    if cur.empty: print("   EMPTY -> would break"); break
