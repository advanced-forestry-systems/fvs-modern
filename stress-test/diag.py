#!/usr/bin/env python3
import os,sys,glob,numpy as np,pandas as pd
PR=os.environ.get("FVS_PROJECT_ROOT","/users/PUOM0008/crsfaaron/fvs-modern")
sys.path.insert(0,PR); sys.path.insert(0,"/users/PUOM0008/crsfaaron/fvs-conus/python")
import perseus_100yr_projection as P
SCR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
# pick a CR stand (CA) that has trees, with moderate BA
tt=pd.read_csv(os.path.join(SCR,"..","FIA_fresh","treeinit_h","CA_FVS_TREEINIT_PLOT.csv"),low_memory=False)
tt["STAND_CN"]=tt["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
# group sizes
g=tt.groupby("STAND_CN")
# pick a stand with ~5-10 trees and reasonable DBH
cand=None
for cn,grp in g:
    grp2=grp[(grp["DIAMETER"]>=1.0)&(grp["TREE_COUNT"]>0)]
    if 4<=len(grp2)<=8 and grp2["DIAMETER"].max()<30:
        cand=(cn,grp2); break
cn,fvs_rows=cand
BA_CONST=np.pi/(4*144)
obs_ba=(BA_CONST*fvs_rows["DIAMETER"]**2*fvs_rows["TREE_COUNT"]).sum()
print(f"STAND_CN={cn} n_trees={len(fvs_rows)}")
print("INPUT treeinit (DIAMETER, TREE_COUNT):")
print(fvs_rows[["SPECIES","DIAMETER","TREE_COUNT"]].to_string(index=False))
print(f"OBSERVED year0 BA = {obs_ba:.3f} ft2/ac (sum pi/4 D^2 TPA/144)")
print(f"  sum TREE_COUNT (=TPA) = {fvs_rows['TREE_COUNT'].sum():.3f}")
sid=f"S{cn}"
plot_data={"INVYR":2014,"LAT":fvs_rows.iloc[0].get("LAT",40),"LON":-120,"ELEV":1000,"SLOPE":10,"ASPECT":180,"STDAGE":50}
sdf=P.build_fvs_standinit(plot_data,sid,"cr")
print("\nstandinit inv_plot_size,num_plots,brk_dbh,baf:",
      float(sdf.iloc[0]["inv_plot_size"]),int(sdf.iloc[0]["num_plots"]),float(sdf.iloc[0]["brk_dbh"]),float(sdf.iloc[0]["basal_area_factor"]))
# build treeinit exactly as v2 driver does
recs=[]
for i,t in enumerate(fvs_rows.itertuples(index=False)):
    d=getattr(t,"DIAMETER",np.nan)
    recs.append({"stand_id":sid,"plot_id":1,"tree_id":i+1,
                 "tree_count":float(getattr(t,"TREE_COUNT",1.0)),
                 "species":int(float(getattr(t,"SPECIES",0))),
                 "diameter":round(float(d),1),
                 "ht":round(float(getattr(t,"HT",0) or 0),0),
                 "crratio":int(float(getattr(t,"CRRATIO",0) or 0))})
tdf=pd.DataFrame(recs)
fr=P.run_fvs_projection(sdf,tdf,sid,"cr",config_version=None,num_cycles=1,cycle_length=5)
tls=fr["treelists"]; print("\nFVS exit:",fr["exit_code"]," treelist years:",sorted(tls.keys()))
y0=min(tls.keys()) if tls else None
if y0 is not None:
    tl=tls[y0]
    cols=[c for c in tl.columns if c.lower() in ("dbh","tpa","speciesfia","species")]
    print("FVS year0 treelist:"); print(tl[cols].to_string(index=False))
    tpa_col="TPA" if "TPA" in tl.columns else "Tpa"
    dbh_col="DBH" if "DBH" in tl.columns else "Dbh"
    eng_ba=(BA_CONST*tl[dbh_col].astype(float)**2*tl[tpa_col].astype(float)).sum()
    print(f"\nENGINE year0 BA = {eng_ba:.3f} ft2/ac   sum FVS TPA={tl[tpa_col].astype(float).sum():.3f}")
    print(f"RATIO engine/obs = {eng_ba/obs_ba:.4f}   TPA ratio = {tl[tpa_col].astype(float).sum()/fvs_rows['TREE_COUNT'].sum():.4f}")
