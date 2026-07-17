#!/usr/bin/env python3
import os,sys,numpy as np,pandas as pd
PR=os.environ.get("FVS_PROJECT_ROOT","/users/PUOM0008/crsfaaron/fvs-modern")
sys.path.insert(0,PR); sys.path.insert(0,"/users/PUOM0008/crsfaaron/fvs-conus/python")
import perseus_100yr_projection as P
SCR="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
tt=pd.read_csv(os.path.join(SCR,"..","FIA_fresh","treeinit_h","CA_FVS_TREEINIT_PLOT.csv"),low_memory=False)
tt["STAND_CN"]=tt["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
g=tt.groupby("STAND_CN")
cand=None
for cn,grp in g:
    grp2=grp[(grp["DIAMETER"]>=1.0)&(grp["TREE_COUNT"]>0)]
    if 4<=len(grp2)<=8 and grp2["DIAMETER"].max()<30:
        cand=(cn,grp2); break
cn,fvs_rows=cand
BA_CONST=np.pi/(4*144)
obs_ba=(BA_CONST*fvs_rows["DIAMETER"]**2*fvs_rows["TREE_COUNT"]).sum()
obs_tpa=fvs_rows["TREE_COUNT"].sum()
sid=f"S{cn}"
plot_data={"INVYR":2014,"LAT":40,"LON":-120,"ELEV":1000,"SLOPE":10,"ASPECT":180,"STDAGE":50}
sdf=P.build_fvs_standinit(plot_data,sid,"cr")
# --- THE FIX: per-acre design, no large-tree breakpoint ---
sdf["inv_plot_size"]=1.0
sdf["num_plots"]=1
sdf["brk_dbh"]=999.0
recs=[]
for i,t in enumerate(fvs_rows.itertuples(index=False)):
    recs.append({"stand_id":sid,"plot_id":1,"tree_id":i+1,
                 "tree_count":float(getattr(t,"TREE_COUNT",1.0)),
                 "species":int(float(getattr(t,"SPECIES",0))),
                 "diameter":round(float(getattr(t,"DIAMETER")),1),
                 "ht":round(float(getattr(t,"HT",0) or 0),0),
                 "crratio":int(float(getattr(t,"CRRATIO",0) or 0))})
tdf=pd.DataFrame(recs)
fr=P.run_fvs_projection(sdf,tdf,sid,"cr",config_version=None,num_cycles=1,cycle_length=5)
tls=fr["treelists"]; y0=min(tls.keys())
tl=tls[y0]; tpa_col="TPA" if "TPA" in tl.columns else "Tpa"; dbh_col="DBH" if "DBH" in tl.columns else "Dbh"
eng_tpa=tl[tpa_col].astype(float).sum()
eng_ba=(BA_CONST*tl[dbh_col].astype(float)**2*tl[tpa_col].astype(float)).sum()
print(f"FIXED standinit: inv_plot_size=1 num_plots=1 brk_dbh=999")
print("FVS year0 treelist (DBH, TPA):")
print(tl[[dbh_col,tpa_col]].to_string(index=False))
print(f"\nOBSERVED/EQ year0 BA = {obs_ba:.3f}  TPA={obs_tpa:.3f}")
print(f"ENGINE(fixed) year0 BA = {eng_ba:.3f}  TPA={eng_tpa:.3f}")
print(f"RATIO engine/obs BA = {eng_ba/obs_ba:.4f}  TPA = {eng_tpa/obs_tpa:.4f}")
