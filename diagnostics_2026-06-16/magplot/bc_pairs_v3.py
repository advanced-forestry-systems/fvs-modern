#!/usr/bin/env python3
"""Correct nested-plot compilation: combine ALL subplots per (site, meas_num).
stem_ha is the per-ha expansion for each tree, so the stand = sum over all live trees."""
import pandas as pd, numpy as np, math
MP="/users/PUOM0008/crsfaaron/magplot"
tr=pd.read_csv(f"{MP}/magp_trees_bc.csv",low_memory=False); tr["magp_site_id"]=tr["magp_site_id"].astype(str)
hdr=pd.read_csv(f"{MP}/magp_tree_header.csv",low_memory=False); hdr["magp_site_id"]=hdr["magp_site_id"].astype(str)
yr=hdr.groupby(["magp_site_id","meas_num"],as_index=False)["meas_year"].first()
tr=tr.merge(yr,on=["magp_site_id","meas_num"],how="left")
live=tr[(tr["tree_status"].astype(str).str.startswith("L"))&(tr["dbh"]>0)&(tr["stem_ha"]>0)&(tr["meas_year"].notna())].copy()
live["ba_tree"]=math.pi/4.0*(live["dbh"]/100.0)**2*live["stem_ha"]
live["d2w"]=live["dbh"]**2*live["stem_ha"]
g=live.groupby(["magp_site_id","meas_num"]).agg(BA=("ba_tree","sum"),TPH=("stem_ha","sum"),
    sumd2w=("d2w","sum"),sumw=("stem_ha","sum"),meas_year=("meas_year","first")).reset_index()
g["QMD"]=np.sqrt(g["sumd2w"]/g["sumw"])
rows=[]
for site,grp in g.groupby("magp_site_id"):
    grp=grp.sort_values("meas_year")
    if grp["meas_num"].nunique()<2: continue
    t1=grp.iloc[0]; t2=grp.iloc[-1]; iv=t2["meas_year"]-t1["meas_year"]
    if not iv>0: continue
    rows.append(dict(STAND=site,SITE=site,interval_years=int(iv),y1=int(t1["meas_year"]),y2=int(t2["meas_year"]),
        meas1=int(t1["meas_num"]),
        BA_t1_m2ha=round(t1["BA"],3),TPH_t1=round(t1["TPH"],1),
        BA_t2_obs=round(t2["BA"],3),TPH_t2_obs=round(t2["TPH"],1),QMD_t2_obs=round(t2["QMD"],3),
        BA_t1_for_chk=round(t1["BA"],2)))
out=pd.DataFrame(rows); out.to_csv(f"{MP}/magplot_ak_bc_pairs_v3.csv",index=False)
print(f"BC site-level pairs: {len(out)}; interval {out['interval_years'].mean():.1f}y; "
      f"BA {out['BA_t1_m2ha'].mean():.1f}->{out['BA_t2_obs'].mean():.1f} m2/ha; "
      f"QMD_t2 {out['QMD_t2_obs'].mean():.1f}cm; TPH_t1 {out['TPH_t1'].mean():.0f}",flush=True)
# growth sanity: distribution of annual BA increment
out["dBA_yr"]=(out["BA_t2_obs"]-out["BA_t1_m2ha"])/out["interval_years"]
print("annual BA increment m2/ha/yr: median %.2f  p10 %.2f p90 %.2f  neg-frac %.2f"%(
    out["dBA_yr"].median(),out["dBA_yr"].quantile(.1),out["dBA_yr"].quantile(.9),(out["dBA_yr"]<0).mean()),flush=True)
print("DONE",flush=True)
