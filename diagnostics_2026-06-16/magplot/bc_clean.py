#!/usr/bin/env python3
"""Protocol-consistent BC pairs: keep sites where t1 and t2 share the same plot_type set
and the same minimum-DBH tag limit, with plausible growth. Cancels protocol-change artifacts."""
import pandas as pd, numpy as np, math
MP="/users/PUOM0008/crsfaaron/magplot"
tr=pd.read_csv(f"{MP}/magp_trees_bc.csv",low_memory=False); tr["magp_site_id"]=tr["magp_site_id"].astype(str)
hdr=pd.read_csv(f"{MP}/magp_tree_header.csv",low_memory=False); hdr["magp_site_id"]=hdr["magp_site_id"].astype(str)
yr=hdr.groupby(["magp_site_id","meas_num"],as_index=False)["meas_year"].first()
tr=tr.merge(yr,on=["magp_site_id","meas_num"],how="left")
live=tr[(tr["tree_status"].astype(str).str.startswith("L"))&(tr["dbh"]>0)&(tr["stem_ha"]>0)&(tr["meas_year"].notna())].copy()
live["ba_tree"]=math.pi/4.0*(live["dbh"]/100.0)**2*live["stem_ha"]
live["d2w"]=live["dbh"]**2*live["stem_ha"]
# per (site,meas): BA, TPH, plot_type set, dbh tag (min dbh), QMD
agg=live.groupby(["magp_site_id","meas_num"]).agg(
    BA=("ba_tree","sum"),TPH=("stem_ha","sum"),sumd2w=("d2w","sum"),sumw=("stem_ha","sum"),
    yr=("meas_year","first"),mindbh=("dbh","min")).reset_index()
agg["QMD"]=np.sqrt(agg["sumd2w"]/agg["sumw"])
ptset=live.groupby(["magp_site_id","meas_num"])["plot_type_id"].apply(lambda s:frozenset(s.unique())).reset_index(name="pts")
agg=agg.merge(ptset,on=["magp_site_id","meas_num"])
rows=[]
for site,grp in agg.groupby("magp_site_id"):
    grp=grp.sort_values("yr")
    if grp["meas_num"].nunique()<2: continue
    t1=grp.iloc[0]; t2=grp.iloc[-1]; iv=t2["yr"]-t1["yr"]
    if not iv>=8: continue                                  # need a real interval
    if t1["pts"]!=t2["pts"]: continue                       # same subplots measured
    if abs(t1["mindbh"]-t2["mindbh"])>1.0: continue         # same tag limit
    ratio=t2["BA"]/t1["BA"] if t1["BA"]>0 else 99
    if not (0.7<=ratio<=1.8): continue                      # plausible BA change
    rows.append(dict(SITE=site,meas1=int(t1["meas_num"]),interval_years=int(iv),y1=int(t1["yr"]),
        BA_t1_m2ha=round(t1["BA"],3),TPH_t1=round(t1["TPH"],1),
        BA_t2_obs=round(t2["BA"],3),TPH_t2_obs=round(t2["TPH"],1),QMD_t2_obs=round(t2["QMD"],3),
        QMD_t1=round(t1["QMD"],3)))
out=pd.DataFrame(rows); out.to_csv(f"{MP}/magplot_ak_bc_pairs_clean.csv",index=False)
print(f"clean BC pairs: {len(out)} (from 6845)")
print(f"  interval {out['interval_years'].mean():.1f}y; BA {out['BA_t1_m2ha'].mean():.1f}->{out['BA_t2_obs'].mean():.1f} m2/ha")
out["dBA_yr"]=(out["BA_t2_obs"]-out["BA_t1_m2ha"])/out["interval_years"]
print(f"  annual BA incr median {out['dBA_yr'].median():.2f} m2/ha/yr; neg-frac {(out['dBA_yr']<0).mean():.2f}")
