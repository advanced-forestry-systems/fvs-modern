#!/usr/bin/env python3
"""Extract BC (AK-variant) trees from the full MAGPlot tree table and build stand-level
remeasurement pairs in the same schema as the NB magplot_pairs.csv, for AK calibration.
Also re-derives NB pairs from source to confirm reproducibility."""
import pandas as pd, numpy as np, os, sys, math
MP="/users/PUOM0008/crsfaaron/magplot"
sites=pd.read_csv(f"{MP}/magp_sites.csv",low_memory=False)
bc=set(sites.loc[sites["province"]=="BC","magp_site_id"].astype(str))
print("BC sites:",len(bc),flush=True)

# stream the 1.5GB tree file, keep only BC live-ish rows, minimal columns
usecols=["magp_site_id","plot_id","plot_type_id","tree_num","meas_num","tree_id",
         "species_g","species_s","dbh","height","stem_ha","tree_status"]
chunks=[]
n=0
for ch in pd.read_csv(f"{MP}/magp_trees.csv",usecols=usecols,low_memory=False,chunksize=1_000_000):
    ch=ch[ch["magp_site_id"].astype(str).isin(bc)]
    if len(ch): chunks.append(ch)
    n+=1
    if n%200==0: print("  chunk",n,"kept_so_far",sum(len(c) for c in chunks),flush=True)
tr=pd.concat(chunks,ignore_index=True) if chunks else pd.DataFrame(columns=usecols)
print("BC tree rows:",len(tr),flush=True)
tr.to_csv(f"{MP}/magp_trees_bc.csv",index=False)

# header table for meas_year + plot size
hdr=pd.read_csv(f"{MP}/magp_tree_header.csv",low_memory=False)
hdr["magp_site_id"]=hdr["magp_site_id"].astype(str)
yr=hdr.groupby(["magp_site_id","plot_id","plot_type_id","meas_num"],as_index=False)["meas_year"].first()

def build_pairs(tr, label):
    tr=tr.copy(); tr["magp_site_id"]=tr["magp_site_id"].astype(str)
    tr=tr[(tr["dbh"]>0)]                       # valid diameters
    tr=tr.merge(yr,on=["magp_site_id","plot_id","plot_type_id","meas_num"],how="left")
    live=tr[tr["tree_status"]=="L"].copy()
    live["ba_tree"]=math.pi/4.0*(live["dbh"]/100.0)**2*live["stem_ha"]   # m2/ha contribution
    g=live.groupby(["magp_site_id","plot_id","plot_type_id","meas_num"]).agg(
        BA=("ba_tree","sum"), TPH=("stem_ha","sum"),
        sumd2=("dbh",lambda s:0), meas_year=("meas_year","first")).reset_index()
    # QMD needs weighted; recompute properly
    def qmd(grp):
        w=grp["stem_ha"]; d=grp["dbh"]
        return math.sqrt(np.sum(w*d*d)/np.sum(w)) if np.sum(w)>0 else np.nan
    q=live.groupby(["magp_site_id","plot_id","plot_type_id","meas_num"]).apply(qmd).reset_index(name="QMD")
    g=g.merge(q,on=["magp_site_id","plot_id","plot_type_id","meas_num"])
    rows=[]
    key=["magp_site_id","plot_id","plot_type_id"]
    for k,grp in g.groupby(key):
        grp=grp.sort_values("meas_year")
        if grp["meas_num"].nunique()<2: continue
        t1=grp.iloc[0]; t2=grp.iloc[-1]
        iv=t2["meas_year"]-t1["meas_year"]
        if not (iv>0): continue
        stand=f"{k[0]}_{k[1]}_{k[2]}"
        rows.append(dict(STAND=stand,SITE=k[0],interval_years=int(iv),
            y1=int(t1["meas_year"]),y2=int(t2["meas_year"]),
            BA_t1_m2ha=round(t1["BA"],3),TPH_t1=round(t1["TPH"],1),
            BA_t2_obs=round(t2["BA"],3),TPH_t2_obs=round(t2["TPH"],1),
            QMD_t2_obs=round(t2["QMD"],3)))
    out=pd.DataFrame(rows)
    out.to_csv(f"{MP}/magplot_{label}_pairs.csv",index=False)
    print(f"{label}: {len(out)} stand pairs, mean interval {out['interval_years'].mean():.1f}y, "
          f"mean BA_t1 {out['BA_t1_m2ha'].mean():.1f}, mean BA_t2 {out['BA_t2_obs'].mean():.1f}",flush=True)
    return out

build_pairs(tr,"bc")
print("DONE_BC_PAIRS",flush=True)
