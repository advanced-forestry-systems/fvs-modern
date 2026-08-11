#!/usr/bin/env python3
import pandas as pd, numpy as np, math
MP="/users/PUOM0008/crsfaaron/magplot"
def build(src,label):
    tr=pd.read_csv(src,low_memory=False); tr["magp_site_id"]=tr["magp_site_id"].astype(str)
    tr=tr[tr["dbh"]>0]
    hdr=pd.read_csv(f"{MP}/magp_tree_header.csv",low_memory=False); hdr["magp_site_id"]=hdr["magp_site_id"].astype(str)
    yr=hdr.groupby(["magp_site_id","meas_num"],as_index=False)["meas_year"].first()   # site-meas year
    tr=tr.merge(yr,on=["magp_site_id","meas_num"],how="left")
    live=tr[(tr["tree_status"].astype(str).str.startswith("L"))&(tr["meas_year"].notna())&(tr["stem_ha"]>0)].copy()
    live["ba_tree"]=math.pi/4.0*(live["dbh"]/100.0)**2*live["stem_ha"]
    live["d2w"]=live["dbh"]**2*live["stem_ha"]
    k=["magp_site_id","plot_type_id","meas_num"]
    g=live.groupby(k).agg(BA=("ba_tree","sum"),TPH=("stem_ha","sum"),
        sumd2w=("d2w","sum"),sumw=("stem_ha","sum"),meas_year=("meas_year","first")).reset_index()
    g["QMD"]=np.sqrt(g["sumd2w"]/g["sumw"])
    rows=[]
    for kk,grp in g.groupby(["magp_site_id","plot_type_id"]):
        grp=grp.sort_values("meas_year")
        if grp["meas_num"].nunique()<2: continue
        t1=grp.iloc[0]; t2=grp.iloc[-1]; iv=t2["meas_year"]-t1["meas_year"]
        if not iv>0: continue
        rows.append(dict(STAND=f"{kk[0]}_{kk[1]}",SITE=kk[0],interval_years=int(iv),
            y1=int(t1["meas_year"]),y2=int(t2["meas_year"]),
            BA_t1_m2ha=round(t1["BA"],3),TPH_t1=round(t1["TPH"],1),
            BA_t2_obs=round(t2["BA"],3),TPH_t2_obs=round(t2["TPH"],1),QMD_t2_obs=round(t2["QMD"],3)))
    out=pd.DataFrame(rows); out.to_csv(f"{MP}/magplot_{label}_pairs.csv",index=False)
    if len(out):
        print(f"{label}: {len(out)} pairs; interval {out['interval_years'].mean():.1f}y; "
              f"BA {out['BA_t1_m2ha'].mean():.1f}->{out['BA_t2_obs'].mean():.1f}; "
              f"QMD_t2 {out['QMD_t2_obs'].mean():.1f}cm; TPH_t1 {out['TPH_t1'].mean():.0f}",flush=True)
    else:
        print(f"{label}: 0 pairs",flush=True)
    return out,live
out,live=build(f"{MP}/magp_trees_bc.csv","bc")
sp=live.groupby(["species_g","species_s"]).size().sort_values(ascending=False).head(18)
print("BC top species:",flush=True)
for (gn,ep),c in sp.items(): print(f"  {gn} {ep}: {c}",flush=True)
print("DONE",flush=True)
