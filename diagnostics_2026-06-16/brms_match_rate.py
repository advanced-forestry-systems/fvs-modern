#!/usr/bin/env python3
"""Per-variant brms SDImax match rate and distribution (red-team master-table item, 2026-06-18).

Match rate = fraction of FIA plots in a variant's states that have a brms site-specific max SDI estimate
(~/overthin_work/brms_SDImax.csv), plus the SDImax distribution per variant (English, metric/2.471). This
quantifies brms coverage per variant and surfaces the Lake States SDImax anomaly the red-team flagged.
Output: brms_match_rate.csv. Env OUT.
"""
import os, csv
import pandas as pd, numpy as np
from pathlib import Path
FIA="/fs/scratch/PUOM0008/crsfaaron/FIA"; SDIc=2.471
OUT=os.environ.get("OUT",os.path.expanduser("~/overthin_work/brms_match_rate.csv"))
FIPS={"AL":1,"CA":6,"CO":8,"CT":9,"GA":13,"ID":16,"IL":17,"IN":18,"ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NY":36,"ND":38,"OR":41,"SC":45,"SD":46,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
VARMAP={"ne":["CT","ME","MA","NH","NY","VT"],"acd":["ME","NH","VT"],"sn":["AL","GA","SC","MS"],"ls":["MI","MN","WI"],
        "cs":["IL","IN","MO"],"ie":["ID","MT"],"kt":["MT","ID"],"ci":["ID"],"cr":["CO","WY"],"ut":["UT","NV"],
        "ca":["CA"],"nc":["CA","OR"],"ec":["OR","WA"],"wc":["OR","WA"],"pn":["OR","WA"]}
brms=pd.read_csv(os.path.expanduser("~/overthin_work/brms_SDImax.csv"))
brms["sdi_eng"]=brms["SDImax.median"]/SDIc
# brms plot keys + per-state stats
brms_by_state={}
for sc,g in brms.groupby("STATECD"):
    brms_by_state[int(sc)]=g
brmskey=set(zip(brms.STATECD.astype(int),brms.UNITCD.astype(int),brms.COUNTYCD.astype(int),brms.PLOT.astype(int)))
rows=[]
for var,states in VARMAP.items():
    nfia=0; nmatch=0; sdivals=[]
    for s in states:
        pf=Path(FIA)/f"{s}_PLOT.csv"
        if not pf.exists(): continue
        pp=pd.read_csv(pf,usecols=lambda c:c in ("STATECD","UNITCD","COUNTYCD","PLOT"),low_memory=False).dropna()
        nfia+=len(pp)
        for _,r in pp.iterrows():
            if (int(r.STATECD),int(r.UNITCD),int(r.COUNTYCD),int(r.PLOT)) in brmskey: nmatch+=1
        sc=FIPS[s]
        if sc in brms_by_state: sdivals+=list(brms_by_state[sc].sdi_eng.values)
    if nfia==0: continue
    sd=np.array(sdivals,float)
    rows.append({"variant":var,"n_fia_plots":nfia,"n_brms_match":nmatch,"match_rate_pct":round(100*nmatch/nfia,1),
                 "SDImax_median":round(float(np.median(sd)),0) if len(sd) else float("nan"),
                 "SDImax_p10":round(float(np.percentile(sd,10)),0) if len(sd) else float("nan"),
                 "SDImax_p90":round(float(np.percentile(sd,90)),0) if len(sd) else float("nan")})
    print("%-4s FIA=%-7d brms_match=%-7d (%.1f%%)  SDImax med=%.0f [%.0f-%.0f]"%(
        var,nfia,nmatch,rows[-1]["match_rate_pct"],rows[-1]["SDImax_median"],rows[-1]["SDImax_p10"],rows[-1]["SDImax_p90"]))
with open(OUT,"w",newline="") as fh: w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote",OUT); print("DONE_MATCHRATE")
