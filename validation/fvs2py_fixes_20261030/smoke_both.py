"""Same 5 NE + 5 ACD stands as smoke_summary.py, run through run_scenario() of BOTH
calibration/python/bakuzis_uncertainty_comparison.py and bakuzis_100yr_comparison.py
(one never-closed FVS object per call). usage: smoke_both.py <out.csv>"""
import sys, logging, warnings
import pandas as pd
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
import bakuzis_uncertainty_comparison as buc
import bakuzis_100yr_comparison as b100
logging.disable(logging.CRITICAL)
from fia_stand_generator import generate_real_stand
rows = []
for v in ("ne", "acd"):
    scen = buc.enumerate_scenarios(buc.get_species_groups(v))
    sp, site, dens, num = scen[13]
    stands = generate_real_stand(variant=v, site_class=site, density_class=dens, n_plots=5, seed=42 + num)
    for rep, (sdf, tdf, _) in enumerate(stands, start=1):
        sid = f"SM{num:02d}R{rep}"; sdf = sdf.copy(); tdf = tdf.copy(); sdf["stand_id"] = sid; tdf["stand_id"] = sid
        for cfg in ("default", "calibrated"):
            for name, mod in (("uncertainty", buc), ("100yr", b100)):
                s = mod.run_scenario(sdf, tdf, sid, v, cfg).copy()
                s["variant"] = v.upper(); s["replicate"] = rep; s["config"] = cfg; s["script"] = name
                rows.append(s)
d = pd.concat(rows, ignore_index=True)
d.to_csv(sys.argv[1], index=False)
print("SMOKE_BOTH_DONE", len(d), "rows", d.shape[1], "columns")
