import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-modern")
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-conus/python")
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import perseus_100yr_projection as P
import run_conus_task_wo1_v2 as V

# pick AK b0, first matched stand_cn from existing out (698882992126144)
variant="AK"; state="AK"; batch_id=0; batch_size=5000
si_all = pd.read_csv("/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant/standinit_AK.csv", low_memory=False)
si = si_all.iloc[0:batch_size].reset_index(drop=True)
si["STAND_CN"] = si["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
tt = pd.read_csv("/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h/AK_FVS_TREEINIT_PLOT.csv", low_memory=False)
tt["STAND_CN"] = tt["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
by_cn = {k:v for k,v in tt.groupby("STAND_CN")}

target="698882992126144"
stand = si[si["STAND_CN"]==target].iloc[0]
cn = stand["STAND_CN"]; sid=f"S{cn}"
inv_year = int(float(stand.get("INV_YEAR") or 2010))
plot_data = {"INVYR": inv_year, "LAT": stand.get("LATITUDE"), "LON": stand.get("LONGITUDE"),
             "ELEV": stand.get("ELEVFT") or 500, "SLOPE": stand.get("SLOPE") or 10,
             "ASPECT": stand.get("ASPECT") or 180, "STDAGE": stand.get("AGE") or 50}
sdf = P.build_fvs_standinit(plot_data, sid, variant.lower())
tdf = V.treeinit_for_stand(by_cn[cn], sid)
nsbe = P.NSBECalculator(P.NSBE_ROOT)
fr = P.run_fvs_projection(sdf, tdf, sid, variant.lower(), config_version=None, num_cycles=20, cycle_length=5)
print("STAND_CN", cn, "inv_year", inv_year)
tl0 = None
for cy, tl in sorted(fr["treelists"].items()):
    py = cy - inv_year
    if py < 0: continue
    agb = P.compute_plot_agb(tl, nsbe)
    sm = V.stand_metrics(tl)
    if tl0 is None: tl0 = (cy, tl, agb, sm)
    print(f"YEAR={cy} PY={py} cols={list(tl.columns)[:6]} AGB={agb:.4f} BA={sm['BA_FT2AC']} QMD={sm['QMD_IN']} TPH={sm['TPH']}")
    if py>=10: break

# hand calc on year-0 treelist using pandas vectorized, independent of helper
cy, tl, agb, sm = tl0
dbh = pd.to_numeric(tl.get("DBH", tl.get("Dbh")), errors="coerce")
tpa = pd.to_numeric(tl.get("TPA", tl.get("Tpa")), errors="coerce")
m = (~dbh.isna()) & (dbh>=1.0) & (~tpa.isna())
d=dbh[m].astype(float); t=tpa[m].astype(float)
ba_hand = (np.pi/(4*144))*(d*d*t).sum()
qmd_hand = np.sqrt((d*d*t).sum()/t.sum())
tph_hand = t.sum()*2.4710538
print("HANDCALC y0:", f"BA={ba_hand:.4f} QMD={qmd_hand:.4f} TPH={tph_hand:.4f}")
print("MATCH:", abs(ba_hand-sm['BA_FT2AC'])<1e-2, abs(qmd_hand-sm['QMD_IN'])<1e-2, abs(tph_hand-sm['TPH'])<1e-2)
# cross-check AGB vs existing out_conus_wo1 row
prev = pd.read_csv("/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_wo1/conus_ak_b0.csv")
prow = prev[(prev["STAND_CN"]==int(cn)) & (prev["YEAR"]==cy) & (prev["CONFIG"]=="default")]
print("PREV AGB(default,y0):", prow["AGB_TONS_AC"].values, "vs NEW:", round(agb,4))
print("FULL treelist cols:", list(tl.columns))
