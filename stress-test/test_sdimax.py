#!/usr/bin/env python3
# test_sdimax.py - run a state sample to ~2100 (no disturbance, calibrated config)
# and report mean aboveground live carbon density. Used to compare the current
# calibrated config (NA max-SDI species revert to high built-in defaults) against
# the SDImax-fix config (NA filled with the calibrated median). Run twice with the
# config/calibrated/<variant>.json swapped between current and the sdifix version.
import os, sys, argparse
import numpy as np, pandas as pd
PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import perseus_100yr_projection as P
import run_conus_task_fvstreeinit as RC
FIPS_INV = {v: k for k, v in RC.FIPS.items()}
TONS_AC_TO_MGHA = 2.241702; C_FRACTION = 0.47

ap = argparse.ArgumentParser()
ap.add_argument("--variant", required=True); ap.add_argument("--state", required=True)
ap.add_argument("--standinit-dir", required=True); ap.add_argument("--treeinit-dir", required=True)
ap.add_argument("--n-plots", type=int, default=40); ap.add_argument("--num-cycles", type=int, default=16)
ap.add_argument("--tag", default="run"); ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
variant, ST = a.variant.lower(), a.state.upper(); fips = FIPS_INV[ST]
nsbe = P.NSBECalculator(P.NSBE_ROOT)

si = pd.read_csv(os.path.join(a.standinit_dir, f"standinit_{variant.upper()}.csv"), low_memory=False)
si["STAND_CN"] = si["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
si = si[si["STATE"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "") == str(fips)]
tt = pd.read_csv(os.path.join(a.treeinit_dir, f"{ST}_FVS_TREEINIT_PLOT.csv"), low_memory=False)
tt["STAND_CN"] = tt["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
by_cn = {k: v for k, v in tt.groupby("STAND_CN")}
si = si[si["STAND_CN"].isin(by_cn)]
rng = np.random.default_rng(a.seed)
si = si.iloc[rng.choice(len(si), size=min(a.n_plots, len(si)), replace=False)]

per_year = {}
n_ok = 0
for _, stand in si.iterrows():
    cn = stand["STAND_CN"]; iy = int(float(stand.get("INV_YEAR") or 2010))
    pdat = {"INVYR": iy, "LAT": stand.get("LATITUDE"), "LON": stand.get("LONGITUDE"),
            "ELEV": stand.get("ELEVFT") or 500, "SLOPE": stand.get("SLOPE") or 10,
            "ASPECT": stand.get("ASPECT") or 180, "STDAGE": stand.get("AGE") or 50}
    sid = f"S{cn}"
    try:
        sdf = P.build_fvs_standinit(pdat, sid, variant); tdf = RC.treeinit_for_stand(by_cn[cn], sid)
        if tdf.empty: continue
        fr = P.run_fvs_projection(sdf, tdf, sid, variant, config_version="calibrated",
                                  num_cycles=a.num_cycles, cycle_length=5)
        for cy, tl in fr["treelists"].items():
            yr = 2025 + (cy - iy)
            agc = P.compute_plot_agb(tl, nsbe) * TONS_AC_TO_MGHA * C_FRACTION
            per_year.setdefault(yr, []).append(float(agc))
        n_ok += 1
    except Exception:
        continue
yrs = sorted(per_year)
def m(y): return float(np.mean(per_year[y])) if y in per_year else float("nan")
y0 = yrs[0]; y100 = max(y for y in yrs if y <= 2100)
print(f"[{a.tag}] n_plots={n_ok}  {y0}={m(y0):.1f}  {y100}={m(y100):.1f}  growth={100*(m(y100)-m(y0))/m(y0):.0f}%  (MgC/ha)")
