#!/usr/bin/env python3
"""fvs_posterior_uncertainty.py -- propagate Bayesian posterior draws through FVS
to get a parametric carbon CI (the parameter-uncertainty layer).

For a calibrated variant, loads its posterior draws JSON
(config/calibrated/<variant>_draws.json, 500 draws x 6 components), draws an
N-sample, and for each draw runs a per-state plot SUBSAMPLE through the standard
projection -- injecting the draw's parameters via UncertaintyEngine.
generate_keywords_for_draw, with no change to production code (we monkeypatch
FvsConfigLoader.generate_keywords to return the draw keywords, so the existing
config_version="calibrated" keyfile path carries the draw).

The ensemble across draws gives the parameter-uncertainty distribution of the
state mean carbon density at each year; percentiles -> CI. This is distinct from
(and complementary to) the plot-sampling CI and the engine-spread (structural)
band.

Output: <out>/posterior_<ST>.csv  (ST,variant,year,mean,p2_5,p50,p97_5,n_draws,
        n_plots)  in Mg C/ha density.

Usage:
  python3 fvs_posterior_uncertainty.py --variant ne --state ME \
     --standinit-dir .../standinit_by_variant --treeinit-dir .../treeinit_h \
     --n-draws 30 --n-plots 80 --out posterior_unc [--seed 1]
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import perseus_100yr_projection as P            # noqa: E402
import config.config_loader as CL               # noqa: E402
from config.uncertainty import UncertaintyEngine  # noqa: E402
import run_conus_task_fvstreeinit as RC          # noqa: E402

FIPS_INV = {v: k for k, v in RC.FIPS.items()}
TONS_AC_TO_MGHA = 2.241702
C_FRACTION = 0.47


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True)
    ap.add_argument("--state", required=True)
    ap.add_argument("--standinit-dir", required=True)
    ap.add_argument("--treeinit-dir", required=True)
    ap.add_argument("--n-draws", type=int, default=30)
    ap.add_argument("--n-plots", type=int, default=80)
    ap.add_argument("--num-cycles", type=int, default=20)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--draw-idx", type=int, default=-1,
                    help="single-draw mode: run only this posterior draw index "
                         "(for SLURM array parallelism); writes posterior_<ST>_d<i>.csv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    variant, ST = a.variant.lower(), a.state.upper()
    fips = FIPS_INV[ST]

    # default config + uncertainty engine
    dcfg_path = os.path.join(PROJECT_ROOT, "config", f"{variant}.json")
    default_config = json.load(open(dcfg_path)) if os.path.exists(dcfg_path) else {}
    eng = UncertaintyEngine(
        variant=variant, config_dir=os.path.join(PROJECT_ROOT, "config"),
        draws_path=os.path.join(PROJECT_ROOT, "config", "calibrated",
                                f"{variant}_draws.json"), seed=a.seed)
    print(f"{variant}/{ST}: {eng.n_draws} draws available; "
          f"sampling {a.n_draws}", flush=True)

    # subsample plots for this state/variant that have trees
    si = pd.read_csv(os.path.join(a.standinit_dir, f"standinit_{variant.upper()}.csv"),
                     low_memory=False)
    si["STAND_CN"] = si["STAND_CN"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else "")
    si = si[si["STATE"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
            == str(fips)]
    tfile = os.path.join(a.treeinit_dir, f"{ST}_FVS_TREEINIT_PLOT.csv")
    tt = pd.read_csv(tfile, low_memory=False)
    tt["STAND_CN"] = tt["STAND_CN"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else "")
    by_cn = {k: v for k, v in tt.groupby("STAND_CN")}
    si = si[si["STAND_CN"].isin(by_cn)]
    rng = np.random.default_rng(a.seed)
    si = si.iloc[rng.choice(len(si), size=min(a.n_plots, len(si)), replace=False)]
    print(f"  {len(si)} subsample plots", flush=True)

    nsbe = P.NSBECalculator(P.NSBE_ROOT)
    orig_gk = CL.FvsConfigLoader.generate_keywords

    # build sdf/tdf once per plot
    plots = []
    for _, stand in si.iterrows():
        cn = stand["STAND_CN"]
        iy = int(float(stand.get("INV_YEAR") or 2010))
        pdat = {"INVYR": iy, "LAT": stand.get("LATITUDE"),
                "LON": stand.get("LONGITUDE"), "ELEV": stand.get("ELEVFT") or 500,
                "SLOPE": stand.get("SLOPE") or 10, "ASPECT": stand.get("ASPECT") or 180,
                "STDAGE": stand.get("AGE") or 50}
        sid = f"S{cn}"
        try:
            sdf = P.build_fvs_standinit(pdat, sid, variant)
            tdf = RC.treeinit_for_stand(by_cn[cn], sid)
        except Exception:
            continue
        if not tdf.empty:
            plots.append((sid, iy, sdf, tdf))
    print(f"  {len(plots)} plots built", flush=True)

    if a.draw_idx >= 0:                          # single-draw (array) mode
        didx = a.draw_idx % eng.n_draws
        kw = eng.generate_keywords_for_draw(eng.get_draw(didx), default_config, didx)
        CL.FvsConfigLoader.generate_keywords = (
            lambda self, include_comments=False, _kw=kw: _kw)
        per_year = {}
        for sid, iy, sdf, tdf in plots:
            try:
                fr = P.run_fvs_projection(sdf, tdf, sid, variant,
                                          config_version="calibrated",
                                          num_cycles=a.num_cycles, cycle_length=5)
                for cy, tl in fr["treelists"].items():
                    yr = 2025 + (cy - iy)
                    agc = P.compute_plot_agb(tl, nsbe) * TONS_AC_TO_MGHA * C_FRACTION
                    per_year.setdefault(yr, []).append(float(agc))
            except Exception:
                continue
        rows = [{"ST": ST, "variant": variant, "draw": didx, "year": yr,
                 "mean_density": round(float(np.mean(v)), 4), "n_plots": len(v)}
                for yr, v in sorted(per_year.items()) if v]
        pd.DataFrame(rows).to_csv(
            os.path.join(a.out, f"posterior_{ST}_d{didx}.csv"), index=False)
        print(f"draw {didx}: wrote {len(rows)} years, {len(plots)} plots",
              flush=True)
        return

    draw_idxs = [eng.sample_draw_index() for _ in range(a.n_draws)]
    # ens[year] -> list over draws of mean density
    ens = {}
    for di, didx in enumerate(draw_idxs):
        kw = eng.generate_keywords_for_draw(eng.get_draw(didx), default_config, didx)
        CL.FvsConfigLoader.generate_keywords = (
            lambda self, include_comments=False, _kw=kw: _kw)
        per_year = {}
        for sid, iy, sdf, tdf in plots:
            try:
                fr = P.run_fvs_projection(sdf, tdf, sid, variant,
                                          config_version="calibrated",
                                          num_cycles=a.num_cycles, cycle_length=5)
                for cy, tl in fr["treelists"].items():
                    yr = 2025 + (cy - iy)
                    agc = P.compute_plot_agb(tl, nsbe) * TONS_AC_TO_MGHA * C_FRACTION
                    per_year.setdefault(yr, []).append(float(agc))
            except Exception:
                continue
        for yr, vals in per_year.items():
            if vals:
                ens.setdefault(yr, []).append(float(np.mean(vals)))
        if (di + 1) % 5 == 0:
            print(f"  draw {di+1}/{a.n_draws}", flush=True)
    CL.FvsConfigLoader.generate_keywords = orig_gk

    rows = []
    for yr in sorted(ens):
        v = np.array(ens[yr])
        rows.append({"ST": ST, "variant": variant, "year": yr,
                     "mean": round(float(v.mean()), 3),
                     "p2_5": round(float(np.percentile(v, 2.5)), 3),
                     "p50": round(float(np.percentile(v, 50)), 3),
                     "p97_5": round(float(np.percentile(v, 97.5)), 3),
                     "n_draws": len(v), "n_plots": len(plots)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.out, f"posterior_{ST}.csv"), index=False)
    print(out.to_string(index=False))
    print(f"wrote {a.out}/posterior_{ST}.csv")


if __name__ == "__main__":
    main()
