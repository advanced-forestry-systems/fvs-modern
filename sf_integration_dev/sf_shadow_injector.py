#!/usr/bin/env python3
"""
sf_shadow_injector.py  --  Stage-1 SHADOW injection of the species-free
equations into the live FVS engine via the fvs2py stop-point-5 hook.

WHAT THIS IS
  The fitted species-free component equations currently live only in the R
  projector and as banked bundles / a conus_sf config block in
  config/config_loader.py (sf_linear_predictor). Nothing applies them inside a
  real FVS run. This script wires the first, zero-risk stage of that injection:
  it boots the real FVS engine, stops every cycle at restart code 5 (after the
  engine has computed diameter growth (dg) and mortality but BEFORE applying
  them), reads the per-tree engine increments via get_tree_attr, and LOGS them
  (no set_tree_attr -> zero behaviour change). It also confirms the conus_sf
  config block is loadable for the variant so Stage 2 (compute species-free
  dg/htg per tree and compare / inject) can wire straight on top.

  Reuses perseus_100yr_projection for keyfile/standinit construction so the
  boot path is identical to the validated engine arms (wo1_v4).

STAGES (staged rollout from 20260614_engine_wiring_design.md)
  1  shadow  : GET only, log sf-vs-engine increments per tree   <-- THIS FILE
  2  single  : SET dg only, benchmark vs observed
  3  full    : SET dg+htg+mort with localized maxSDI
  4  blended : per-species shrinkage blend, tune kappa

USAGE (smoke, one synthetic NE stand)
  python3 sf_shadow_injector.py --variant ne --smoke --out shadow_ne_smoke.csv
"""
from __future__ import annotations
import argparse, os, sys, tempfile, sqlite3, json
import numpy as np, pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.expanduser("~/fvs-conus/python"))
import perseus_100yr_projection as P   # KEYFILE_TEMPLATE, build_fvs_standinit, CONFIG_DIR, FVS_LIB_DIR


# ---- per-acre seeding fix (identical to run_conus_task_wo1_v4.py) ----------
_orig_build = P.build_fvs_standinit
def build_standinit_peracre(plot_data, stand_id, variant):
    sdf = _orig_build(plot_data, stand_id, variant)
    sdf["inv_plot_size"] = 1.0
    sdf["num_plots"] = 1
    sdf["brk_dbh"] = 999.0
    sdf["basal_area_factor"] = 0.0
    return sdf


def synthetic_ne_stand():
    """A small, plausible NE mixed stand (per-acre TPA) for the boot smoke."""
    plot_data = {"INVYR": 2010, "LAT": 45.0, "LON": -69.0, "ELEV": 500,
                 "SLOPE": 10, "ASPECT": 180, "STDAGE": 50}
    # SPCD: 12 balsam fir, 97 red spruce, 316 red maple, 318 sugar maple, 371 yellow birch
    trees = [
        # species, dbh_in, ht_ft, cr%, tpa
        (12,  6.2, 38, 45, 30.0), (12, 9.1, 52, 40, 18.0),
        (97,  7.8, 46, 50, 24.0), (97, 12.3, 64, 45, 10.0),
        (316, 5.4, 34, 55, 28.0), (316, 10.7, 58, 50, 12.0),
        (318, 8.9, 55, 60, 14.0), (371, 11.5, 61, 48, 9.0),
    ]
    rows = []
    for i, (sp, dbh, ht, cr, tpa) in enumerate(trees, start=1):
        rows.append({"stand_id": "S_SMOKE", "plot_id": 1, "tree_id": i,
                     "tree_count": tpa, "species": sp, "diameter": dbh,
                     "ht": ht, "crratio": cr})
    return plot_data, pd.DataFrame(rows)


def maybe_load_sf(variant):
    """Confirm the conus_sf config block is loadable (Stage-2 readiness)."""
    try:
        from config.config_loader import FvsConfigLoader
        for ver in ("conus_sf", "conus_hybrid"):
            try:
                L = FvsConfigLoader(variant.lower(), version=ver, config_dir=P.CONFIG_DIR)
                if getattr(L, "has_conus_sf_block", False):
                    comps = L.conus_sf_components_present()
                    return {"version": ver, "components": comps, "loader": L}
            except Exception as e:
                last = str(e)
        return {"version": None, "error": locals().get("last", "no conus_sf block")}
    except Exception as e:
        return {"version": None, "error": f"loader import failed: {e}"}


def run_shadow(variant, plot_data, tree_df, num_cycles=10, cycle_length=5,
               config_version="calibrated", out_csv=None):
    """Boot FVS, loop at stop point 5, log engine increments per cycle (no SET)."""
    from fvs2py import FVS
    sdf = build_standinit_peracre(plot_data, "S_SMOKE", variant.lower())

    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "FVS_Data.db")
        conn = sqlite3.connect(db)
        sdf.to_sql("fvs_standinit", conn, if_exists="replace", index=False)
        tree_df.to_sql("fvs_treeinit", conn, if_exists="replace", index=False)
        conn.close()

        cal_kw = "** DEFAULT PARAMETERS"
        if config_version == "calibrated":
            try:
                from config.config_loader import FvsConfigLoader
                cal_kw = FvsConfigLoader(variant.lower(), version="calibrated",
                                         config_dir=P.CONFIG_DIR
                                         ).generate_keywords(include_comments=False)
            except Exception as e:
                print(f"[warn] calibrated keywords unavailable: {e}")

        key = P.KEYFILE_TEMPLATE.format(stand_id="S_SMOKE", db_path=db,
                                        calibration_keywords=cal_kw,
                                        num_cycles=num_cycles,
                                        cycle_length=cycle_length)
        kpath = os.path.join(tmp, f"{variant}_smoke.key")
        with open(kpath, "w") as f:
            f.write(key)

        lib = os.path.join(P.FVS_LIB_DIR, f"FVS{variant.lower()}.so")
        if not os.path.exists(lib):
            raise FileNotFoundError(f"variant lib not found: {lib}")
        print(f"[boot] FVS lib={lib}  config={config_version}")

        fvs = FVS(lib_path=lib, config_version=config_version, config_dir=P.CONFIG_DIR)
        fvs.load_keyfile(kpath)

        log = []
        stop = 0
        fvs.run(stop_point_code=5, stop_point_year=-1)
        while getattr(fvs, "restart_code", 100) == 5:
            stop += 1
            def ga(name):
                try:
                    return np.asarray(fvs.get_tree_attr(name), dtype=float)
                except Exception:
                    return np.array([])
            dbh = ga("dbh"); dg = ga("dg"); htg = ga("htg")
            ht = ga("ht"); spp = ga("species"); tpa = ga("tpa")
            n = len(dbh)
            row = {"stop": stop, "ntrees": n,
                   "mean_dbh_in": float(np.nanmean(dbh)) if n else np.nan,
                   "mean_engine_dg_in": float(np.nanmean(dg)) if len(dg) else np.nan,
                   "mean_engine_htg_ft": float(np.nanmean(htg)) if len(htg) else np.nan,
                   "mean_ht_ft": float(np.nanmean(ht)) if len(ht) else np.nan,
                   "sum_tpa": float(np.nansum(tpa)) if len(tpa) else np.nan}
            log.append(row)
            print(f"[stop {stop:2d}] ntrees={n:3d}  mean dbh={row['mean_dbh_in']:.2f}in  "
                  f"engine dg={row['mean_engine_dg_in']:.4f}in  htg={row['mean_engine_htg_ft']:.3f}ft")
            # ---- SHADOW: Stage-2 hook would compute sf dg/htg here and compare;
            #      no set_tree_attr in Stage 1, so engine behaviour is unchanged.
            fvs.run(stop_point_code=5, stop_point_year=-1)

        # final flush so FVS finishes writing
        try:
            fvs.run()
        except Exception:
            pass

        df = pd.DataFrame(log)
        if out_csv:
            df.to_csv(out_csv, index=False)
            print(f"[done] {stop} stop-point-5 hooks logged -> {out_csv}")
        return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="ne")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="calibrated")
    ap.add_argument("--num-cycles", type=int, default=10)
    ap.add_argument("--cycle-length", type=int, default=5)
    ap.add_argument("--out", default="shadow_log.csv")
    a = ap.parse_args()

    sfinfo = maybe_load_sf(a.variant)
    print(f"[sf-block] variant={a.variant}: {json.dumps({k:v for k,v in sfinfo.items() if k!='loader'})}")

    if not a.smoke:
        print("Only --smoke (synthetic NE stand) is implemented in Stage 1.")
        sys.exit(0)
    plot_data, tdf = synthetic_ne_stand()
    print(f"[stand] synthetic NE: {len(tdf)} trees, sum TPA={tdf['tree_count'].sum():.0f}")
    run_shadow(a.variant, plot_data, tdf, num_cycles=a.num_cycles,
               cycle_length=a.cycle_length, config_version=a.config, out_csv=a.out)


if __name__ == "__main__":
    main()
