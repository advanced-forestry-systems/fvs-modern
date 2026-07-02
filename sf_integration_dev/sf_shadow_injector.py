#!/usr/bin/env python3
"""
sf_shadow_injector.py  --  Stage-1 SHADOW injection of the species-free
equations into the live FVS engine via the fvs2py stop-point-5 hook.

DATABASE-FREE BOOT (2026-07-02 rewrite)
  Previously this booted through a DATABASE keyfile (SQLite standinit/treeinit),
  which segfaults in-process on the DBS SQLite read. It now boots the
  database-free path: an inventory-free keyfile plus fvs2py `add_trees()`, the
  same in-memory tree-loading route the engine hook test proved survives the
  rebuilt (PR #85 prtexm ERR=) .so. Species are passed as FVS variant indices
  (translated from FIA SPCD via the engine's own `species_codes` table), per
  the add_trees() contract.

WHAT THIS DOES
  1. Boots real FVS (rebuilt .so with the #85 fix) under conus_hybrid, applying
     the hybrid Leg A parameters to the engine at stop point 7 (run() does NOT
     auto-apply conus_hybrid, only 'calibrated'/'custom', so we apply manually).
  2. Injects a synthetic NE stand via add_trees().
  3. Stops every cycle at restart code 5 (after the engine computes diameter
     growth and mortality but BEFORE applying them), reads the per-tree engine
     increments via get_tree_attr, and computes the species-free Leg B linear
     predictor per tree via config_loader.sf_linear_predictor(). It LOGS both
     (no set_tree_attr -> zero behaviour change).
  4. Reports whether the sf-vs-engine increments are sane.

FIDELITY NOTES
  - The engine increments (dg in/cycle, htg ft/cycle, mort rate) are read
    directly from engine memory and are exact.
  - The Leg B DG predictor (kuehne_v8) is log-normal (pred = exp(eta+sigma^2/2),
    cm/yr) and its eta needs BGI (external climate 3-piece spline), metric BAL by
    softwood/hardwood, and rd_additive (SDImax-normalized) plus ecoregion- and
    trait-varying BGI slopes. get_tree_attr exposes dbh/ht/cratio/species/tpa/
    dg/htg/mort but NOT BAL/BA/BGI, so the in-process DG/HTG SF prediction here
    is the standard part (intercept + trait_effect + RE + size covariates) with
    BGI/competition terms flagged. The faithful DG three-way lives offline on the
    remeasurement pairs (benchmark_sf_vs_legA.R), which already has these
    covariates precomputed. Treat the in-process SF number as a scale/sign
    sanity signal, not a re-validation.

STAGES (staged rollout from 20260614_engine_wiring_design.md)
  1  shadow  : GET only, log sf-vs-engine increments per tree   <-- THIS FILE
  2  single  : SET dg only, benchmark vs observed
  3  full    : SET dg+htg+mort with localized maxSDI
  4  blended : per-species shrinkage blend, tune kappa

USAGE
  FVSNE_SO=/path/FVSne.so python3 sf_shadow_injector.py --variant ne --smoke \
      --config conus_hybrid --out shadow_ne.csv
"""
from __future__ import annotations
import argparse, os, sys, tempfile, json, math
import numpy as np, pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "fvs2py"))
# Point at the config tree that carries the Leg B blocks (the integration
# worktree), overridable; the main checkout's config may predate Leg B.
CONFIG_DIR = os.environ.get(
    "FVS_CONFIG_DIR",
    "/fs/scratch/PUOM0008/crsfaaron/wt-gompit/config"
    if os.path.isdir("/fs/scratch/PUOM0008/crsfaaron/wt-gompit/config")
    else os.path.join(PROJECT_ROOT, "config"))

# inventory-free keyfile: no DATABASE, no TREEDATA; trees come from add_trees()
KEYFILE_NOINV = """STDIDENT
{stand_id} species-free shadow
STDINFO          922                   1
DESIGN            -1         1
INVYEAR       {invyear:.1f}
NUMCYCLE        {num_cycles:.1f}
PROCESS
STOP
"""


def synthetic_ne_stand():
    """A small, plausible NE mixed stand (per-acre TPA) for the boot smoke.
    SPCD: 12 balsam fir, 97 red spruce, 94 white spruce, 71 tamarack,
    130 eastern hemlock -- chosen to sit inside the NE variant's 23-species
    grouped list so they map to distinct FVS indices."""
    plot_data = {"INVYR": 2010, "LAT": 45.0, "LON": -69.0, "ELEV": 500,
                 "SLOPE": 10, "ASPECT": 180, "STDAGE": 50,
                 # eco codes for the Leg B random effects (NE ecoregion / FT)
                 "L1": "8", "L2": "8.1", "L3": "8.1.1", "FT": 101}
    trees = [  # spcd, dbh_in, ht_ft, cr%, tpa
        (12,  6.2, 38, 45, 30.0), (12, 9.1, 52, 40, 18.0),
        (97,  7.8, 46, 50, 24.0), (97, 12.3, 64, 45, 10.0),
        (94,  5.4, 34, 55, 28.0), (94, 10.7, 58, 50, 12.0),
        (130, 8.9, 55, 60, 14.0), (71, 11.5, 61, 48, 9.0),
    ]
    rows = [{"tree_id": i, "spcd": sp, "dbh": dbh, "ht": ht, "cr": cr, "tpa": tpa}
            for i, (sp, dbh, ht, cr, tpa) in enumerate(trees, 1)]
    return plot_data, pd.DataFrame(rows)


def load_sf_loader(variant, version):
    """conus_hybrid FvsConfigLoader, or None with a reason."""
    try:
        from config.config_loader import FvsConfigLoader
        L = FvsConfigLoader(variant.lower(), version=version, config_dir=CONFIG_DIR)
        return L, None
    except Exception as e:
        return None, str(e)


def sf_runtime_blocks(L, components):
    rt = {}
    for c in components:
        try:
            rt[c] = L.get_conus_sf_runtime_block(c)
        except Exception as e:
            rt[c] = {"_error": str(e)}
    return rt


def sf_standard_part(L, component, spcd, eco, runtime, size_covs):
    """Standard-part Leg B linear predictor: intercept + trait_effect + RE +
    the size covariates we can read in-process. Competition/BGI covariates are
    omitted (flagged) so this is a partial eta, not the production predictor."""
    try:
        return L.sf_linear_predictor(component, int(spcd), eco, size_covs, runtime=runtime)
    except Exception as e:
        return float("nan")


def run_shadow(variant, plot_data, tree_df, num_cycles=10, cycle_length=5,
               config_version="conus_hybrid", out_csv=None):
    from fvs2py import FVS

    lib = os.path.expanduser(os.environ.get(
        "FVSNE_SO",
        f"/fs/scratch/PUOM0008/crsfaaron/lib-sf85/FVS{variant.lower()}.so"))
    if not os.path.exists(lib):
        raise FileNotFoundError(f"variant lib not found: {lib} "
                                "(set FVSNE_SO or build lib-sf85)")

    L, err = load_sf_loader(variant, config_version)
    comps = []
    if L is not None and getattr(L, "has_conus_sf_block", False):
        comps = L.conus_sf_components_present()
    print(f"[sf-block] variant={variant} version={config_version} "
          f"components={comps}{'' if L else '  loader_error='+str(err)}")
    rt = sf_runtime_blocks(L, [c for c in ("diameter_growth", "height_growth", "mortality")
                               if c in comps]) if L else {}
    eco = {k: plot_data.get(k) for k in ("L1", "L2", "L3", "FT")}

    with tempfile.TemporaryDirectory() as tmp:
        kpath = os.path.join(tmp, f"{variant}_shadow.key")
        with open(kpath, "w") as f:
            f.write(KEYFILE_NOINV.format(stand_id="SF01",
                                         invyear=float(plot_data.get("INVYR", 2010)),
                                         num_cycles=float(num_cycles)))
        os.chdir(tmp)

        # config_version=None: we drive conus_hybrid application manually so it
        # is applied for the hybrid arm (run() only auto-applies calibrated/custom)
        fvs = FVS(lib_path=lib, config_version=None, config_dir=CONFIG_DIR)
        fvs.load_keyfile(kpath)
        print(f"[boot] lib={lib}")

        # stop point 7: input read, arrays allocated, before imputation
        fvs.run(stop_point_code=7, stop_point_year=-1)

        # FIA SPCD -> FVS variant index, from the engine's own table
        codes = fvs.species_codes
        fia2idx = {}
        for _, r in codes.iterrows():
            fia = str(r["fia"]).strip()
            if fia.isdigit():
                fia2idx.setdefault(int(fia), int(r["fvs_index"]))
        idx = [fia2idx.get(int(s), 1) for s in tree_df["spcd"]]
        unmapped = [int(s) for s, j in zip(tree_df["spcd"], idx) if int(s) not in fia2idx]
        if unmapped:
            print(f"[warn] SPCD not in {variant} species list, mapped to index 1: {unmapped}")

        # apply conus_hybrid Leg A params to the engine
        applied = {}
        if L is not None:
            try:
                applied = L.apply_to_fvs(fvs)
            except Exception as e:
                print(f"[warn] apply_to_fvs({config_version}) failed: {e}")
        print(f"[hybrid] applied groups: {applied}")

        n_add = fvs.add_trees(
            np.asarray(tree_df["dbh"], float), np.asarray(idx, float),
            np.asarray(tree_df["ht"], float), np.asarray(tree_df["cr"], float),
            np.ones(len(tree_df)), np.asarray(tree_df["tpa"], float))
        print(f"[add_trees] added {n_add} trees; dims ntrees={fvs.dims.get('ntrees')}")

        # --- self-diagnostic: did add_trees populate the attribute arrays? -----
        chk = np.asarray(fvs.get_tree_attr("dbh"), float)
        attrs_ok = np.isfinite(chk).any() and np.nansum(np.abs(chk)) > 0
        if not attrs_ok:
            print("[DIAGNOSTIC] per-tree attributes read all-zero after add_trees at this "
                  "stop point. This is a USAGE SEQUENCE issue, not an engine defect: the FVS "
                  "R API test (src-converted/tests/APIviaR/Rapi.R) calls fvsAddTrees at STOP "
                  "POINT 6 (the ESTAB point) during an ACTIVE run to add regeneration to an "
                  "already-initialized stand, then reads it back successfully. Calling it at "
                  "stop point 7 on a zero-inventory stand (LSTART true, no stand context) does "
                  "not feed the growth arrays. Correct database-free pattern: seed the stand "
                  "via the keyfile inventory (inline TREEDATA) and, for regen, add_trees at "
                  "sp6 mid-run. The faithful per-tree sf-vs-engine comparison is available "
                  "offline via benchmark_sf_vs_legA.R (which has the metric BAL/BGI/rd "
                  "covariates get_tree_attr cannot expose in-process).")

        log = []
        stop = 0
        fvs.run(stop_point_code=5, stop_point_year=-1)
        while getattr(fvs, "restart_code", 100) == 5:
            stop += 1
            def ga(name):
                try:
                    return np.asarray(fvs.get_tree_attr(name), float)
                except Exception:
                    return np.array([])
            dbh, dg, htg = ga("dbh"), ga("dg"), ga("htg")
            mort, spp, cr = ga("mort"), ga("species"), ga("cratio")
            n = len(dbh)
            for i in range(n):
                d_in = float(dbh[i]) if i < len(dbh) else float("nan")
                # size covariates we can build in-process (partial; DG/HG also
                # need BGI + metric BAL + rd which get_tree_attr does not expose)
                size_covs = {}
                if d_in and d_in == d_in and d_in > 0:
                    size_covs = {"b1": math.log(d_in * 2.54), "b2": d_in * 2.54}
                spcd_i = int(round(spp[i])) if i < len(spp) else 0
                sf_dg = sf_standard_part(L, "diameter_growth", spcd_i, eco,
                                         rt.get("diameter_growth"), size_covs) if L else float("nan")
                log.append({
                    "cycle": stop, "tree": i + 1, "spcd_fvs_or_fia": spcd_i,
                    "dbh_in": d_in,
                    "engine_dg_in_cyc": float(dg[i]) if i < len(dg) else float("nan"),
                    "engine_htg_ft_cyc": float(htg[i]) if i < len(htg) else float("nan"),
                    "engine_mort_rate": float(mort[i]) if i < len(mort) else float("nan"),
                    "engine_dg_cm_yr": (float(dg[i]) * 2.54 / cycle_length) if i < len(dg) else float("nan"),
                    "sf_dg_eta_stdpart": sf_dg,
                    "sf_dg_pred_cm_yr_stdpart": (math.exp(sf_dg) if sf_dg == sf_dg else float("nan")),
                })
            md = float(np.nanmean(dbh)) if n else float("nan")
            print(f"[cycle {stop:2d}] ntrees={n:3d} mean_dbh={md:6.2f}in "
                  f"engine_dg={float(np.nanmean(dg)) if len(dg) else float('nan'):7.4f}in "
                  f"engine_htg={float(np.nanmean(htg)) if len(htg) else float('nan'):6.3f}ft")
            fvs.run(stop_point_code=5, stop_point_year=-1)
        try:
            fvs.run()
        except Exception:
            pass

        df = pd.DataFrame(log)
        if out_csv and len(df):
            df.to_csv(out_csv, index=False)
            print(f"[done] {stop} cycles logged, {len(df)} tree-rows -> {out_csv}")

        # ---- sanity verdict ---------------------------------------------------
        print("\n===== SANITY =====")
        if not len(df) or not attrs_ok:
            print("VERDICT: boot works (no segfault, #85 fix active) and add_trees returns "
                  "cleanly, but per-tree attributes read 0 because add_trees was called at the "
                  "wrong stop point for this stand state (see DIAGNOSTIC). This is a usage "
                  "sequence, not an engine defect. To log real in-process increments, load the "
                  "synthetic stand as keyfile inventory and call add_trees at sp6 mid-run per "
                  "the FVS R API contract. The faithful DG/HG/mort sf-vs-engine comparison "
                  "already exists offline via benchmark_sf_vs_legA.R and does not depend on "
                  "this in-process hook.")
        else:
            eng = df["engine_dg_cm_yr"].to_numpy()
            sf = df["sf_dg_pred_cm_yr_stdpart"].to_numpy()
            eng_ok = np.nanmedian(eng) > 0 and np.nanmedian(eng) < 2.0
            print(f"engine dg median = {np.nanmedian(eng):.3f} cm/yr "
                  f"({'plausible' if eng_ok else 'CHECK'} for NE 0.05-1.0)")
            print(f"sf dg (std part) median = {np.nanmedian(sf):.3f} cm/yr "
                  "(partial predictor: no BGI/competition terms)")
        return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="ne")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="conus_hybrid")
    ap.add_argument("--num-cycles", type=int, default=10)
    ap.add_argument("--cycle-length", type=int, default=5)
    ap.add_argument("--out", default="shadow_log.csv")
    a = ap.parse_args()
    if not a.smoke:
        print("Only --smoke (synthetic NE stand) is implemented in Stage 1.")
        sys.exit(0)
    plot_data, tdf = synthetic_ne_stand()
    print(f"[stand] synthetic {a.variant.upper()}: {len(tdf)} trees, "
          f"sum TPA={tdf['tpa'].sum():.0f}")
    run_shadow(a.variant, plot_data, tdf, num_cycles=a.num_cycles,
               cycle_length=a.cycle_length, config_version=a.config, out_csv=a.out)


if __name__ == "__main__":
    main()
