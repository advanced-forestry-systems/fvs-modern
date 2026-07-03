#!/usr/bin/env python3
"""
arm_uncertainty.py -- full posterior uncertainty at BOTH tree and stand level
for the three CONUS arms (conus / conus_sf / conus_organon) plus the shared
Bayesian modifier layer.

The methodological core: trees within a stand SHARE the drawn global fixed
effects, the drawn species / trait effects, the drawn ecoregion / forest-type
random effects, and the drawn modifier coefficients. Only the residual is drawn
independently per tree. That shared structure induces the correct positive
correlation among trees in a stand, so the stand-level predictive interval
reflects it instead of the too-narrow interval a naive independent-tree
aggregation would give.

Monte Carlo, per draw d:
  theta_d      ~ N(fixed_mean, fixed_sd)           shared across the stand
  z_level_d    ~ N(re_mean_level, sigma_level)      shared across the stand (per L1/L2/L3/FT)
  sp_effect_d  ~ N(sp_mean, sp_sd)                  shared across trees of that species
  mod_coef_d   ~ N(mod_mean, mod_sd)                shared across the stand
  eps_id       ~ N(0, sigma_resid)                  independent per tree
  y_id = backtransform(eta_id + eps_id) * exp(mod_eta_id)

Tree-level interval  = quantiles of y_id over d.
Stand-level interval = quantiles of aggregate_d = sum_i (y_id contribution) over d.

Reference implementation for diameter growth (dg). The same machinery applies to
the other components (each stores the same posterior pieces) and to conus /
conus_organon (per-species intercept posterior / organon coefficient posterior).
"""
from __future__ import annotations
import os, sys, math
import numpy as np

sys.path.insert(0, os.environ.get("FVS_WT", "/fs/scratch/PUOM0008/crsfaaron/wt-gompit"))
from config.config_loader import FvsConfigLoader

RNG = np.random.default_rng(20260702)


def _fixed_draws(rt_fixed_mean, rt_fixed_sd, params, n_draws):
    """Draw shared global fixed effects: {param: array(n_draws)}."""
    out = {}
    for p in params:
        m = rt_fixed_mean.get(p, 0.0); s = rt_fixed_sd.get(p, 0.0)
        out[p] = RNG.normal(m, max(s, 0.0), n_draws) if s > 0 else np.full(n_draws, m)
    return out


def dg_uncertainty_leg_b(L, stand_trees, eco, drivers=None, n_draws=1000):
    """Full tree+stand uncertainty for the Leg B (conus_sf) DG arm + modifier.

    stand_trees: list of dicts with spcd and the dg covariates already built
                 (ln_dbh, DBH1(cm), ln_cr_adj, ln_bal_sw_adj, bal_hw, bgi, ...
                  matching the kuehne_v8 fixed-effect param names).
    Returns dict with tree_pi (per tree q05/q50/q95 of dg cm/yr) and stand_pi
    (q05/q50/q95 of stand mean dg), plus the naive-independent stand SD for
    contrast.
    """
    b = L.get_conus_sf_block("diameter_growth")
    fe = b["fixed_effects"]
    fmean = dict(zip(fe["param"], fe["mean"]))
    fsd = dict(zip(fe["param"], fe.get("sd", [0.0] * len(fe["param"]))))
    params = [p for p in fe["param"] if not p.startswith(("sigma", "phi"))]
    sigma_resid = float(fmean.get("sigma", 1.0))
    # variance components for shared RE draws
    sig = {lvl: float(fmean.get(f"sigma_{lvl}", 0.0)) for lvl in ("L1", "L2", "L3", "FT")}
    rt = L.get_conus_sf_runtime_block("diameter_growth")

    theta = _fixed_draws(fmean, fsd, params, n_draws)          # shared fixed effects
    # shared ecoregion / FT random effects (drawn once per stand per draw)
    z = {}
    for lvl in ("L1", "L2", "L3", "FT"):
        code = (eco or {}).get(lvl)
        re_mean = rt.get(f"re_{lvl}", {}).get(str(code), 0.0)
        z[lvl] = RNG.normal(re_mean, sig[lvl], n_draws) if sig[lvl] > 0 else np.full(n_draws, re_mean)
    # shared modifier coefficient draws (mean +/- sd if present, else fixed point)
    mrt = L.get_modifier_runtime_block("diameter_growth")
    drivers = drivers or {}
    mblk = (L.config.get("categories_conus_mod", {}) or {}).get("diameter_growth", {})
    def _mdraw(mean, sd):
        return RNG.normal(mean, sd, n_draws) if (sd and sd > 0) else np.full(n_draws, mean)
    if mrt.get("_present"):
        cm = _mdraw(mrt["coef_mgmt"], float((mblk.get("management") or {}).get("coef_sd", 0.0)))
        cd = _mdraw(mrt["coef_dstrb"], float((mblk.get("disturbance") or {}).get("coef_sd", 0.0)))
        cb1 = _mdraw(mrt["b_bgi1"], float((mblk.get("driver_bgi") or {}).get("b1_sd", 0.0)))
        cb2 = _mdraw(mrt["b_bgi2"], float((mblk.get("driver_bgi") or {}).get("b2_sd", 0.0)))

    # per-species trait-effect draws (shared across trees of a species): use the
    # bundle trait_effect mean with sigma_sp as the shared species-level scale
    sigma_sp = float(fmean.get("sigma_sp", 0.0))

    def sp_effect(spcd):
        te = L.sf_trait_effect(rt, int(spcd))
        return RNG.normal(te, sigma_sp, n_draws) if sigma_sp > 0 else np.full(n_draws, te)

    intercept = theta.get(rt["intercept_name"], np.full(n_draws, rt["intercept"]))
    tree_draws = np.zeros((len(stand_trees), n_draws))
    for i, t in enumerate(stand_trees):
        eta = intercept + sp_effect(t["spcd"]) + z["L1"] + z["L2"] + z["L3"] + z["FT"]
        for p in params:
            if p == rt["intercept_name"]:
                continue
            eta = eta + theta[p] * float(t.get(p, 0.0))
        eps = RNG.normal(0.0, sigma_resid, n_draws)             # independent per tree
        dg = np.exp(np.clip(eta + eps, -30, 20))                # log-normal backtransform, cm/yr
        # shared modifier multiplier per draw
        mult = np.ones(n_draws)
        if mrt.get("_present"):
            meta = np.zeros(n_draws)
            yst = drivers.get("years_since_trt"); ysd = drivers.get("years_since_dstrb"); bgi = drivers.get("bgi")
            if yst is not None and yst >= 0:
                meta = meta + cm * math.exp(-min(yst, 100) / mrt["tau_m"])
            if ysd is not None and ysd >= 0:
                meta = meta + cd * math.exp(-min(ysd, 100) / mrt["tau_d"])
            if bgi is not None:
                meta = meta + cb1 * bgi + cb2 * max(bgi - mrt["bgi_knot"], 0.0)
            mult = np.exp(meta)
        tree_draws[i] = dg * mult

    q = lambda a, p: float(np.quantile(a, p))
    tree_pi = [{"spcd": t["spcd"], "q05": q(tree_draws[i], .05),
                "q50": q(tree_draws[i], .5), "q95": q(tree_draws[i], .95)}
               for i, t in enumerate(stand_trees)]
    stand_by_draw = tree_draws.mean(axis=0)                     # shared-draw stand mean dg
    # naive independent contrast: shuffle each tree's draws independently, re-aggregate
    naive = np.array([RNG.permutation(tree_draws[i]) for i in range(len(stand_trees))]).mean(axis=0)
    return {
        "tree_pi": tree_pi,
        "stand_pi": {"q05": q(stand_by_draw, .05), "q50": q(stand_by_draw, .5), "q95": q(stand_by_draw, .95),
                     "sd": float(stand_by_draw.std())},
        "stand_naive_sd": float(naive.std()),
        "n_draws": n_draws,
    }


if __name__ == "__main__":
    CFG = os.environ.get("FVS_CONFIG_DIR", "/fs/scratch/PUOM0008/crsfaaron/wt-gompit/config")
    L = FvsConfigLoader("ne", version="conus_sf", config_dir=CFG)
    # synthetic NE stand: build kuehne_v8 dg covariates per tree
    def tree(spcd, dbh_cm, cr, bal_sw, bal_hw, bgi=6.0):
        return {"spcd": spcd, "b1": math.log(dbh_cm), "b2": dbh_cm,
                "b3": math.log((cr + 0.2) / 1.2), "b4": math.log(bal_sw + 0.01),
                "b5": bal_hw, "b6": bgi}
    stand = [tree(12, 15, .5, 5, 2), tree(97, 20, .45, 8, 3),
             tree(12, 25, .4, 12, 4), tree(97, 30, .5, 15, 5)]
    eco = {"L1": "8", "L2": "8.1", "L3": "8.1.1", "FT": 101}
    res = dg_uncertainty_leg_b(L, stand, eco, drivers={"bgi": 6.0}, n_draws=2000)
    print("== DG full uncertainty (conus_sf), synthetic NE stand ==")
    for t in res["tree_pi"]:
        print(f"  tree spcd {t['spcd']:3d}: dg cm/yr  q05={t['q05']:.3f}  q50={t['q50']:.3f}  q95={t['q95']:.3f}")
    sp = res["stand_pi"]
    print(f"  STAND mean dg: q05={sp['q05']:.3f} q50={sp['q50']:.3f} q95={sp['q95']:.3f}  sd={sp['sd']:.4f}")
    print(f"  stand SD with shared-RE correlation: {sp['sd']:.4f}")
    print(f"  stand SD naive independent (WRONG):  {res['stand_naive_sd']:.4f}")
    print(f"  correlation inflation factor: {sp['sd']/res['stand_naive_sd']:.2f}x")
    print("UNCERTAINTY_OK")
