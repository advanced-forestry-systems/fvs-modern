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


def _modifier_draws(L, component, n_draws):
    """Shared modifier-coefficient draws for a stand (one draw set per stand,
    shared across every tree in it). Returns (mrt, mult_fn) where mult_fn(drivers)
    -> array(n_draws) multiplier. Coefs drawn ~ N(mean, sd) when the bundle
    carries coef_sd, else held at the point estimate. Identical across arms so
    every arm applies the SAME shared multiplicative modifier layer."""
    mrt = L.get_modifier_runtime_block(component)
    mblk = (L.config.get("categories_conus_mod", {}) or {}).get(component, {})

    def _mdraw(mean, sd):
        return RNG.normal(mean, sd, n_draws) if (sd and sd > 0) else np.full(n_draws, mean)

    if mrt.get("_present"):
        cm = _mdraw(mrt["coef_mgmt"], float((mblk.get("management") or {}).get("coef_sd", 0.0)))
        cd = _mdraw(mrt["coef_dstrb"], float((mblk.get("disturbance") or {}).get("coef_sd", 0.0)))
        cb1 = _mdraw(mrt["b_bgi1"], float((mblk.get("driver_bgi") or {}).get("b1_sd", 0.0)))
        cb2 = _mdraw(mrt["b_bgi2"], float((mblk.get("driver_bgi") or {}).get("b2_sd", 0.0)))
    else:
        cm = cd = cb1 = cb2 = np.zeros(n_draws)

    def mult_fn(drivers):
        if not mrt.get("_present"):
            return np.ones(n_draws)
        drivers = drivers or {}
        meta = np.zeros(n_draws)
        yst = drivers.get("years_since_trt"); ysd = drivers.get("years_since_dstrb"); bgi = drivers.get("bgi")
        if yst is not None and yst >= 0:
            meta = meta + cm * math.exp(-min(yst, 100) / mrt["tau_m"])
        if ysd is not None and ysd >= 0:
            meta = meta + cd * math.exp(-min(ysd, 100) / mrt["tau_d"])
        if bgi is not None:
            meta = meta + cb1 * bgi + cb2 * max(bgi - mrt["bgi_knot"], 0.0)
        return np.exp(meta)

    return mrt, mult_fn


def _aggregate(stand_trees, tree_draws, n_draws):
    """Shared-draw stand aggregation + naive-independent contrast.
    Returns the tree_pi / stand_pi / stand_naive_sd structure common to all
    three arms."""
    q = lambda a, p: float(np.quantile(a, p))
    tree_pi = [{"spcd": t["spcd"], "q05": q(tree_draws[i], .05),
                "q50": q(tree_draws[i], .5), "q95": q(tree_draws[i], .95)}
               for i, t in enumerate(stand_trees)]
    stand_by_draw = tree_draws.mean(axis=0)                     # shared-draw stand mean
    # naive independent contrast: shuffle each tree's draws independently, re-aggregate
    naive = np.array([RNG.permutation(tree_draws[i]) for i in range(len(stand_trees))]).mean(axis=0)
    return {
        "tree_pi": tree_pi,
        "stand_pi": {"q05": q(stand_by_draw, .05), "q50": q(stand_by_draw, .5),
                     "q95": q(stand_by_draw, .95), "sd": float(stand_by_draw.std())},
        "stand_naive_sd": float(naive.std()),
        "n_draws": n_draws,
    }


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


def dg_uncertainty_leg_a(L, stand_trees, eco, drivers=None, n_draws=1000,
                         sigma_resid=None):
    """Full tree+stand uncertainty for the Leg A (conus) DG arm + modifier.

    The Leg A DG block is categories_conus.diameter_growth. Its per-species
    intercepts (species_intercepts: SPCD, mean, sd) are the ONLY drawn quantity
    here: each species intercept ~ N(mean, sd), SHARED across all trees of that
    species in the stand. The FVS-native base linear predictor (global fixed
    effects b0..b8 over the tree covariates + the weighted ecodivision RE) is
    treated as fixed/known, matching the Leg A operational convention. The
    residual sigma is per-tree independent.

    stand_trees: list of dicts with 'spcd' and the Leg A DG covariates keyed by
                 the fixed-effect param names (b1..b8), e.g. b1=ln_dbh, ...
    Returns the same tree_pi / stand_pi / stand_naive_sd structure as Leg B.
    """
    block = L.get_conus_block("diameter_growth")           # raw block (has species_intercepts sd)
    rt = L.get_conus_runtime_block("diameter_growth")      # fixed / species_re / eco_re
    fixed = rt["fixed"]
    eco_re = float(rt.get("eco_re", 0.0))
    intercept = float(fixed.get("b0", 0.0))
    # residual sigma: Leg A carries its own 'sigma'; use it, else fall back to
    # the Leg B sf sigma as a tagged stand-in.
    resid_tag = "leg_a_sigma"
    if sigma_resid is None:
        sigma_resid = float(fixed.get("sigma", 0.0))
        if not (sigma_resid > 0):
            try:
                sf = L.get_conus_sf_block("diameter_growth")
                sfe = sf["fixed_effects"]
                sf_fixed = dict(zip(sfe["param"], sfe["mean"]))
                sigma_resid = float(sf_fixed.get("sigma", 1.0)); resid_tag = "leg_b_sigma_standin"
            except Exception:
                sigma_resid = 1.0; resid_tag = "default_sigma_standin"

    # per-species shared intercept draws ~ N(mean, sd) from species_intercepts
    si = block.get("species_intercepts", {}) or {}
    si_spcd = [int(s) for s in (si.get("SPCD") or [])]
    si_mean = {s: float(m) for s, m in zip(si_spcd, si.get("mean", []))}
    si_sd = {s: float(v) for s, v in zip(si_spcd, si.get("sd", []))}

    def sp_intercept(spcd):
        spcd = int(spcd)
        m = si_mean.get(spcd, 0.0)              # species outside Leg A -> 0 (base only)
        s = si_sd.get(spcd, 0.0)
        return RNG.normal(m, s, n_draws) if s > 0 else np.full(n_draws, m)

    # covariate params = fixed b1..bK (exclude intercept + scale/sigma params)
    cov_params = [p for p in fixed if p != "b0" and not p.startswith(("sigma", "phi", "lp__"))]
    _mrt, mult_fn = _modifier_draws(L, "diameter_growth", n_draws)

    tree_draws = np.zeros((len(stand_trees), n_draws))
    for i, t in enumerate(stand_trees):
        base_cov = sum(float(fixed[p]) * float(t.get(p, 0.0)) for p in cov_params)  # fixed/known base
        eta = intercept + eco_re + base_cov + sp_intercept(t["spcd"])              # shared sp draw
        eps = RNG.normal(0.0, sigma_resid, n_draws)                                # per-tree residual
        dg = np.exp(np.clip(eta + eps, -30, 20))                                   # log-normal cm/yr
        tree_draws[i] = dg * mult_fn(drivers)                                      # shared modifier

    out = _aggregate(stand_trees, tree_draws, n_draws)
    out["arm"] = "conus (Leg A)"
    out["resid_source"] = resid_tag
    out["note"] = ("Leg A DG: base linear predictor treated as fixed/known; "
                   "only species intercepts drawn ~ N(mean, sd) shared across "
                   "trees of a species; residual per-tree independent.")
    return out


def dg_uncertainty_organon(L, stand_trees, eco, drivers=None, n_draws=1000,
                           nominal_coef_sd=0.05, sigma_resid=None):
    """Full tree+stand uncertainty for the conus_organon arm + modifier.

    Uses categories_conus_organon.diameter_growth (organon_linear_predictor):
    global fixed_effects (a0..a6) with a per-species calibration multiplier.
    The scaffold carries NO posterior SDs yet, so the fixed coefficients are
    drawn ~ N(mean, nominal_coef_sd*|mean|) SHARED across the stand -- clearly
    TAGGED as a placeholder that BACKFILLS from the CONUS ORGANON recalibration
    once posteriors exist. The per-species calib is applied (as fixed). Same
    shared-draw stand aggregation + shared modifier.

    stand_trees: list of dicts with 'spcd' and the organon covariates keyed by
                 the native param names (a1..a6: ln_dbh, dbh_sq, ln_cr_adj,
                 ln_si, bal_over_lndbh, sqrt_sba).
    Returns the same tree_pi / stand_pi / stand_naive_sd structure.
    """
    rt = L.get_conus_organon_runtime_block("diameter_growth")
    fixed = rt["fixed"]                       # a0..a6 means
    sp_coef = rt["species_coef"]              # {spcd: {calib: ...}}
    cov_params = [p for p in fixed if p != "a0"]

    # shared fixed-coef draws ~ N(mean, nominal sd) -- BACKFILL PLACEHOLDER
    def _coef_draw(p):
        m = float(fixed[p]); s = abs(m) * nominal_coef_sd
        return RNG.normal(m, s, n_draws) if s > 0 else np.full(n_draws, m)
    theta = {p: _coef_draw(p) for p in fixed}

    # residual sigma: organon scaffold has none -> Leg B sf sigma stand-in, tagged
    resid_tag = "leg_b_sigma_standin"
    if sigma_resid is None:
        try:
            sf = L.get_conus_sf_block("diameter_growth")
            sfe = sf["fixed_effects"]
            sf_fixed = dict(zip(sfe["param"], sfe["mean"]))
            sigma_resid = float(sf_fixed.get("sigma", 1.0))
        except Exception:
            sigma_resid = 1.0; resid_tag = "default_sigma_standin"

    _mrt, mult_fn = _modifier_draws(L, "diameter_growth", n_draws)

    tree_draws = np.zeros((len(stand_trees), n_draws))
    for i, t in enumerate(stand_trees):
        calib = float((sp_coef.get(int(t["spcd"]), {}) or {}).get("calib", 1.0))  # per-species calib (fixed)
        eta = theta["a0"] + sum(theta[p] * float(t.get(p, 0.0)) for p in cov_params)
        eta = eta * calib
        eps = RNG.normal(0.0, sigma_resid, n_draws)
        dg = np.exp(np.clip(eta + eps, -30, 20))                                  # log-normal cm/yr
        tree_draws[i] = dg * mult_fn(drivers)                                     # shared modifier

    out = _aggregate(stand_trees, tree_draws, n_draws)
    out["arm"] = "conus_organon"
    out["resid_source"] = resid_tag
    out["coef_posterior"] = (f"PLACEHOLDER: fixed coefs drawn ~ N(mean, "
                             f"{nominal_coef_sd:.0%}*|mean|); BACKFILLS from the "
                             f"CONUS ORGANON recalibration when posteriors land.")
    return out


def stand_constrained_mortality_uncertainty(
        L, stand_trees, eco, T_years=10.0, n_draws=1000,
        drivers=None,
        stand_survival_bundle_path="/fs/scratch/PUOM0008/crsfaaron/"
        "fvs-conus_output_conus/stand_survival/stand_survival_bundle.json",
        sdimax_samples_path=None, fixed_sdimax=None):
    """
    Full-uncertainty stand-constrained MORTALITY path over a synthetic stand.

    Per Monte-Carlo draw d (trees SHARE the drawn stand-level structure, so the
    stand correlation from the shared draw is preserved):

      (a) Density cap: draw SDIMAX from the self-thinning posterior
          (sdimax_draws_from_posterior on stand_density_samples.rds); if the rds
          is unavailable fall back to a fixed SDIMAX (TAGGED). Apply the
          probabilistic Reineke cap (sdimax_density_cap_draws) to the per-tree
          TPA so the density ceiling carries the SDIMAX posterior uncertainty.

      (b) Mortality disaggregation: if stand_survival_bundle.json exists, compute
          the stand mortality target M_stand from it (cloglog + log(YEARS)
          exposure, using rd/ln_qmd/ba_metric/bgi), then call
          stand_disaggregate_mortality to solve kappa and rescale the per-tree
          hazards so the (capped-)TPA-weighted tree mortality equals M_stand.
          If the bundle is absent, skip disaggregation gracefully (TAGGED
          "stand_survival pending"): report the unconstrained tree mortality.

    stand_trees: list of dicts each with 'spcd', 'dbh_cm', 'tpa', 'base_hazard'
                 (per-tree ANNUAL base hazard from the arm's survival model),
                 and the stand summaries carried on the first tree or passed via
                 the returned dict. QMD (cm) per tree via 'dbh_cm'.
    Returns a dict with unconstrained vs constrained stand mortality (mean +
    q05/q95 over draws), the density-cap effect, kappa summary, and tags.
    """
    import json
    import os
    # import the stand-constraint helpers (support running as a script or module)
    try:
        from stand_constraint import (
            sdimax_density_cap_draws, sdimax_draws_from_posterior,
            stand_disaggregate_mortality, stand_mortality_target,
        )
    except Exception:
        from sf_integration_dev.stand_constraint import (  # type: ignore
            sdimax_density_cap_draws, sdimax_draws_from_posterior,
            stand_disaggregate_mortality, stand_mortality_target,
        )

    drivers = drivers or {}
    T = float(T_years)
    dbh = np.array([float(t["dbh_cm"]) for t in stand_trees])
    tpa0 = np.array([float(t["tpa"]) for t in stand_trees])
    h_base = np.array([float(t["base_hazard"]) for t in stand_trees])

    # stand summaries (TPA-weighted) -- QMD from the per-tree DBH, BA metric,
    # relative density rd = SDI/SDIMAX (canonical rd from the design note).
    W = float(np.sum(tpa0))
    qmd = float(np.sqrt(np.sum(tpa0 * dbh ** 2) / W))            # cm
    ba_metric = float(np.sum(tpa0 * (math.pi / 4.0) * (dbh / 100.0) ** 2))  # m2/ha
    ln_qmd = math.log(qmd)
    bgi = float(drivers.get("bgi", 6.0))
    sdi_stand = float(np.sum(tpa0 * (dbh / 10.0) ** 1.605))

    tags = []

    # ---- SDIMAX posterior draws (or fixed fallback, tagged) ------------------
    sdimax_draws = None
    if sdimax_samples_path is not None:
        sdimax_draws = sdimax_draws_from_posterior(sdimax_samples_path,
                                                   n_draws=n_draws, rng=RNG)
    if sdimax_draws is None:
        fx = fixed_sdimax if fixed_sdimax is not None else max(sdi_stand * 1.05, 1.0)
        sdimax_draws = np.full(n_draws, float(fx))
        tags.append(f"SDIMAX fixed fallback={float(fx):.1f} (posterior rds unavailable)")
    else:
        if sdimax_draws.size < n_draws:
            sdimax_draws = RNG.choice(sdimax_draws, size=n_draws, replace=True)
        else:
            sdimax_draws = sdimax_draws[:n_draws]
        tags.append(f"SDIMAX drawn from self-thinning posterior "
                    f"(median={float(np.median(sdimax_draws)):.1f})")

    # ---- density cap per draw (probabilistic Reineke) -----------------------
    cap = sdimax_density_cap_draws(tpa0, dbh, sdimax_draws)     # QMD=dbh here
    capped_tpa = cap["capped_tpa"]                              # (n_draws, n_trees)
    sdi_uncapped = cap["sdi"]
    frac_binds = float(np.mean(cap["binds"]))

    # rd uses the capped SDI vs drawn SDIMAX (density state after the cap)
    sdi_capped_draws = np.sum(capped_tpa * (dbh[None, :] / 10.0) ** 1.605, axis=1)
    rd_draws = sdi_capped_draws / np.maximum(sdimax_draws, 1e-9)

    # ---- stand mortality target from the bundle (if present) ----------------
    bundle = None
    if os.path.exists(stand_survival_bundle_path):
        try:
            with open(stand_survival_bundle_path) as fh:
                bundle = json.load(fh)
        except Exception:
            bundle = None
    have_bundle = bundle is not None

    # per-draw unconstrained (capped-TPA-weighted) and constrained stand mortality
    unconstrained = np.zeros(n_draws)
    constrained = np.zeros(n_draws)
    kappas = np.zeros(n_draws)
    m_target_draws = np.zeros(n_draws)
    for d in range(n_draws):
        w = capped_tpa[d]
        surv = np.exp(-h_base * T)
        m_unc = float(np.sum(w * (1.0 - surv)) / np.sum(w))
        unconstrained[d] = m_unc
        if have_bundle:
            m_tgt = stand_mortality_target(
                bundle, T, rd=float(rd_draws[d]), ln_qmd=ln_qmd,
                ba_metric=ba_metric, bgi=bgi,
                trt_decay=float(drivers.get("trt_decay", 0.0)),
                dstrb_decay=float(drivers.get("dstrb_decay", 0.0)),
                L1=(eco or {}).get("L1"))
            m_target_draws[d] = m_tgt
            kappa, h_prime = stand_disaggregate_mortality(h_base, w, T, m_tgt)
            kappas[d] = kappa
            constrained[d] = float(np.sum(w * (1.0 - np.exp(-h_prime * T)))
                                   / np.sum(w))
        else:
            constrained[d] = m_unc
            kappas[d] = 1.0

    if not have_bundle:
        tags.append("stand_survival pending (bundle absent): mortality "
                    "disaggregation skipped, constrained == unconstrained")

    q = lambda a, p: float(np.quantile(a, p))
    uncapped_total = float(W)
    capped_total_mean = float(capped_tpa.sum(axis=1).mean())
    return {
        "have_bundle": have_bundle,
        "T_years": T,
        "stand_summary": {"qmd_cm": qmd, "ba_m2ha": ba_metric,
                          "sdi_uncapped": sdi_uncapped,
                          "rd_mean": float(rd_draws.mean())},
        "density_cap": {
            "sdi_uncapped": sdi_uncapped,
            "sdimax_median": float(np.median(sdimax_draws)),
            "frac_draws_bind": frac_binds,
            "tpa_uncapped": uncapped_total,
            "tpa_capped_mean": capped_total_mean,
            "tpa_reduction_pct": 100.0 * (1.0 - capped_total_mean / uncapped_total),
        },
        "stand_mortality_unconstrained": {
            "mean": float(unconstrained.mean()),
            "q05": q(unconstrained, .05), "q95": q(unconstrained, .95)},
        "stand_mortality_constrained": {
            "mean": float(constrained.mean()),
            "q05": q(constrained, .05), "q95": q(constrained, .95)},
        "stand_mortality_target": {
            "mean": float(m_target_draws.mean()) if have_bundle else None,
            "q05": q(m_target_draws, .05) if have_bundle else None,
            "q95": q(m_target_draws, .95) if have_bundle else None},
        "kappa": {"mean": float(kappas.mean()), "q05": q(kappas, .05),
                  "q95": q(kappas, .95)},
        "n_draws": n_draws,
        "tags": tags,
    }


def _report_mortality(res):
    print("\n== Stand-constrained MORTALITY (full uncertainty), synthetic NE stand ==")
    ss = res["stand_summary"]
    print(f"  stand: QMD={ss['qmd_cm']:.1f} cm  BA={ss['ba_m2ha']:.2f} m2/ha  "
          f"SDI={ss['sdi_uncapped']:.1f}  rd(mean)={ss['rd_mean']:.3f}")
    dc = res["density_cap"]
    print(f"  DENSITY CAP: SDIMAX(median)={dc['sdimax_median']:.1f}  "
          f"binds in {100*dc['frac_draws_bind']:.0f}% of draws  "
          f"TPA {dc['tpa_uncapped']:.1f} -> {dc['tpa_capped_mean']:.1f} "
          f"({dc['tpa_reduction_pct']:.1f}% mean reduction)")
    mu = res["stand_mortality_unconstrained"]; mc = res["stand_mortality_constrained"]
    print(f"  UNCONSTRAINED stand mortality: mean={mu['mean']:.4f} "
          f"[q05={mu['q05']:.4f}, q95={mu['q95']:.4f}]")
    mt = res["stand_mortality_target"]
    if res["have_bundle"]:
        print(f"  STAND TARGET   (from bundle):  mean={mt['mean']:.4f} "
              f"[q05={mt['q05']:.4f}, q95={mt['q95']:.4f}]")
        print(f"  CONSTRAINED stand mortality:   mean={mc['mean']:.4f} "
              f"[q05={mc['q05']:.4f}, q95={mc['q95']:.4f}]  "
              f"(matches target: {abs(mc['mean']-mt['mean']):.2e})")
        kp = res["kappa"]
        print(f"  kappa (proportional-hazard): mean={kp['mean']:.4f} "
              f"[q05={kp['q05']:.4f}, q95={kp['q95']:.4f}]")
    else:
        print(f"  CONSTRAINED == UNCONSTRAINED (bundle absent)")
    for t in res["tags"]:
        print(f"  [tag] {t}")


def _report(tag, res):
    print(f"\n== DG full uncertainty ({tag}), synthetic NE stand ==")
    for t in res["tree_pi"]:
        print(f"  tree spcd {t['spcd']:3d}: dg cm/yr  q05={t['q05']:.3f}  "
              f"q50={t['q50']:.3f}  q95={t['q95']:.3f}")
    sp = res["stand_pi"]
    print(f"  STAND mean dg: q05={sp['q05']:.3f} q50={sp['q50']:.3f} "
          f"q95={sp['q95']:.3f}  sd={sp['sd']:.4f}")
    print(f"  stand SD with shared-RE correlation: {sp['sd']:.4f}")
    print(f"  stand SD naive independent (WRONG):  {res['stand_naive_sd']:.4f}")
    infl = sp['sd'] / res['stand_naive_sd'] if res['stand_naive_sd'] > 0 else float('nan')
    print(f"  correlation inflation factor: {infl:.2f}x")
    for k in ("resid_source", "coef_posterior", "note"):
        if res.get(k):
            print(f"  [{k}] {res[k]}")


if __name__ == "__main__":
    CFG = os.environ.get("FVS_CONFIG_DIR", "/fs/scratch/PUOM0008/crsfaaron/wt-gompit/config")
    eco = {"L1": "8", "L2": "8.1", "L3": "8.1.1", "FT": 101}
    drivers = {"bgi": 6.0}
    ND = 2000

    # --- Leg B (conus_sf): kuehne_v8 dg covariates per tree ------------------
    Lb = FvsConfigLoader("ne", version="conus_sf", config_dir=CFG)
    def tree_b(spcd, dbh_cm, cr, bal_sw, bal_hw, bgi=6.0):
        return {"spcd": spcd, "b1": math.log(dbh_cm), "b2": dbh_cm,
                "b3": math.log((cr + 0.2) / 1.2), "b4": math.log(bal_sw + 0.01),
                "b5": bal_hw, "b6": bgi}
    stand_b = [tree_b(12, 15, .5, 5, 2), tree_b(97, 20, .45, 8, 3),
               tree_b(12, 25, .4, 12, 4), tree_b(97, 30, .5, 15, 5)]
    _report("conus_sf / Leg B", dg_uncertainty_leg_b(Lb, stand_b, eco, drivers=drivers, n_draws=ND))

    # --- Leg A (conus): per-species intercept posterior over a fixed base ----
    La = FvsConfigLoader("ne", version="conus", config_dir=CFG)
    def tree_a(spcd, dbh_cm, cr, bal_sw, bal_hw):
        # Leg A kuehne base covariates keyed by fixed param names b1..b8
        return {"spcd": spcd,
                "b1": math.log(dbh_cm), "b2": dbh_cm,
                "b3": math.log((cr + 0.2) / 1.2),
                "b4": bal_sw + bal_hw, "b5": dbh_cm ** 2,
                "b6": math.log(bal_sw + 0.01), "b7": bal_hw, "b8": cr}
    stand_a = [tree_a(12, 15, .5, 5, 2), tree_a(97, 20, .45, 8, 3),
               tree_a(12, 25, .4, 12, 4), tree_a(97, 30, .5, 15, 5)]
    _report("conus / Leg A", dg_uncertainty_leg_a(La, stand_a, eco, drivers=drivers, n_draws=ND))

    # --- ORGANON (conus_organon): native a1..a6 covariates -------------------
    Lo = FvsConfigLoader("ne", version="conus_organon", config_dir=CFG)
    def tree_o(spcd, dbh_cm, cr, sba, bal):
        return {"spcd": spcd,
                "a1": math.log(dbh_cm), "a2": dbh_cm ** 2,
                "a3": math.log((cr + 0.2) / 1.2), "a4": math.log(max(sba, 0.1)),
                "a5": bal / math.log(dbh_cm + 1), "a6": math.sqrt(max(sba, 0.0))}
    stand_o = [tree_o(202, 15, .5, 20, 7), tree_o(17, 20, .45, 25, 11),
               tree_o(202, 25, .4, 30, 15), tree_o(17, 30, .5, 35, 19)]
    _report("conus_organon", dg_uncertainty_organon(Lo, stand_o, eco, drivers=drivers, n_draws=ND))

    # --- Stand-constrained MORTALITY path (TASK 1) ---------------------------
    # Synthetic NE stand: per-tree base ANNUAL hazard from a plausible arm
    # survival model (smaller / more suppressed trees -> higher hazard). This
    # is the arm's own per-tree hazard the stand model reconciles.
    def tree_m(spcd, dbh_cm, cr, bal, tpa):
        # base hazard: increases with competition (bal) and decreases with size
        # and crown ratio. eta form mirrors the tree survival scale (h=exp(-eta)).
        eta = 2.6 + 0.9 * math.log(dbh_cm) + 1.2 * cr - 0.05 * bal
        return {"spcd": spcd, "dbh_cm": dbh_cm, "tpa": tpa,
                "base_hazard": math.exp(-eta)}
    # a moderately dense stand so the density cap is exercised
    stand_m = [tree_m(12, 12, .35, 14, 180), tree_m(97, 18, .40, 10, 120),
               tree_m(12, 24, .45, 7, 70),  tree_m(97, 34, .55, 4, 30)]
    SAMPLES = "/users/PUOM0008/crsfaaron/fvs-conus/output/variants/ne/stand_density_samples.rds"
    REAL_BUNDLE = ("/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/"
                   "stand_survival/stand_survival_bundle.json")
    STUB_BUNDLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "stand_survival_bundle_STUB.json")

    # (1) REAL bundle path: absent while the fit is pending -> graceful skip.
    print("\n--- stand-constrained mortality: REAL bundle path ---")
    res_real = stand_constrained_mortality_uncertainty(
        Lb, stand_m, eco, T_years=10.0, n_draws=ND, drivers={"bgi": 6.0},
        stand_survival_bundle_path=REAL_BUNDLE, sdimax_samples_path=SAMPLES)
    _report_mortality(res_real)

    # (2) STUB bundle path: demonstrates constrained-matches-target + density cap
    #     (STUB coefficients; replace with the real bundle when 71_* lands).
    print("\n--- stand-constrained mortality: STUB bundle path (demo of match) ---")
    res_stub = stand_constrained_mortality_uncertainty(
        Lb, stand_m, eco, T_years=10.0, n_draws=ND, drivers={"bgi": 6.0},
        stand_survival_bundle_path=STUB_BUNDLE, sdimax_samples_path=SAMPLES)
    _report_mortality(res_stub)

    print("\nUNCERTAINTY_OK")
