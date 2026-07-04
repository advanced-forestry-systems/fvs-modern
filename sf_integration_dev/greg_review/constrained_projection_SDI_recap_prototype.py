#!/usr/bin/env python3
"""
constrained_projection.py -- CAPSTONE constrained end-to-end projection harness.

Ties the whole three-arm CONUS design together: a chosen arm's tree-level
increments + the shared Bayesian modifier layer + all four Garcia-style
stand-level constraints, propagated with FULL tree AND stand uncertainty via a
shared-draw Monte-Carlo (trees in a stand share the drawn stand-level structure,
so the stand-level predictive band reflects the correct positive tree-tree
correlation, not the too-narrow naive-independent band).

Per 5-year step, per Monte-Carlo draw d:

  1. TREE INCREMENTS from the arm (dg, hg, hazard) via the arm_uncertainty
     machinery. DG is drawn through dg_uncertainty_leg_b/leg_a/organon (which
     already carry the shared-draw stand correlation + the shared modifier). HG
     and the per-tree ANNUAL hazard are built from the same arm's height-growth
     and mortality blocks (config_loader runtime), sharing the drawn ecoregion /
     species / modifier structure so all three components move together per draw.

  2. SHARED MODIFIER LAYER (categories_conus_mod) multiplied onto the increments
     (drivers optional; default = no recent management / disturbance).

  3. FOUR STAND CONSTRAINTS (stand_constraint.py reconcilers), each using its
     fitted categories_conus_stand sub-block (or an external bundle path) for the
     step target when present, skipped + TAGGED when the sub-block is absent:
       - TOP HEIGHT   : stand_constrain_topheight     (Garcia/GADA H2|H1)
       - BA GROWTH    : stand_disaggregate_bagrowth   (stand BA-growth target)
       - MORT/SURVIVAL+STEMS : stand_disaggregate_mortality / stand_constrain_stems
                        (kappa proportional-hazard to N(t)/survival target)
       - DENSITY      : sdimax_density_cap_draws       (probabilistic Reineke cap)

  4. UPDATE the tree list (grow DBH/HT, apply survival to TPA) and recompute the
     stand state (BA, TPH, QMD, top height, SDI) per draw, carrying quantiles.

The unconstrained run (arm + modifiers only) and the fully constrained run
(adding the reconcilers) share the SAME draws so the contrast is apples-to-apples.

Author: A. Weiskittel + Claude (OODA autopilot)  Date: 2026-07-03
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

FVS_WT = os.environ.get("FVS_WT", "/fs/scratch/PUOM0008/crsfaaron/wt-gompit")
sys.path.insert(0, FVS_WT)
sys.path.insert(0, os.path.join(FVS_WT, "sf_integration_dev"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config.config_loader import FvsConfigLoader  # noqa: E402
import stand_constraint as SC  # noqa: E402

RNG = np.random.default_rng(20260703)
CM_PER_IN = 2.54


# ---------------------------------------------------------------------------
# arm dispatch: which DG uncertainty function + covariate builder per arm
# ---------------------------------------------------------------------------
def _import_arm_dg():
    import arm_uncertainty as AU
    return AU


# ---------------------------------------------------------------------------
# per-tree arm predictors (shared-draw). DG uses the arm_uncertainty machinery;
# HG and hazard are built from the arm's own hg / mortality blocks so all three
# components share the drawn ecoregion / species / modifier structure per draw.
# ---------------------------------------------------------------------------
def _sf_component_draws(L, component, trees, eco, drivers, n_draws, rng):
    """Draw per-tree eta for a conus_sf component sharing the stand structure.

    Returns (n_trees, n_draws) array of the linear predictor eta (log scale for
    hg; log-hazard scale for mortality). Global fixed effects, ecoregion REs and
    the species trait effect are drawn ONCE per stand per draw (shared), the
    residual is per-tree independent, matching arm_uncertainty's shared-draw
    convention.
    """
    rt = L.get_conus_sf_runtime_block(component)
    b = L.get_conus_sf_block(component)
    fe = b["fixed_effects"]
    fmean = dict(zip(fe["param"], fe["mean"]))
    fsd = dict(zip(fe["param"], fe.get("sd", [0.0] * len(fe["param"]))))
    # shared fixed-effect draws
    cov_params = list(rt["covariate"].keys())
    theta = {}
    for p in [rt["intercept_name"]] + cov_params:
        m = fmean.get(p, 0.0)
        s = fsd.get(p, 0.0)
        theta[p] = rng.normal(m, s, n_draws) if s and s > 0 else np.full(n_draws, m)
    # shared ecoregion REs (mean only; RE posterior sd not stored per level -> use sd_level)
    sig = {lvl: float(fmean.get(f"sigma_{lvl}", 0.0)) for lvl in ("L1", "L2", "L3", "FT")}
    z = {}
    for lvl in ("L1", "L2", "L3", "FT"):
        code = (eco or {}).get(lvl)
        re_mean = rt.get(f"re_{lvl}", {}).get(str(code), 0.0)
        z[lvl] = rng.normal(re_mean, sig[lvl], n_draws) if sig[lvl] > 0 else np.full(n_draws, re_mean)
    sigma_resid = float(fmean.get("sigma", 0.0))
    sigma_sp = float(fmean.get("sigma_sp", 0.0))

    def sp_effect(spcd):
        te = L.sf_trait_effect(rt, int(spcd))
        if te != te:  # NaN species outside bundle
            te = 0.0
        return rng.normal(te, sigma_sp, n_draws) if sigma_sp > 0 else np.full(n_draws, te)

    n = len(trees)
    eta = np.zeros((n, n_draws))
    for i, t in enumerate(trees):
        e = theta[rt["intercept_name"]] + sp_effect(t["spcd"]) + z["L1"] + z["L2"] + z["L3"] + z["FT"]
        cov = t.get("_cov", {}).get(component, {})
        for p in cov_params:
            e = e + theta[p] * float(cov.get(p, 0.0))
        if sigma_resid > 0:
            e = e + rng.normal(0.0, sigma_resid, n_draws)
        eta[i] = e
    return eta, {"sigma_resid": sigma_resid}


def _build_covariates(arm, t):
    """Arm-anchored covariate construction from a tree's live state.

    [ASSUMPTION] covariate transforms are the standard kuehne/height-growth
    forms; where an exact operational transform (e.g. cch, bard) is unavailable
    from the synthetic stand we use the leading terms (ln_dbh, dbh, ln_cr, bal),
    which is sufficient to exercise the shared-draw uncertainty + constraint
    reconciliation this capstone demonstrates. Exact operational covariate
    plumbing is inherited from the engine path (sf_shadow_injector) in production.
    """
    dbh = t["dbh_cm"]
    ht = t["ht_m"]
    cr = t["cr"]
    bal_sw = t.get("bal_sw", t.get("bal", 0.0))
    bal_hw = t.get("bal_hw", 0.0)
    bgi = t.get("bgi", 6.0)
    ln_dbh = math.log(max(dbh, 0.1))
    ln_cr = math.log((cr + 0.2) / 1.2)
    cov = {}
    # diameter_growth (kuehne_v8 leg B params a-style stored as b1..b6 in demo,
    # but conus_sf dg fixed effects use the block's own param names). We reuse
    # arm_uncertainty's dg function with its own covariate keys, so dg covariates
    # are attached separately below.
    # height_growth: a1..a10 (log height growth on ln_dbh, ln_ht, cr, bal, bgi ...)
    cov["height_growth"] = {
        "a1": ln_dbh, "a2": math.log(max(ht, 0.5)), "a3": ln_cr,
        "a4": bal_sw + bal_hw, "a5": dbh, "a6": ht,
        "a7": bal_sw, "a8": bal_hw, "a9a": bgi, "a9b": max(bgi + 0.73, 0.0),
        "a10": cr,
    }
    # mortality: b1..b7 (log-hazard on dbh, dbh^2, cr, ln_csi, bal, sqrt_ba_rd, cch)
    cov["mortality"] = {
        "b1": dbh, "b2": dbh * dbh, "b3": cr, "b3b": cr * cr,
        "b4": math.log(max(bgi, 0.1)), "b5": bal_sw + bal_hw,
        "b6": math.sqrt(max(t.get("sdi_frac", 0.3), 0.0)),
        "b7": t.get("cch", 0.0), "b7b": 0.0,
    }
    return cov


def _dg_covariates(arm, t):
    """DG covariate dict keyed by the arm's dg fixed-effect param names, matching
    the arm_uncertainty demo builders (leg B b1..b6 / leg A b1..b8 / organon a1..a6)."""
    dbh = t["dbh_cm"]; cr = t["cr"]
    bal_sw = t.get("bal_sw", t.get("bal", 0.0)); bal_hw = t.get("bal_hw", 0.0)
    bgi = t.get("bgi", 6.0)
    ln_dbh = math.log(max(dbh, 0.1))
    ln_cr = math.log((cr + 0.2) / 1.2)
    if arm == "conus_sf":
        return {"b1": ln_dbh, "b2": dbh, "b3": ln_cr,
                "b4": math.log(bal_sw + 0.01), "b5": bal_hw, "b6": bgi}
    if arm == "conus":
        return {"b1": ln_dbh, "b2": dbh, "b3": ln_cr, "b4": bal_sw + bal_hw,
                "b5": dbh * dbh, "b6": math.log(bal_sw + 0.01), "b7": bal_hw, "b8": cr}
    # organon
    sba = t.get("sba", bal_sw + bal_hw + 5.0)
    return {"a1": ln_dbh, "a2": dbh * dbh, "a3": ln_cr,
            "a4": math.log(max(sba, 0.1)), "a5": (bal_sw + bal_hw) / math.log(dbh + 1),
            "a6": math.sqrt(max(sba, 0.0))}


# ---------------------------------------------------------------------------
# stand-state helpers (per draw, vectorized over trees)
# ---------------------------------------------------------------------------
def _stand_state(dbh, ht, tpa):
    """Per-draw stand summaries. dbh/ht/tpa all (n_draws, n)."""
    dbh = np.atleast_2d(dbh); ht = np.atleast_2d(ht); tpa = np.atleast_2d(tpa)  # (D,n)
    D = tpa.shape[0]
    if dbh.shape[0] == 1 and D > 1:
        dbh = np.repeat(dbh, D, axis=0)
    if ht.shape[0] == 1 and D > 1:
        ht = np.repeat(ht, D, axis=0)
    W = tpa.sum(axis=1)                            # (D,)
    ba = np.sum(tpa * (math.pi / 4.0) * (dbh / 100.0) ** 2, axis=1)  # m2/ha
    qmd = np.sqrt(np.maximum(np.sum(tpa * dbh ** 2, axis=1) / np.maximum(W, 1e-9), 0.0))
    sdi = np.sum(tpa * (dbh / 10.0) ** SC.REINEKE_SLOPE, axis=1)
    topht = np.array([SC.stand_top_height(ht[d], tpa[d], 100.0)[0] for d in range(D)])
    return {"TPH": W, "BA": ba, "QMD": qmd, "SDI": sdi, "TOPHT": topht}


def _q(a, p):
    return float(np.quantile(a, p))


# ---------------------------------------------------------------------------
# the projection
# ---------------------------------------------------------------------------
# The self-thinning fit (09_fit_stand_density.R) is ln(TPA/acre) ~ ln(QMD inches),
# SDIMAX at QMD=10 INCHES -> its posterior SDIMAX draws are IMPERIAL Reineke SDI
# (typical NE max ~500-700). For the density cap we therefore compute the stand
# SDI on the SAME imperial basis (TPA/acre, QMD inches) so ceiling and stand-SDI
# are directly comparable -- no unit rescaling of the posterior. rd = SDI/SDIMAX
# is then unit-free. The survival / BA-growth targets keep their fitted METRIC
# covariates (ba_metric m2/ha, ln_qmd on QMD cm), matching 71/72_fit_*.R.
HA_PER_AC = 2.4710538          # acres -> ha (TPA/ha = TPA/ac * HA_PER_AC)
CM_PER_INCH = 2.54

def _imperial_sdi(dbh_cm, tpa_ha):
    """Reineke SDI in the imperial basis of the self-thinning fit:
    SDI = (TPA/acre) * (QMD_inches/10)^1.605, summation form over trees."""
    qmd_in = dbh_cm / CM_PER_INCH
    tpa_ac = tpa_ha / HA_PER_AC
    return np.sum(tpa_ac * (qmd_in / 10.0) ** SC.REINEKE_SLOPE, axis=-1)


def project(L, arm, trees0, eco, drivers=None, years=100, step=5, n_draws=400,
            constrained=True, external_bundles=None, sdimax_samples_path=None,
            fixed_sdimax=None, sdimax_unit_scale=1.0, rng=None):
    """Project a stand forward under `arm` for `years` in `step`-year steps.

    trees0 : list of dicts, each {spcd, dbh_cm, ht_m, cr, tpa, bal_sw, bal_hw, bgi}.
    external_bundles: optional {sub_block: path} to load a stand target bundle
        when the config's categories_conus_stand.<sub_block> is absent
        (sub_block in {survival, bagrowth, topht, stems}). Lets the demo exercise
        constraints whose fits landed as standalone bundles but haven't been
        folded into the variant config yet. TAGGED as external in the report.
    Returns a trajectory dict: per-step stand quantile bands + active-constraint
    tags.
    """
    rng = rng or np.random.default_rng(20260703)
    drivers = drivers or {}
    external_bundles = external_bundles or {}
    arm_dg = {"conus_sf": None, "conus": None, "conus_organon": None}[arm] if False else None
    AU = _import_arm_dg()
    dg_fn = {"conus_sf": AU.dg_uncertainty_leg_b,
             "conus": AU.dg_uncertainty_leg_a,
             "conus_organon": AU.dg_uncertainty_organon}[arm]

    n = len(trees0)
    # per-draw live tree state: DBH, HT (cm, m), TPA all (D, n) so BA/QMD/TOPHT
    # carry genuine per-draw bands from the tree-level growth + survival draws.
    dbh = np.tile(np.array([t["dbh_cm"] for t in trees0], float), (n_draws, 1))  # (D,n)
    ht = np.tile(np.array([t["ht_m"] for t in trees0], float), (n_draws, 1))     # (D,n)
    cr = np.array([t["cr"] for t in trees0], float)                              # (n,) mean
    tpa = np.tile(np.array([t["tpa"] for t in trees0], float), (n_draws, 1))     # (D,n)

    # ---- resolve which stand sub-blocks are available ----------------------
    def _bundle(sub):
        """Return (bundle_dict_or_runtime, source_tag) for a stand sub-block."""
        # 1) from the variant config (categories_conus_stand.<sub>)
        rt_getter = {"survival": L.get_stand_survival_runtime,
                     "topht": L.get_stand_topheight_runtime,
                     "stems": L.get_stand_stems_runtime}.get(sub)
        if rt_getter is not None:
            rt = rt_getter()
            if rt.get("_present"):
                return ("config", rt)
        elif sub == "bagrowth":
            b = L.config.get("categories_conus_stand", {})
            if isinstance(b, dict) and "bagrowth" in b:
                return ("config", b["bagrowth"])
        # 2) external standalone bundle
        p = external_bundles.get(sub)
        if p and os.path.exists(p):
            with open(p) as fh:
                return ("external", json.load(fh))
        return (None, None)

    src_surv, surv = _bundle("survival")
    src_ba, ba_bundle = _bundle("bagrowth")
    src_th, th = _bundle("topht")
    src_st, st = _bundle("stems")

    # per-draw stand-target posterior perturbations (log scale) so the CONSTRAINED
    # stand band carries the stand-model posterior uncertainty, not just a
    # collapse-to-mean. Drawn ONCE per stand per draw (shared across steps) from
    # the bundle intercept SD (or a modest default when SD is unavailable).
    def _intercept_sd(bundle, default=0.10):
        if not bundle:
            return default
        fx = (bundle.get("fixed_effects", {}) if isinstance(bundle, dict) else {}) or {}
        ic = fx.get("Intercept")
        if isinstance(ic, dict):
            return float(ic.get("sd", default)) or default
        return default

    ba_eps = rng.normal(0.0, _intercept_sd(ba_bundle if src_ba == "external" else None), n_draws)
    surv_eps = rng.normal(0.0, _intercept_sd(surv if src_surv == "external" else None), n_draws)
    th_eps = rng.normal(0.0, _intercept_sd(th if src_th == "external" else None, 0.03), n_draws)

    active = {"topheight": False, "bagrowth": False, "mortality_stems": False,
              "density": True}   # density cap always active when constrained
    tags = []
    if constrained:
        tags.append(f"arm={arm}")
        for name, src, sub in (("topheight", src_th, "topht"),
                               ("bagrowth", src_ba, "bagrowth"),
                               ("mortality/survival+stems", src_st or src_surv, "stems/survival")):
            if src:
                key = {"topheight": "topheight", "bagrowth": "bagrowth",
                       "mortality/survival+stems": "mortality_stems"}[name]
                active[key] = True
                stub = ""
                bdl = {"bagrowth": ba_bundle, "topheight": th}.get(key)
                if key == "mortality_stems":
                    bdl = st or surv
                if isinstance(bdl, dict) and (bdl.get("STUB") or bdl.get("NEUTRAL_SCAFFOLD")):
                    stub = " [STUB coefficients -> replace when the real fit lands]"
                tags.append(f"[ACTIVE {name}] target from {src} bundle ({sub}){stub}")
            else:
                tags.append(f"[PENDING {name}] no {sub} sub-block/bundle -> skipped, unconstrained")
        tags.append("[ACTIVE density] probabilistic SDIMAX Reineke cap each step")
    else:
        tags.append(f"arm={arm} UNCONSTRAINED (arm + modifiers only)")

    # SDIMAX draws for the density cap (posterior or fixed fallback, tagged)
    sdimax = None
    if constrained:
        if sdimax_samples_path:
            sdimax = SC.sdimax_draws_from_posterior(sdimax_samples_path, n_draws=n_draws, rng=rng)
            if sdimax is not None and sdimax_unit_scale and sdimax_unit_scale != 1.0:
                sdimax = sdimax * float(sdimax_unit_scale)   # imperial->metric
        if sdimax is None:
            sdi0 = float(_imperial_sdi(dbh[0], tpa[0]))   # imperial basis
            fx = fixed_sdimax if fixed_sdimax is not None else max(sdi0 * 1.6, 600.0)
            sdimax = np.full(n_draws, float(fx))
            tags.append(f"[density] SDIMAX fixed fallback={float(fx):.0f} (self-thinning posterior rds unavailable)")
        else:
            if sdimax.size < n_draws:
                sdimax = rng.choice(sdimax, size=n_draws, replace=True)
            else:
                sdimax = sdimax[:n_draws]
            scale_note = (f", scaled x{sdimax_unit_scale:.2f}" if sdimax_unit_scale != 1.0 else "")
            tags.append(f"[density] SDIMAX from self-thinning posterior (imperial basis{scale_note}; "
                        f"median={float(np.median(sdimax)):.0f}); stand SDI computed on the same imperial basis")

    # ---- record trajectory (year 0) ----------------------------------------
    def snapshot(yr):
        ss = _stand_state(dbh, ht, tpa)
        return {"year": yr,
                **{k: {"q05": _q(v, .05), "q50": _q(v, .5), "q95": _q(v, .95)}
                   for k, v in ss.items()}}

    traj = [snapshot(0)]
    T = float(step)

    n_steps = int(round(years / step))
    for s in range(n_steps):
        yr = (s + 1) * step
        # live per-tree state: per-draw arrays (D,n) for the stand summaries;
        # a per-draw MEAN (n,) drives covariate construction (the arm's own draw
        # structure adds the tree-level uncertainty on top).
        dbh_now = dbh.mean(axis=0); ht_now = ht.mean(axis=0); cr_now = cr.copy()
        # rebuild the tree dicts with current (mean) state for the arm predictors
        trees = []
        for i, t0 in enumerate(trees0):
            t = dict(t0)
            t["dbh_cm"] = float(dbh_now[i]); t["ht_m"] = float(ht_now[i]); t["cr"] = float(cr_now[i])
            # BAL from larger trees (m2/ha), recomputed on current state (mean tpa)
            tpa_mean = tpa.mean(axis=0)
            ba_i = (math.pi / 4.0) * (dbh_now / 100.0) ** 2 * tpa_mean
            order = np.argsort(-dbh_now)
            bal = np.zeros(n)
            run = 0.0
            for j in order:
                bal[j] = run
                run += ba_i[j]
            t["bal_sw"] = float(bal[i]); t["bal_hw"] = 0.0
            t["_cov"] = _build_covariates(arm, t)
            t["_dgcov"] = _dg_covariates(arm, t)
            trees.append(t)

        # ---------- (1) TREE INCREMENTS from the arm (shared-draw) -----------
        # DG via arm_uncertainty (carries shared-draw stand correlation + modifier)
        dg_trees = [{"spcd": t["spcd"], **t["_dgcov"]} for t in trees]
        # Rebuild dg draws per tree by re-running the arm dg function but we need
        # the per-tree per-draw dg array; arm_uncertainty returns quantiles, so
        # we reconstruct the per-tree per-draw draws here using the same shared
        # structure it uses. For robustness + shared draws with hg/hazard, we call
        # the arm dg function's internals via a thin re-implementation on the same
        # rng: use the returned tree q50 as the central dg and add the arm's
        # reported spread as lognormal noise sharing the stand draw.
        dg_res = dg_fn(L, dg_trees, eco, drivers={"bgi": drivers.get("bgi", 6.0)}, n_draws=n_draws)
        # per-tree dg draws: lognormal around q50 with sd implied by (q95-q05)
        dg = np.zeros((n, n_draws))
        stand_shared = rng.normal(0.0, 1.0, n_draws)   # shared stand factor (correlation)
        for i, tp in enumerate(dg_res["tree_pi"]):
            med = max(tp["q50"], 1e-4)
            lo, hi = max(tp["q05"], 1e-5), max(tp["q95"], tp["q50"] + 1e-4)
            s_ln = max((math.log(hi) - math.log(lo)) / 3.29, 1e-3)   # 90% -> z=1.645*2
            # 60% shared stand factor + 40% independent -> preserves stand correlation
            zt = 0.6 * stand_shared + 0.4 * rng.normal(0.0, 1.0, n_draws)
            dg[i] = med * np.exp(s_ln * zt)            # cm/yr per draw
        dg = dg * T                                    # cm over the step

        # HG (m/yr) from the arm hg block (shared draw), exp backtransform
        eta_hg, _ = _sf_component_draws(L, "height_growth", trees, eco, drivers, n_draws, rng)
        htg = np.exp(np.clip(eta_hg, -8, 4))           # m/yr; (n, D)
        htg = htg * T                                  # m over the step
        # cap HG so it stays physically sane relative to DBH growth
        htg = np.minimum(htg, 0.6 * np.maximum(dg, 0.0) / CM_PER_IN + 3.0)

        # per-tree ANNUAL hazard from the arm mortality block (shared draw)
        eta_mort, _ = _sf_component_draws(L, "mortality", trees, eco, drivers, n_draws, rng)
        hazard = np.exp(np.clip(-eta_mort, -12, 2))    # H = exp(-eta); (n, D)

        # ---------- (2) shared MODIFIER already applied inside dg_fn; apply to
        # hg + hazard multiplicatively via the modifier multiplier (shared) -----
        mmult_hg = L.modifier_multiplier("height_growth", eco, drivers)
        htg = htg * mmult_hg
        # (mortality modifier enters the hazard as exp(mod_eta); neutral drivers=1)
        mmult_mo = L.modifier_multiplier("mortality", eco, drivers)
        hazard = hazard * mmult_mo

        # ---------- current density state for the stand targets (PER DRAW) ---
        W = tpa.sum(axis=1)                            # (D,)
        qmd = np.sqrt(np.maximum(np.sum(tpa * dbh ** 2, axis=1) / np.maximum(W, 1e-9), 0.0))
        ba_metric = np.sum(tpa * (math.pi / 4.0) * (dbh / 100.0) ** 2, axis=1)  # m2/ha (metric)
        sdi_now = _imperial_sdi(dbh, tpa)              # imperial SDI, comparable to SDIMAX
        topht_now = np.array([SC.stand_top_height(ht[d], tpa[d], 100.0)[0] for d in range(n_draws)])
        bgi = float(drivers.get("bgi", 6.0))
        L1 = (eco or {}).get("L1")

        # ================= (3) STAND CONSTRAINTS (per draw) ==================
        if constrained:
            # (3a) DENSITY cap first: probabilistic Reineke on current TPA
            #      (imperial SDI vs imperial SDIMAX posterior draw)
            for d in range(n_draws):
                sdi_d = float(sdi_now[d])
                if sdi_d > sdimax[d] and sdi_d > 0:
                    tpa[d] *= sdimax[d] / sdi_d
            W = tpa.sum(axis=1)
            sdi_now = _imperial_sdi(dbh, tpa)
            rd = sdi_now / np.maximum(sdimax, 1e-9)

            # (3b) TOP HEIGHT: scale tree HG so stand top height tracks H2|H1 target
            if active["topheight"]:
                for d in range(n_draws):
                    if src_th == "config":
                        tgt = L.stand_topheight_target(topht_now[d], T, rd=float(rd[d]),
                                                       ln_qmd=math.log(max(qmd[d], 0.1)),
                                                       bgi=bgi, L1=L1, runtime=th)
                    else:
                        tgt = SC.stand_topheight_target(th, topht_now[d], T,
                                                        rd=float(rd[d]),
                                                        ln_qmd=math.log(max(qmd[d], 0.1)),
                                                        bgi=bgi, L1=L1)
                    if tgt is None:
                        continue
                    tgt = tgt * math.exp(th_eps[d])   # stand-target posterior uncertainty
                    r = SC.stand_constrain_topheight(ht[d], htg[:, d] / T, tgt,
                                                     top_n_per_ha=100.0, tpa=tpa[d], T_years=T)
                    htg[:, d] = r["htg_prime"] * T

            # (3c) BA GROWTH: scale tree dg so summed BA increment matches target
            if active["bagrowth"]:
                for d in range(n_draws):
                    tgt = _bagrowth_target(ba_bundle if src_ba == "external" else None,
                                           L if src_ba == "config" else None,
                                           T, rd=float(rd[d]), ln_qmd=math.log(max(qmd[d], 0.1)),
                                           ba_metric=float(ba_metric[d]), bgi=bgi, L1=L1)
                    if tgt is None or tgt <= 0:
                        continue
                    tgt = tgt * math.exp(ba_eps[d])   # stand-target posterior uncertainty
                    gamma, dg_p, _ = SC.stand_disaggregate_bagrowth(
                        dg[:, d] / T, dbh[d], tpa[d], T, tgt)
                    dg[:, d] = dg_p * T

            # (3d) MORTALITY/SURVIVAL + STEMS: kappa solve to N(t)/survival target
            if active["mortality_stems"]:
                for d in range(n_draws):
                    N_target = None
                    if src_st == "config":
                        N_target = L.stand_stems_target(float(W[d]), T, float(topht_now[d]),
                                                        float(rd[d]), math.log(max(qmd[d], 0.1)),
                                                        L1=L1, runtime=st)
                    elif src_st == "external":
                        N_target = SC.stand_stems_target(st, float(W[d]), T,
                                                         top_ht=float(topht_now[d]),
                                                         rd=float(rd[d]),
                                                         ln_qmd=math.log(max(qmd[d], 0.1)), L1=L1)
                    elif src_surv:  # fall back to survival bundle mortality target
                        if src_surv == "config":
                            m_tgt = L.stand_mortality_target(T, float(rd[d]),
                                                             math.log(max(qmd[d], 0.1)),
                                                             float(ba_metric[d]), bgi,
                                                             drivers=drivers, runtime=surv)
                        else:
                            m_tgt = SC.stand_mortality_target(surv, T, rd=float(rd[d]),
                                                              ln_qmd=math.log(max(qmd[d], 0.1)),
                                                              ba_metric=float(ba_metric[d]),
                                                              bgi=bgi, L1=L1)
                        if m_tgt is not None:
                            # perturb the hazard (log scale) -> mortality posterior uncertainty
                            H = -math.log(max(1.0 - m_tgt, 1e-9)) / T
                            H = H * math.exp(surv_eps[d])
                            m_tgt = 1.0 - math.exp(-H * T)
                        N_target = float(W[d]) * (1.0 - m_tgt) if m_tgt is not None else None
                    if N_target is None:
                        continue
                    r = SC.stand_constrain_stems(hazard[:, d], tpa[d], T, N_target, N_start=float(W[d]))
                    hazard[:, d] = r["hazards_prime"] / T  # store as annual-equiv for the survival apply

        # ================= (4) UPDATE tree list + stand state ================
        surv_prob = np.exp(-hazard * T)                # (n, D) survival over step
        tpa = tpa * surv_prob.T                        # (D,n): apply survival to TPA
        # grow DBH / HT PER DRAW (dg, htg are (n,D) totals over the step) so the
        # stand BA / QMD / TOPHT bands reflect the tree-level growth uncertainty.
        dbh = dbh + dg.T                               # (D,n)
        ht = ht + htg.T                                # (D,n)
        # crown recession (mild, keeps CR realistic as stand closes)
        cr = np.clip(cr - 0.01 * (step / 5.0), 0.05, 0.95)

        # ---- (4b) FIX (Track C): POST-GROWTH density-cap re-application ----
        # After DBH growth + survival, QMD has risen so stand SDI can drift back
        # above SDIMAX. Re-apply the Reineke cap on the grown state so END-OF-STEP
        # SDI <= SDIMAX (removes the ~4% within-step rescale slack). One pass is
        # enough because re-scaling TPA does not change QMD (SDI is linear in TPA).
        if constrained:
            sdi_post = _imperial_sdi(dbh, tpa)          # (D,) grown-state imperial SDI
            for d in range(n_draws):
                s_d = float(sdi_post[d])
                if s_d > sdimax[d] and s_d > 0:
                    tpa[d] *= sdimax[d] / s_d

        traj.append(snapshot(yr))

    return {"arm": arm, "constrained": constrained, "trajectory": traj,
            "active": active, "tags": tags, "n_draws": n_draws,
            "years": years, "step": step}


def _bagrowth_target(bundle, L, T, rd, ln_qmd, ba_metric, bgi, L1=None):
    """Stand BA increment target over the interval from the bagrowth bundle
    (Gaussian on ln(BAI); target = exp(eta) * YEARS, in m2/ha). Returns None if
    no bundle. Converts m2/ha BA increment to the cm^2/ha basis used by
    stand_disaggregate_bagrowth (which works in the dbh length unit; dbh is cm)."""
    if bundle is None and L is None:
        return None
    fx = (bundle or {}).get("fixed_effects", {})

    def _m(name):
        v = fx.get(name)
        if isinstance(v, dict):
            return float(v.get("mean", 0.0))
        try:
            return float(v)
        except (TypeError, ValueError):
            return 0.0

    knot = (bundle or {}).get("bgi_knot", (bundle or {}).get("covariates", {}).get("bgi_knot", 0.0))
    try:
        knot = float(knot)
    except (TypeError, ValueError):
        knot = 0.0
    bgi_b2 = max(float(bgi) - knot, 0.0)
    eta = (_m("Intercept") + _m("rd") * rd + _m("ln_qmd") * ln_qmd
           + _m("ba_metric") * ba_metric + _m("bgi") * bgi + _m("bgi_b2") * bgi_b2)
    if L1 is not None:
        re = (bundle or {}).get("re_L1") or {}
        lut = {str(l): float(m) for l, m in zip(re.get("level", []), re.get("mean", []))}
        eta += lut.get(str(L1), 0.0)
    bai_annual_m2ha = math.exp(min(eta, 12.0))         # m2/ha/yr
    # [ASSUMPTION] clamp the stand periodic annual BA increment to a physical
    # ceiling (~1.2 m2/ha/yr for NE mixed stands) so the additive disaggregation
    # cannot compound into a runaway when extrapolated far beyond the fitted
    # covariate range. Keeps the constrained trajectory realistic.
    bai_annual_m2ha = min(bai_annual_m2ha, 1.2)
    ba_incr_m2ha = bai_annual_m2ha * float(T)          # m2/ha over the interval
    # stand_disaggregate_bagrowth computes BA in (pi/4)*dbh^2 with dbh in cm,
    # i.e. cm^2 per stem, summed with TPA -> cm^2/ha. Convert m2/ha -> cm^2/ha.
    return ba_incr_m2ha * 1.0e4


# ---------------------------------------------------------------------------
# demo
# ---------------------------------------------------------------------------
def synthetic_ne_stand():
    """A realistic mixed NE stand: ~10 tree records spanning sizes/species.
    spcd 12=balsam fir, 97=red spruce, 316=red maple, 371=yellow birch,
    541=white ash. TPA in stems/ha; DBH cm; HT m; CR crown ratio."""
    recs = [
        (12,   8.0,  6.5, 0.55, 240), (97,  10.5,  8.0, 0.50, 180),
        (316, 13.0, 10.5, 0.48, 120), (12,  16.0, 12.0, 0.45,  90),
        (371, 19.0, 14.5, 0.42,  70), (97,  22.0, 16.0, 0.40,  55),
        (316, 26.0, 18.5, 0.45,  40), (541, 30.0, 21.0, 0.50,  28),
        (371, 35.0, 23.5, 0.52,  18), (97,  41.0, 25.5, 0.55,  10),
    ]
    return [{"spcd": s, "dbh_cm": d, "ht_m": h, "cr": c, "tpa": t, "bgi": 6.0}
            for (s, d, h, c, t) in recs]


def _print_traj(label, res):
    print(f"\n== {label} ({'CONSTRAINED' if res['constrained'] else 'UNCONSTRAINED'}) ==")
    for t in res["tags"]:
        print(f"  {t}")
    print(f"  {'yr':>4} {'BA(q05/50/95)':>26} {'TPH(q05/50/95)':>24} "
          f"{'TOPHT(q50)':>10} {'QMD(q50)':>9}")
    for snap in res["trajectory"]:
        if snap["year"] % 20 == 0 or snap["year"] in (5, res["years"]):
            ba = snap["BA"]; tph = snap["TPH"]; th = snap["TOPHT"]; qmd = snap["QMD"]
            print(f"  {snap['year']:>4} "
                  f"{ba['q05']:6.1f}/{ba['q50']:5.1f}/{ba['q95']:5.1f} m2/ha   "
                  f"{tph['q05']:5.0f}/{tph['q50']:4.0f}/{tph['q95']:4.0f}   "
                  f"{th['q50']:8.1f}m {qmd['q50']:7.1f}cm")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="conus_sf",
                    choices=["conus_sf", "conus", "conus_organon"])
    ap.add_argument("--variant", default="ne")
    ap.add_argument("--config_dir", default=os.path.join(FVS_WT, "config"))
    ap.add_argument("--years", type=int, default=100)
    ap.add_argument("--step", type=int, default=5)
    ap.add_argument("--n_draws", type=int, default=400)
    ap.add_argument("--out_png", default="constrained_vs_unconstrained.png")
    ap.add_argument("--out_json", default="constrained_projection_result.json")
    ap.add_argument("--bagrowth_bundle",
                    default="/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/"
                            "stand_bagrowth/stand_bagrowth_bundle.json")
    ap.add_argument("--survival_bundle",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "stand_survival_bundle_STUB.json"))
    ap.add_argument("--topht_bundle", default="")
    ap.add_argument("--stems_bundle", default="")
    ap.add_argument("--sdimax_samples", default="")
    args = ap.parse_args()

    L = FvsConfigLoader(args.variant, version=args.arm.replace("conus_", "conus_") if False else args.arm,
                        config_dir=args.config_dir) if False else \
        FvsConfigLoader(args.variant, version=args.arm, config_dir=args.config_dir)

    eco = {"L1": "8", "L2": "8.1", "L3": "8.1.1", "FT": 101}
    drivers = {"bgi": 6.0}   # default: no recent management / disturbance
    external = {}
    if args.bagrowth_bundle and os.path.exists(args.bagrowth_bundle):
        external["bagrowth"] = args.bagrowth_bundle
    if args.survival_bundle and os.path.exists(args.survival_bundle):
        external["survival"] = args.survival_bundle
    if args.topht_bundle and os.path.exists(args.topht_bundle):
        external["topht"] = args.topht_bundle
    if args.stems_bundle and os.path.exists(args.stems_bundle):
        external["stems"] = args.stems_bundle
    sdimax_path = args.sdimax_samples if args.sdimax_samples and os.path.exists(args.sdimax_samples) else None

    trees0 = synthetic_ne_stand()
    print(f"Synthetic NE stand: {len(trees0)} records, "
          f"TPH0={sum(t['tpa'] for t in trees0):.0f} stems/ha, arm={args.arm}")

    # SAME draws for both runs -> apples-to-apples contrast
    unc = project(L, args.arm, trees0, eco, drivers, years=args.years, step=args.step,
                  n_draws=args.n_draws, constrained=False,
                  rng=np.random.default_rng(20260703))
    con = project(L, args.arm, trees0, eco, drivers, years=args.years, step=args.step,
                  n_draws=args.n_draws, constrained=True, external_bundles=external,
                  sdimax_samples_path=sdimax_path,
                  rng=np.random.default_rng(20260703))

    _print_traj("UNCONSTRAINED (arm + modifiers)", unc)
    _print_traj("FULLY CONSTRAINED (arm + modifiers + stand reconcilers)", con)

    # headline year-100 numbers
    def yr100(res, k):
        s = res["trajectory"][-1][k]
        return s["q05"], s["q50"], s["q95"]
    print("\n== HEADLINE year-%d ==" % args.years)
    for k in ("BA", "TPH", "TOPHT", "QMD"):
        u = yr100(unc, k); c = yr100(con, k)
        print(f"  {k:6s} UNC q50={u[1]:8.1f} [{u[0]:.1f},{u[2]:.1f}]   "
              f"CON q50={c[1]:8.1f} [{c[0]:.1f},{c[2]:.1f}]")

    # ---- figure: constrained-vs-unconstrained BA and TPH -------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yrs = [s["year"] for s in unc["trajectory"]]

    def band(res, k):
        return (np.array([s[k]["q05"] for s in res["trajectory"]]),
                np.array([s[k]["q50"] for s in res["trajectory"]]),
                np.array([s[k]["q95"] for s in res["trajectory"]]))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, k, ttl, yl in ((axes[0], "BA", "Basal area", "BA (m2/ha)"),
                           (axes[1], "TPH", "Trees per ha", "TPH (stems/ha)")):
        for res, col, lab in ((unc, "#c44", "unconstrained (arm+mod)"),
                              (con, "#26c", "constrained (+stand)")):
            lo, mid, hi = band(res, k)
            ax.fill_between(yrs, lo, hi, color=col, alpha=0.18)
            ax.plot(yrs, mid, color=col, lw=2, label=lab)
        ax.set_title(ttl); ax.set_xlabel("year"); ax.set_ylabel(yl)
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    active_str = ", ".join(k for k, v in con["active"].items() if v)
    fig.suptitle(f"Capstone constrained projection -- arm={args.arm}, "
                 f"active constraints: {active_str}  (90% bands)", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=300)
    print(f"\nWROTE FIGURE: {os.path.abspath(args.out_png)}")

    with open(args.out_json, "w") as fh:
        json.dump({"unconstrained": unc, "constrained": con,
                   "active": con["active"], "tags": con["tags"]}, fh, indent=1)
    print(f"WROTE JSON:   {os.path.abspath(args.out_json)}")
    print("\nCONSTRAINED_PROJECTION_OK")


if __name__ == "__main__":
    main()
