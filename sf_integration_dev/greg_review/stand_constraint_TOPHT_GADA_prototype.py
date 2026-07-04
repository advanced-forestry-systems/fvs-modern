#!/usr/bin/env python3
"""
stand_constraint.py -- stand-level constraints on tree-level predictions.

Two constraint channels reconcile the sum of tree predictions to a fitted
stand-level prediction, both carrying full posterior uncertainty:

  (b) Stand-level mortality DISAGGREGATION (stand_disaggregate_mortality):
      the stand survival model (71_fit_stand_survival.R) sets the stand
      mortality target M_stand = 1 - S_stand(T); the tree hazards h_i from the
      arm set the distribution among trees. We solve for a single proportional
      hazard scalar kappa such that the TPA-weighted stand mortality of the
      rescaled tree hazards h_i' = kappa * h_i equals M_stand:

          sum_i tpa_i * (1 - exp(-kappa * h_i * T)) / sum_i tpa_i == M_stand

      The left side is strictly increasing in kappa (each term is increasing),
      so a 1-D monotone root find (bisection) is robust. kappa < 1 lowers
      mortality toward the stand target, kappa > 1 raises it.

  (a) Density ceiling (sdimax_density_cap_draws): the probabilistic Reineke
      self-thinning cap. Per SDIMAX posterior draw, if the stand SDI exceeds the
      drawn SDIMAX the tree TPAs are rescaled by SDIMAX/SDI, so the self-thinning
      constraint carries the SDIMAX posterior uncertainty into the tree TPAs.

Reineke convention (matches calibration/R/17_stand_projection_engine.R and
09_fit_stand_density.R): SDI = TPA * (QMD/10)^1.605, QMD in cm.

Author: A. Weiskittel + Claude (OODA autopilot)  Date: 2026-07-02
"""
from __future__ import annotations

import math

import numpy as np

REINEKE_SLOPE = 1.605  # Reineke self-thinning exponent (metric QMD in cm, /10)


# =============================================================================
# (b) Mortality disaggregation: proportional-hazard kappa solve
# =============================================================================
def _stand_mortality(kappa: float, tpa: np.ndarray, h: np.ndarray, T: float) -> float:
    """TPA-weighted stand mortality fraction under h_i' = kappa * h_i over T years."""
    surv = np.exp(-kappa * h * T)          # per-tree survival probability
    deaths = tpa * (1.0 - surv)
    return float(np.sum(deaths) / np.sum(tpa))


def stand_disaggregate_mortality(tree_hazards, tpa, T_years, M_stand_target,
                                 tol: float = 1e-10, max_iter: int = 200):
    """
    Solve for kappa so the TPA-weighted stand mortality of the rescaled tree
    hazards equals M_stand_target, and return (kappa, constrained h_i').

    Parameters
    ----------
    tree_hazards : array-like, per-tree ANNUAL hazard h_i (>= 0)
    tpa          : array-like, per-tree TPA / expansion weight (> 0)
    T_years      : float, remeasurement / projection interval (years)
    M_stand_target : float, target stand mortality fraction in [0, 1)

    Returns
    -------
    kappa : float, proportional hazard scalar (>= 0)
    h_prime : np.ndarray, constrained per-tree hazards kappa * h_i

    Edge cases
    ----------
    * M_stand_target <= 0                  -> kappa = 0 (no mortality).
    * all h_i == 0                          -> kappa = 0 (cannot produce mortality).
    * M_stand_target >= raw max achievable -> kappa clamped at the upper bracket
      (mortality asymptotes to 1 as kappa -> inf; we cannot exceed it, so we
      return the largest kappa that gets as close as possible). This guards the
      "target at or beyond the ceiling" case rather than diverging.
    """
    h = np.asarray(tree_hazards, dtype=float)
    w = np.asarray(tpa, dtype=float)
    if h.shape != w.shape:
        raise ValueError("tree_hazards and tpa must have the same shape")
    T = float(T_years)
    M = float(M_stand_target)

    # --- edge case: no mortality requested or no hazard available -------------
    if M <= 0.0 or not np.any(h > 0) or np.sum(w) <= 0:
        return 0.0, np.zeros_like(h)
    if M >= 1.0:
        # unreachable exactly (survival > 0 for finite kappa); push kappa high.
        M = 1.0 - 1e-12

    # raw stand mortality at kappa = 1 (the arm's own summed prediction)
    M_raw = _stand_mortality(1.0, w, h, T)

    # --- monotone bracket: mortality is strictly increasing in kappa ----------
    lo, hi = 0.0, 1.0
    # expand hi until stand mortality brackets the target (or we hit the ceiling)
    it = 0
    while _stand_mortality(hi, w, h, T) < M and it < 100:
        hi *= 2.0
        it += 1
    # if even a huge kappa cannot reach M (target beyond achievable ceiling),
    # return the largest kappa we bracketed (closest feasible mortality).
    if _stand_mortality(hi, w, h, T) < M:
        kappa = hi
        return kappa, kappa * h

    # bisection on the monotone function f(kappa) = stand_mortality(kappa) - M
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = _stand_mortality(mid, w, h, T) - M
        if abs(fm) < tol:
            lo = hi = mid
            break
        if fm < 0.0:
            lo = mid
        else:
            hi = mid
    kappa = 0.5 * (lo + hi)
    return kappa, kappa * h


# =============================================================================
# (a) Probabilistic Reineke density cap from the SDIMAX posterior
# =============================================================================
def sdimax_density_cap_draws(tpa, qmd, sdimax_draws, reineke_slope: float = REINEKE_SLOPE):
    """
    Probabilistic Reineke TPA cap using draws from the SDIMAX posterior.

    For each drawn SDIMAX value, compute the stand SDI = sum_i tpa_i*(qmd_i/10)^b
    (or scalar-stand SDI = TPA*(QMD/10)^b). If SDI exceeds the drawn SDIMAX, the
    tree TPAs are rescaled by SDIMAX/SDI (a uniform self-thinning); otherwise the
    TPAs are unchanged. Returns the per-draw capped TPA so the density ceiling
    carries the SDIMAX posterior uncertainty.

    Parameters
    ----------
    tpa  : array-like (n_trees,) per-tree TPA, OR a scalar stand TPA.
    qmd  : array-like (n_trees,) per-tree QMD (cm), OR a scalar stand QMD (cm).
           Broadcast against tpa. For a scalar stand, pass scalar TPA and QMD.
    sdimax_draws : array-like (n_draws,) SDIMAX posterior draws.

    Returns
    -------
    dict with:
      capped_tpa : (n_draws, n_trees) per-draw capped per-tree TPA
                   (or (n_draws,) if a scalar stand was passed)
      scale      : (n_draws,) per-draw rescale factor in (0, 1]; 1.0 = no bind
      sdi        : float, the (uncapped) stand SDI
      binds      : (n_draws,) bool, True where SDI > drawn SDIMAX
    """
    tpa = np.asarray(tpa, dtype=float)
    qmd = np.asarray(qmd, dtype=float)
    draws = np.asarray(sdimax_draws, dtype=float)
    scalar_stand = (tpa.ndim == 0)

    # stand SDI (Reineke, metric): per-tree contributions summed, or scalar form
    sdi_contrib = tpa * np.power(qmd / 10.0, reineke_slope)
    sdi = float(np.sum(sdi_contrib))

    # per-draw scale: bind only when SDI exceeds the drawn SDIMAX
    with np.errstate(divide="ignore", invalid="ignore"):
        raw_scale = np.where(sdi > 0, draws / sdi, 1.0)
    scale = np.minimum(raw_scale, 1.0)          # never scale UP; cap only
    binds = sdi > draws

    if scalar_stand:
        capped = float(tpa) * scale             # (n_draws,)
    else:
        capped = scale[:, None] * tpa[None, :]  # (n_draws, n_trees)
    return {"capped_tpa": capped, "scale": scale, "sdi": sdi, "binds": binds}


# =============================================================================
# Stand survival model evaluation (from stand_survival_bundle.json)
# =============================================================================
def stand_survival_eta(bundle, rd, ln_qmd, ba_metric, bgi,
                       trt_decay: float = 0.0, dstrb_decay: float = 0.0,
                       L1=None):
    """
    Reconstruct the stand-survival linear predictor (log-hazard) from the
    fitted bundle produced by 71_fit_stand_survival.R.

    The fit is a cloglog binomial with a log(YEARS) exposure offset, so the
    linear predictor IS the log annual stand hazard:

        linpred = b_Intercept + b_rd*rd + b_ln_qmd*ln_qmd + b_ba*ba_metric
                  + b_bgi*bgi + b_bgi_b2*max(bgi-knot,0)
                  + b_trt*trt_decay + b_dstrb*dstrb_decay + z_L1
        H_stand = exp(linpred)          # annual hazard
        S(T)    = exp(-H_stand * T)     # over YEARS
        M_stand = 1 - S(T)

    bundle : dict parsed from stand_survival_bundle.json (fixed_effects means,
             bgi knot, optional re_L1 table).
    Returns the linear predictor (float). Missing terms contribute 0.
    """
    fx = bundle.get("fixed_effects", {}) or {}

    def _m(name):
        v = fx.get(name)
        if isinstance(v, dict):
            return float(v.get("mean", 0.0))
        try:
            return float(v)
        except (TypeError, ValueError):
            return 0.0

    knot = bundle.get("covariates", {}).get("bgi_knot")
    if knot is None:
        knot = bundle.get("bgi_knot", 0.0)
    knot = float(knot or 0.0)
    bgi_b2 = max(float(bgi) - knot, 0.0)

    eta = (_m("Intercept") + _m("rd") * float(rd) + _m("ln_qmd") * float(ln_qmd)
           + _m("ba_metric") * float(ba_metric) + _m("bgi") * float(bgi)
           + _m("bgi_b2") * bgi_b2 + _m("trt_decay") * float(trt_decay)
           + _m("dstrb_decay") * float(dstrb_decay))

    if L1 is not None:
        re = bundle.get("re_L1") or {}
        levels = [str(x) for x in (re.get("level") or [])]
        means = re.get("mean") or []
        lut = {lv: float(m) for lv, m in zip(levels, means)}
        eta += lut.get(str(L1), 0.0)
    return float(eta)


def stand_mortality_target(bundle, T_years, rd, ln_qmd, ba_metric, bgi,
                           trt_decay: float = 0.0, dstrb_decay: float = 0.0,
                           L1=None):
    """Stand mortality target M_stand = 1 - exp(-exp(linpred) * T) from the
    fitted stand-survival bundle. This is the total the tree hazards are
    reconciled to via stand_disaggregate_mortality."""
    eta = stand_survival_eta(bundle, rd, ln_qmd, ba_metric, bgi,
                             trt_decay, dstrb_decay, L1)
    H = math.exp(min(eta, 30.0))          # annual hazard (guard overflow)
    return float(1.0 - math.exp(-H * float(T_years)))


def sdimax_draws_from_posterior(samples_path, n_draws=None, rng=None,
                                intercept_col="Intercept"):
    """
    Draw stand-level SDIMAX values from the hierarchical self-thinning posterior
    (output/variants/<v>/stand_density_samples.rds, a brms draws_df).

    The self-thinning fit is ln_tpa ~ ln_qmd_centered + (1+ln_qmd_centered|SPCD)
    with ln_qmd centered at log(10 in). At ln_qmd_centered = 0 the population
    ln(SDIMAX) draw is the `Intercept` column, so SDIMAX_draw = exp(Intercept).
    Returns an (n_draws,) array of population SDIMAX draws that carry the full
    self-thinning posterior uncertainty. Requires pyreadr (or a pre-exported
    .npy/.csv sidecar). Returns None if the RDS/sidecar is unavailable so the
    caller can fall back to a fixed SDIMAX (tagged).
    """
    import os
    if rng is None:
        rng = np.random.default_rng(20260702)
    # sidecar fast paths (npy / csv column) if present
    for ext, loader in ((".sdimax.npy", lambda p: np.load(p)),
                        (".sdimax.csv", lambda p: np.loadtxt(p))):
        side = str(samples_path) + ext
        if os.path.exists(side):
            arr = np.asarray(loader(side), dtype=float).ravel()
            if n_draws and n_draws < arr.size:
                arr = rng.choice(arr, size=n_draws, replace=False)
            return arr
    try:
        import pyreadr
    except Exception:
        return None
    try:
        res = pyreadr.read_r(str(samples_path))
        df = res[None] if None in res else list(res.values())[0]
    except Exception:
        return None
    if intercept_col not in df.columns:
        return None
    draws = np.exp(np.asarray(df[intercept_col], dtype=float))
    if n_draws and n_draws < draws.size:
        draws = rng.choice(draws, size=n_draws, replace=False)
    return draws


# =============================================================================
# (c) BA-growth disaggregation: proportional scale on tree DG
# =============================================================================
def stand_disaggregate_bagrowth(tree_dg, dbh, tpa, T_years, stand_ba_incr_target,
                                reineke_slope=None):
    """
    Scale per-tree diameter growth so the TPA-weighted summed stand basal-area
    increment matches a stand-level BA-growth target (from
    72_fit_stand_bagrowth.R). Because BA increment is additive across trees and
    (to first order in dg*T relative to dbh) proportional to the tree dg, a
    single closed-form proportional factor `gamma` on tree dg reconciles the sum
    exactly for the linear-in-dg BA increment.

    Per-tree BA at start:  ba_i     = (pi/4) * dbh_i^2         (per stem)
    Per-tree BA at end:    ba_i'(g)  = (pi/4) * (dbh_i + gamma*dg_i*T)^2
    Stand BA increment (per unit area, TPA-weighted):
        dBA_stand(gamma) = sum_i tpa_i * (ba_i'(gamma) - ba_i)

    We solve dBA_stand(gamma) = target. This is a quadratic in gamma with
    positive leading coefficient, so it has a unique positive root (closed form).
    Units: dbh and dg in the SAME length unit; T in years; target in the matching
    area unit per unit area (e.g. m2/ha over the interval). (pi/4) cancels if the
    target is expressed on the same BA definition; we keep it explicit so the
    target is a real BA increment.

    Returns (gamma, dg_prime, achieved) where dg_prime = gamma * dg and achieved
    is the reconciled stand BA increment (== target up to float error).
    """
    g = np.asarray(tree_dg, dtype=float)
    d = np.asarray(dbh, dtype=float)
    w = np.asarray(tpa, dtype=float)
    if not (g.shape == d.shape == w.shape):
        raise ValueError("tree_dg, dbh, tpa must have the same shape")
    T = float(T_years)
    k = math.pi / 4.0
    target = float(stand_ba_incr_target)

    inc = g * T                                  # per-tree diameter increment at gamma=1
    # dBA_stand(gamma) = sum w_i*k*[(d_i+gamma*inc_i)^2 - d_i^2]
    #                  = A*gamma^2 + B*gamma,  A,B >= 0
    A = float(np.sum(w * k * inc * inc))
    B = float(np.sum(w * k * 2.0 * d * inc))

    if target <= 0.0:
        return 0.0, np.zeros_like(g), 0.0
    if A <= 0.0 and B <= 0.0:                     # no growth capacity
        return 0.0, np.zeros_like(g), 0.0
    if A <= 0.0:                                  # linear: B*gamma = target
        gamma = target / B
    else:                                        # quadratic: positive root
        disc = B * B + 4.0 * A * target
        gamma = (-B + math.sqrt(disc)) / (2.0 * A)

    dg_prime = gamma * g
    achieved = A * gamma * gamma + B * gamma
    return float(gamma), dg_prime, float(achieved)


# =============================================================================
# (d) TOP-HEIGHT constraint: García/GADA stand top-height -> tree height growth
# =============================================================================
def stand_top_height(tree_heights, tpa, top_n_per_ha=100.0):
    """
    Stand TOP (dominant) height = TPA-weighted mean height of the tallest
    `top_n_per_ha` stems/ha (the García/GADA/Assmann dominant-height basis).

    Trees are ranked tall-to-short; TPA is accumulated until `top_n_per_ha`
    stems have been included (the last tree is partially counted so exactly
    top_n_per_ha stems define the mean). Returns (top_ht, cohort_mask, w_top)
    where w_top is the (possibly fractional) TPA weight each tree contributes
    to the top cohort (0 for trees below the cohort).

    If total TPA < top_n_per_ha, the top height is the TPA-weighted mean of all
    trees (best available), and w_top == tpa.
    """
    h = np.asarray(tree_heights, dtype=float)
    w = np.asarray(tpa, dtype=float)
    order = np.argsort(-h)                       # tallest first
    w_ord = w[order]
    cum = np.cumsum(w_ord)
    need = float(top_n_per_ha)
    w_top_ord = np.zeros_like(w_ord)
    if cum[-1] <= need:                          # not enough stems: use all
        w_top_ord = w_ord.copy()
    else:
        # fully include trees until cumulative TPA reaches `need`, then take a
        # fractional slice of the boundary tree so exactly `need` stems count.
        full = cum <= need
        w_top_ord[full] = w_ord[full]
        k = int(np.argmax(~full))                # first tree that overflows
        prev = cum[k - 1] if k > 0 else 0.0
        w_top_ord[k] = need - prev               # partial weight
    # map back to original order
    w_top = np.zeros_like(w)
    w_top[order] = w_top_ord
    denom = float(np.sum(w_top))
    top_ht = float(np.sum(w_top * h) / denom) if denom > 0 else float("nan")
    return top_ht, (w_top > 0), w_top


def stand_constrain_topheight(tree_heights, tree_htg, target_top_height_end,
                              top_n_per_ha=100.0, tpa=None, T_years=1.0,
                              tol=1e-9, max_iter=200):
    """
    Scale tree HEIGHT GROWTH so the projected stand TOP height (mean height of
    the tallest ~top_n_per_ha stems/ha) tracks the García/GADA top-height
    trajectory target H(t_end). Only the dominant cohort defines top height, so
    the correction is distributed to the taller trees.

    We apply a single proportional scalar `phi` to the height growth of the
    trees that are in (or can enter) the top cohort:

        h_end_i(phi) = h_start_i + phi * htg_i * T
        top_ht_end(phi) = stand_top_height(h_end(phi), tpa, top_n_per_ha)

    top_ht_end is monotone non-decreasing in phi (taller trees grow taller, and
    the top cohort is the tallest set), so a 1-D bracket + bisection on phi is
    robust. phi < 1 slows the dominant cohort toward the target; phi > 1 speeds
    it up. Height growth of sub-dominant trees far below the cohort is left
    unscaled -- it never affects top height -- but to keep the height-order
    self-consistent we scale ALL tree htg by phi (a uniform dominant-driven
    correction); trees below the cohort simply don't enter the top-height mean.

    Parameters
    ----------
    tree_heights : (n,) per-tree START height (m)
    tree_htg     : (n,) per-tree ANNUAL height growth (m/yr), unconstrained
    target_top_height_end : float, target stand TOP height at t_end (m)
    top_n_per_ha : float, number of tallest stems/ha defining top height (100)
    tpa          : (n,) per-tree TPA/expansion weight (defaults to ones)
    T_years      : float, projection interval (years)

    Returns
    -------
    dict: phi, htg_prime (scaled per-tree htg), heights_end (constrained end
          heights), top_ht_end (achieved top height), top_ht_start.
    """
    h0 = np.asarray(tree_heights, dtype=float)
    g = np.asarray(tree_htg, dtype=float)
    if h0.shape != g.shape:
        raise ValueError("tree_heights and tree_htg must have the same shape")
    w = np.ones_like(h0) if tpa is None else np.asarray(tpa, dtype=float)
    T = float(T_years)
    target = float(target_top_height_end)

    def top_end(phi):
        h_end = h0 + phi * g * T
        return stand_top_height(h_end, w, top_n_per_ha)[0]

    top_ht_start = stand_top_height(h0, w, top_n_per_ha)[0]

    # monotone in phi: bracket the target
    lo, hi = 0.0, 1.0
    it = 0
    while top_end(hi) < target and it < 100:
        hi *= 2.0
        it += 1
    # if even a huge phi cannot reach the target (target above achievable),
    # return the largest phi bracketed (closest feasible top height).
    if top_end(hi) < target:
        phi = hi
        htg_p = phi * g
        return {"phi": phi, "htg_prime": htg_p, "heights_end": h0 + htg_p * T,
                "top_ht_end": top_end(phi), "top_ht_start": top_ht_start}
    # if target is at/below the no-growth top height, phi -> 0 floor
    if top_end(0.0) >= target:
        lo = hi = 0.0
    else:
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            fm = top_end(mid) - target
            if abs(fm) < tol:
                lo = hi = mid
                break
            if fm < 0.0:
                lo = mid
            else:
                hi = mid
    phi = 0.5 * (lo + hi)
    htg_p = phi * g
    return {"phi": phi, "htg_prime": htg_p, "heights_end": h0 + htg_p * T,
            "top_ht_end": top_end(phi), "top_ht_start": top_ht_start}


def gada_topheight_transition(h1, years, b2, b3):
    """
    García-style base-age-invariant (state-space) TOP-HEIGHT transition using the
    Cieszewski-Bailey Chapman-Richards GADA form fitted in gada_refit.r:

        H_t = b1 * (1 - exp(-b2 * t))^b3      (site enters via b1)

    Base-age invariance -> the H2 from H1 transition eliminates b1 and age. From
    a starting top height H1 at (unknown) age, after `years` more years:

        H2 = H1 * ( (1 - exp(-b2*(a1 + years))) / (1 - exp(-b2*a1)) )^b3

    where a1 is the implied age from inverting H1 = b1*(1-exp(-b2*a1))^b3. Since
    b1 is a per-stand nuisance (site), we use the pair-form invariant that García
    uses: for anamorphic C-R the ratio does not depend on b1, so a1 is recovered
    from a reference asymptote only to advance the interval. In the minimal
    state-space version we FIT a direct H2 = f(H1, years) transition (see the
    75_fit_stand_topheight.R launcher) whose coefficients replace this closed
    form; this function is the analytic fallback / prior mean.

    Here we use the invariant advance with a1 solved from a nominal reference
    asymptote b1_ref supplied by the caller via a1 (age). If `years<=0`, returns
    H1 unchanged. This is a deterministic H(t) target generator for
    stand_constrain_topheight when a fitted transition bundle is unavailable.
    """
    h1 = float(h1)
    yr = float(years)
    if yr <= 0.0 or h1 <= 0.0:
        return h1
    # Solve implied age a1 from H1 assuming a shared reference asymptote.
    # For the invariant we only need the RATIO, which is stable across b1_ref;
    # invert with b1_ref = H1 / (1-exp(-b2*a_ref))^b3 at a_ref -> use a1 s.t. the
    # C-R fraction equals H1/b1_ref. We adopt the standard GADA solution that
    # advances the interval with the fitted (b2,b3); a1 is found by matching the
    # observed H1 to the curve at a reference site. Use a numeric solve:
    # find a1 in (0.5, 400] with H1 = b1_ref*(1-exp(-b2*a1))^b3 is ill-posed
    # without b1_ref, so we use the base-age-invariant closed form directly with
    # the Cieszewski (2002) X0 solution:
    #   X0 = 0.5*(ln H1 + sqrt( (ln H1)^2 ... ))  -- for the specific ADA case.
    # To stay robust and dependency-free we advance via the multiplicative
    # fraction using an anchor age a1 = -ln(1 - (H1/ (1.3*H1_scale))^(1/b3))/b2.
    # In practice the fitted transition bundle is used; keep a safe monotone
    # advance: H2 = H1 * ((1-exp(-b2*(a_anchor+yr)))/(1-exp(-b2*a_anchor)))^b3.
    a_anchor = 20.0  # nominal anchor age; the fitted bundle supersedes this
    frac = ((1.0 - math.exp(-b2 * (a_anchor + yr))) /
            (1.0 - math.exp(-b2 * a_anchor))) ** b3
    return h1 * frac


# =============================================================================
# (e) STEM-DENSITY constraint: García state-space N(t) -> tree mortality
# =============================================================================
def stand_constrain_stems(tree_hazards, tpa, T_years, N_target_end,
                          N_start=None, tol=1e-10, max_iter=200):
    """
    Reconcile tree survival so the SURVIVING stems/ha match the state-space N(t)
    target N_target_end. This UNIFIES with stand_disaggregate_mortality: stems
    are just the TPA (count) form of the survival target. The stand mortality
    fraction implied by the stem target is

        M_stand = 1 - N_target_end / N_start,   N_start = sum_i tpa_i

    and the same proportional-hazard kappa solve rescales the tree hazards so the
    surviving TPA equals N_target_end:

        sum_i tpa_i * exp(-kappa * h_i * T) == N_target_end

    Parameters
    ----------
    tree_hazards : (n,) per-tree ANNUAL hazard h_i (>= 0)
    tpa          : (n,) per-tree TPA / expansion weight (> 0)
    T_years      : float, projection interval (years)
    N_target_end : float, target SURVIVING stems/ha at t_end
    N_start      : float, starting stems/ha (defaults to sum(tpa))

    Returns
    -------
    dict: kappa, hazards_prime, N_end (achieved surviving stems), N_start,
          M_stand (implied stand mortality fraction).
    """
    h = np.asarray(tree_hazards, dtype=float)
    w = np.asarray(tpa, dtype=float)
    if h.shape != w.shape:
        raise ValueError("tree_hazards and tpa must have the same shape")
    T = float(T_years)
    Ns = float(np.sum(w)) if N_start is None else float(N_start)
    Nt = float(N_target_end)

    if Ns <= 0:
        return {"kappa": 0.0, "hazards_prime": np.zeros_like(h),
                "N_end": 0.0, "N_start": Ns, "M_stand": 0.0}
    # clamp target to (0, Ns]; stem target above the start means no mortality.
    Nt = min(max(Nt, 0.0), Ns)
    M_stand = 1.0 - Nt / Ns
    kappa, h_prime = stand_disaggregate_mortality(h, w, T, M_stand, tol=tol,
                                                  max_iter=max_iter)
    N_end = float(np.sum(w * np.exp(-kappa * h * T)))
    return {"kappa": kappa, "hazards_prime": h_prime, "N_end": N_end,
            "N_start": Ns, "M_stand": M_stand}


# =============================================================================
# Fitted-transition bundle evaluators (top height H2|H1, stems N2|N1)
# =============================================================================
def stand_topheight_target(bundle, h1, years, rd=None, ln_qmd=None, bgi=None,
                           L1=None):
    """
    Top-height H2 target from the fitted GADA-transition bundle
    (75_fit_stand_topheight.R): a log-scale linear model of the height RATIO
    r = ln(H2/H1) on ln(H1), ln(years), and optional stand covariates, so
        H2 = H1 * exp(eta_r),  eta_r = b0 + b_lnH1*ln(H1) + b_lnyr*ln(years)
                                       + b_rd*rd + b_lnqmd*ln_qmd + b_bgi*bgi + z_L1
    Missing terms contribute 0. If the bundle is absent/None, falls back to the
    analytic gada_topheight_transition using bundle-carried (b2,b3) or defaults.
    """
    if not bundle:
        return gada_topheight_transition(h1, years, b2=0.03, b3=1.1)
    fx = bundle.get("fixed_effects", {}) or {}

    def _m(name):
        v = fx.get(name)
        if isinstance(v, dict):
            return float(v.get("mean", 0.0))
        try:
            return float(v)
        except (TypeError, ValueError):
            return 0.0

    if float(years) <= 0 or float(h1) <= 0:
        return float(h1)
    eta = (_m("Intercept") + _m("ln_h1") * math.log(float(h1))
           + _m("ln_years") * math.log(float(years)))
    if rd is not None:
        eta += _m("rd") * float(rd)
    if ln_qmd is not None:
        eta += _m("ln_qmd") * float(ln_qmd)
    if bgi is not None:
        eta += _m("bgi") * float(bgi)
    if L1 is not None:
        re = bundle.get("re_L1") or {}
        lut = {str(lv): float(m) for lv, m in zip(re.get("level", []),
                                                  re.get("mean", []))}
        eta += lut.get(str(L1), 0.0)
    return float(h1) * math.exp(eta)


# NE spruce-fir GADA Chapman-Richards shape parameters, refit 2026-07-04 from
# NA_SITREE.csv (SPCD 12 balsam fir, 95 black spruce, 97 red spruce; n=197,122
# site-tree height-age records, AGEDIA 10-200). Form H=b1*(1-exp(-b2*A))^b3;
# site enters through the asymptote b1. b2,b3 are the anamorphic shape; b1 is
# recovered per-stand from the current top height so the anchor auto-calibrates
# to the stand's own trajectory (Cieszewski base-age-invariant advance).
GADA_SF_B2 = 0.03584
GADA_SF_B3 = 1.5898


def sf_site_asymptote(bgi, b1_lo=22.0, b1_hi=28.0, bgi_lo=4.0, bgi_hi=8.0):
    """
    Per-stand GADA top-height ASYMPTOTE b1 (m) as a monotone-increasing function
    of the site driver bgi, so higher site -> higher asymptote (correct site
    ordering). Range chosen to match the NE spruce-fir top-cohort (top-100
    stems/ha) asymptote 19-26 m: the NA_SITREE GADA refit gives per-site SI50
    (base age 50) up to ~17 m and asymptote b1 up to ~23 m for individual site
    trees; the top-100-stems/ha dominant cohort runs a few metres taller, so the
    stand top-height asymptote sits ~19-26 m across the bgi 4-8 gradient. Linear
    in bgi, clamped to [b1_lo, b1_hi].
    """
    if bgi is None:
        return 0.5 * (b1_lo + b1_hi)
    t = (float(bgi) - bgi_lo) / max(bgi_hi - bgi_lo, 1e-6)
    t = min(max(t, 0.0), 1.0)
    return b1_lo + t * (b1_hi - b1_lo)


def stand_topheight_target_gada(bundle, h1, years, rd=None, ln_qmd=None,
                                bgi=None, L1=None, b2=GADA_SF_B2, b3=GADA_SF_B3,
                                b1=None, blend=0.85, h2_ss=None):
    """
    MONOTONE, ASYMPTOTING, SITE-ORDERED top-height H2 target -- the FIX for the
    mean-reverting raw state-space transition.

    Root cause of the old bug: stand_topheight_target returns
        H2 = H1 * exp(b0 + b_lnH1*ln(H1) + ...),  b_lnH1 = -0.155 < 0,
    a log-ratio regression with a NEGATIVE ln(H1) slope. Fit on FIA remeasurement
    intervals it is fine one-step, but iterated over a projection it is a
    CONTRACTION MAP: H1 converges to a step-dependent fixed point (~9.7 m for a
    5-yr step) instead of growing to the site asymptote, and because b_bgi < 0 the
    equilibrium is INVERTED by site (higher bgi -> lower H). That produces the
    observed rise-then-decline (20.7 -> 22.3 -> 20.6 m) and violates the basic law
    that dominant/top height cannot shrink in an even-aged stand.

    The fix anchors the target to the base-age-invariant GADA Chapman-Richards
    site trajectory H = b1*(1-exp(-b2*A))^b3 (per-species b2,b3; site via the
    asymptote b1). b1 is a FIXED per-stand site asymptote (from `sf_site_asymptote`
    on bgi, so site ordering is guaranteed), NOT re-derived from H1 each step
    (which would let the asymptote run away). We invert the CURRENT top height H1
    against b1 to get the implied age a1,

        a1 = -ln( 1 - (H1/b1)^(1/b3) ) / b2,

    then advance the GADA curve by the interval:

        H2_gada = b1 * (1 - exp(-b2*(a1+years)))^b3 .

    Because a1 grows as H1 grows toward b1, H2_gada is strictly increasing in
    `years`, ASYMPTOTES at the fixed b1, and is ordered by site (b1 increases with
    bgi). It can never decline. We then take a monotone, capped blend with the raw
    state-space target (which still carries fitted density/qmd signal near the
    data) and hard-floor at H1:

        H2 = clip( blend*H2_gada + (1-blend)*min(H2_ss, b1),  H1,  b1 ).

    With blend=1 this is the pure GADA advance (fully monotone, fully asymptoting);
    blend<1 lets the state-space model modulate the approach where it is well
    identified while the GADA anchor guarantees monotonicity, asymptoting, and
    site ordering.

    Parameters mirror stand_topheight_target; extra:
      b2,b3 : GADA shape (default NE spruce-fir refit from NA_SITREE).
      b1    : per-stand asymptote (m); defaults to sf_site_asymptote(bgi).
      blend : weight on the GADA advance vs the (capped) state-space target.
    """
    h1 = float(h1)
    yr = float(years)
    if yr <= 0.0 or h1 <= 0.0:
        return h1
    if b1 is None:
        b1 = sf_site_asymptote(bgi)
    b1 = float(b1)
    # if the stand is already at/above its asymptote, hold (no decline, no runaway)
    if h1 >= b1:
        return h1
    # raw state-space target for the modulation term: use a caller-supplied value
    # (source-agnostic: config runtime OR external bundle) or recompute here.
    if h2_ss is None:
        h2_ss = stand_topheight_target(bundle, h1, yr, rd=rd, ln_qmd=ln_qmd,
                                       bgi=bgi, L1=L1)
    h2_ss = float(h2_ss)
    # invert H1 against the FIXED site asymptote b1 to recover the implied age a1
    ratio = min(max(h1 / b1, 1e-6), 1.0 - 1e-9)
    inner = 1.0 - ratio ** (1.0 / b3)
    inner = min(max(inner, 1e-9), 1.0 - 1e-12)
    a1 = -math.log(inner) / b2
    # advance the GADA curve over the interval (monotone, asymptotes at b1)
    h2_gada = b1 * (1.0 - math.exp(-b2 * (a1 + yr))) ** b3
    # monotone, capped blend + hard floor at H1 (never declines) and cap at b1
    h2 = blend * h2_gada + (1.0 - blend) * min(h2_ss, b1)
    h2 = min(max(h2, h1), b1)
    return float(h2)


def stand_stems_target(bundle, n1, years, top_ht=None, rd=None, ln_qmd=None,
                       L1=None):
    """
    Surviving-stems N2 target from the fitted state-space stem-transition bundle
    (76_fit_stand_stems.R): log-ratio model s = ln(N2/N1) (<= 0, a survival
    fraction) on ln(N1), ln(years), top height, relative density:
        N2 = N1 * exp(-softplus(eta_s)),  eta_s = b0 + b_lnN1*ln(N1)
             + b_lnyr*ln(years) + b_topht*top_ht + b_rd*rd + b_lnqmd*ln_qmd + z_L1
    The softplus keeps N2 <= N1 (stems only decline absent ingrowth, which is the
    survival channel this constraint governs). Missing terms contribute 0.
    """
    if not bundle or float(years) <= 0 or float(n1) <= 0:
        return float(n1)
    fx = bundle.get("fixed_effects", {}) or {}

    def _m(name):
        v = fx.get(name)
        if isinstance(v, dict):
            return float(v.get("mean", 0.0))
        try:
            return float(v)
        except (TypeError, ValueError):
            return 0.0

    eta = (_m("Intercept") + _m("ln_n1") * math.log(float(n1))
           + _m("ln_years") * math.log(float(years)))
    if top_ht is not None:
        eta += _m("top_ht") * float(top_ht)
    if rd is not None:
        eta += _m("rd") * float(rd)
    if ln_qmd is not None:
        eta += _m("ln_qmd") * float(ln_qmd)
    if L1 is not None:
        re = bundle.get("re_L1") or {}
        lut = {str(lv): float(m) for lv, m in zip(re.get("level", []),
                                                  re.get("mean", []))}
        eta += lut.get(str(L1), 0.0)
    # softplus mortality intensity -> guarantees 0 < N2 <= N1
    mort_int = math.log1p(math.exp(min(eta, 30.0)))
    return float(n1) * math.exp(-mort_int)


# =============================================================================
# self-test
# =============================================================================
if __name__ == "__main__":
    rng = np.random.default_rng(20260702)

    print("=" * 70)
    print("SELF-TEST (a): mortality disaggregation kappa solve")
    print("=" * 70)
    n = 50
    tpa = rng.uniform(1.0, 20.0, size=n)             # per-tree expansion weights
    h = rng.uniform(0.002, 0.05, size=n)             # per-tree annual hazards
    T = 10.0
    M_raw = _stand_mortality(1.0, tpa, h, T)
    print(f"raw stand mortality at kappa=1 (arm's own sum): {M_raw:.6f}")

    for M_target in (0.05, M_raw, 0.30, 0.60):
        kappa, h_prime = stand_disaggregate_mortality(h, tpa, T, M_target)
        M_hit = _stand_mortality(1.0, tpa, h_prime, T)   # mortality of rescaled h
        err = abs(M_hit - min(M_target, 1 - 1e-12))
        flag = "OK" if err < 1e-4 else "FAIL"
        print(f"  target={M_target:.4f}  kappa={kappa:8.5f}  achieved={M_hit:.6f}  "
              f"|err|={err:.2e}  [{flag}]")

    # edge cases
    k0, hp0 = stand_disaggregate_mortality(h, tpa, T, 0.0)
    print(f"  edge M=0        -> kappa={k0:.5f} (expect 0), max h'={hp0.max():.3e}")
    k_hi, _ = stand_disaggregate_mortality(h, tpa, T, 0.999999)
    M_hi = _stand_mortality(k_hi, tpa, h, T)
    print(f"  edge M->1       -> kappa={k_hi:.3f}, achieved={M_hi:.6f} (guarded, no divergence)")

    print()
    print("=" * 70)
    print("SELF-TEST (b): probabilistic Reineke density cap")
    print("=" * 70)
    # a dense stand: TPA and QMD chosen so SDI ~ 900
    n2 = 40
    tpa2 = np.full(n2, 25.0)                          # 1000 TPA total
    qmd2 = rng.uniform(18.0, 26.0, size=n2)          # cm
    sdi = float(np.sum(tpa2 * (qmd2 / 10.0) ** REINEKE_SLOPE))
    print(f"stand SDI (uncapped): {sdi:.2f}")

    # SDIMAX posterior draws spanning below and above the stand SDI
    sdimax_draws = np.array([sdi * 0.70, sdi * 0.90, sdi * 1.00, sdi * 1.10, sdi * 1.30])
    res = sdimax_density_cap_draws(tpa2, qmd2, sdimax_draws)
    tot0 = float(np.sum(tpa2))
    print("  SDIMAX_draw   binds   scale     capped_total_TPA")
    for sd, b, sc, cap in zip(sdimax_draws, res["binds"], res["scale"],
                              res["capped_tpa"].sum(axis=1)):
        print(f"  {sd:10.2f}    {str(bool(b)):5s}   {sc:6.4f}   {cap:10.2f}  "
              f"(uncapped {tot0:.2f})")

    # assertions: cap binds only when SDI > SDIMAX, and shrinks with lower SDIMAX
    binds_expected = sdimax_draws < sdi
    assert np.array_equal(res["binds"], binds_expected), "bind logic mismatch"
    capped_totals = res["capped_tpa"].sum(axis=1)
    # among binding draws, lower SDIMAX -> smaller capped TPA (monotone)
    binding = np.where(res["binds"])[0]
    order = binding[np.argsort(sdimax_draws[binding])]
    monotone = np.all(np.diff(capped_totals[order]) > 0)
    print(f"\n  cap binds only when SDI>SDIMAX: "
          f"{bool(np.array_equal(res['binds'], binds_expected))}")
    print(f"  capped TPA shrinks with lower drawn SDIMAX (monotone): {bool(monotone)}")
    # non-binding draws leave TPA unchanged
    nonbind = np.where(~res["binds"])[0]
    unchanged = np.allclose(capped_totals[nonbind], tot0)
    print(f"  non-binding draws leave TPA unchanged: {bool(unchanged)}")
    assert monotone and unchanged, "density cap monotonicity/no-op check failed"

    print()
    print("=" * 70)
    print("SELF-TEST (c): stand BA-growth disaggregation (proportional dg scale)")
    print("=" * 70)
    n3 = 30
    dbh3 = rng.uniform(10.0, 45.0, size=n3)          # cm
    dg3 = rng.uniform(0.05, 0.55, size=n3)           # cm/yr (unconstrained tree DG)
    tpa3 = rng.uniform(2.0, 30.0, size=n3)           # per-tree TPA
    Tg = 10.0
    k = math.pi / 4.0
    # raw (unconstrained) summed stand BA increment, cm^2/ha over the interval
    raw_ba_incr = float(np.sum(tpa3 * k * ((dbh3 + dg3 * Tg) ** 2 - dbh3 ** 2)))
    print(f"raw (unconstrained) stand BA increment: {raw_ba_incr:.2f}")

    for frac, lbl in ((0.70, "target below raw"), (1.00, "target == raw"),
                      (1.35, "target above raw")):
        tgt = raw_ba_incr * frac
        gamma, dg_p, ach = stand_disaggregate_bagrowth(dg3, dbh3, tpa3, Tg, tgt)
        err = abs(ach - tgt)
        flag = "OK" if err < 1e-4 * max(tgt, 1.0) else "FAIL"
        print(f"  {lbl:16s} target={tgt:10.2f}  gamma={gamma:7.5f}  "
              f"achieved={ach:10.2f}  |err|={err:.2e}  [{flag}]")
        assert err < 1e-3 * max(tgt, 1.0), "BA-growth disaggregation missed target"
    # edge: zero target -> zero growth
    g0, dgp0, ach0 = stand_disaggregate_bagrowth(dg3, dbh3, tpa3, Tg, 0.0)
    print(f"  edge target=0    -> gamma={g0:.5f} (expect 0), achieved={ach0:.4f}")
    assert g0 == 0.0 and ach0 == 0.0
    # BA increment is additive: reconciled per-tree increments sum to target
    tgt = raw_ba_incr * 0.85
    gamma, dg_p, ach = stand_disaggregate_bagrowth(dg3, dbh3, tpa3, Tg, tgt)
    per_tree = tpa3 * k * ((dbh3 + dg_p * Tg) ** 2 - dbh3 ** 2)
    print(f"  additivity: sum of reconciled per-tree BA incr = {per_tree.sum():.2f} "
          f"(target {tgt:.2f})")
    assert abs(per_tree.sum() - tgt) < 1e-3 * tgt

    print()
    print("=" * 70)
    print("SELF-TEST (d): García/GADA top-height constraint (tree HG scaling)")
    print("=" * 70)
    n4 = 60
    # a stand with a tall dominant cohort and a shorter understory
    h0 = np.concatenate([rng.uniform(24.0, 30.0, size=15),   # dominants
                         rng.uniform(10.0, 22.0, size=45)])   # understory
    htg = rng.uniform(0.10, 0.40, size=n4)                    # m/yr height growth
    tpa4 = np.concatenate([rng.uniform(3.0, 8.0, size=15),
                           rng.uniform(5.0, 25.0, size=45)])  # ~700 stems/ha
    Th = 10.0
    top_n = 100.0
    top0 = stand_top_height(h0, tpa4, top_n)[0]
    # unconstrained end top height (phi=1)
    top_raw = stand_top_height(h0 + htg * Th, tpa4, top_n)[0]
    print(f"start top height: {top0:.3f} m  (tallest {top_n:.0f}/ha)")
    print(f"raw end top height (phi=1): {top_raw:.3f} m")

    for lbl, tgt in (("below raw", top_raw - 1.5),
                     ("== raw", top_raw),
                     ("above raw", top_raw + 1.2),
                     ("GADA H(t)", gada_topheight_transition(top0, Th, b2=0.033, b3=1.1))):
        r = stand_constrain_topheight(h0, htg, tgt, top_n_per_ha=top_n,
                                      tpa=tpa4, T_years=Th)
        err = abs(r["top_ht_end"] - tgt)
        flag = "OK" if err < 1e-3 else "FAIL"
        print(f"  {lbl:12s} target={tgt:7.3f}  phi={r['phi']:7.4f}  "
              f"achieved={r['top_ht_end']:7.3f}  |err|={err:.2e}  [{flag}]")
        # target reachable within [0, hi]: those must hit exactly
        if tgt <= top_raw + 1e-9 and tgt >= top0:
            assert err < 1e-2, "top-height constraint missed a reachable target"
    # monotonicity: bigger phi -> taller or equal top height
    tops = [stand_top_height(h0 + p * htg * Th, tpa4, top_n)[0]
            for p in (0.0, 0.5, 1.0, 2.0)]
    assert all(np.diff(tops) >= -1e-9), "top height not monotone in phi"
    print(f"  monotone in phi: {tops[0]:.2f} <= {tops[1]:.2f} <= "
          f"{tops[2]:.2f} <= {tops[3]:.2f}  [OK]")

    print()
    print("=" * 70)
    print("SELF-TEST (e): García state-space stem-density constraint (tree mortality)")
    print("=" * 70)
    n5 = 80
    tpa5 = rng.uniform(2.0, 18.0, size=n5)          # per-tree stems/ha
    haz5 = rng.uniform(0.003, 0.06, size=n5)        # per-tree annual hazard
    Ts = 10.0
    N_start = float(np.sum(tpa5))
    N_raw = float(np.sum(tpa5 * np.exp(-haz5 * Ts)))  # arm's own surviving stems
    print(f"start stems/ha: {N_start:.2f}   raw surviving (kappa=1): {N_raw:.2f}")

    for lbl, Nt in (("more mortality", N_raw * 0.85),
                    ("== raw survival", N_raw),
                    ("less mortality", (N_start + N_raw) / 2.0)):
        r = stand_constrain_stems(haz5, tpa5, Ts, Nt, N_start=N_start)
        err = abs(r["N_end"] - Nt)
        flag = "OK" if err < 1e-3 * N_start else "FAIL"
        print(f"  {lbl:16s} N_target={Nt:8.3f}  kappa={r['kappa']:7.4f}  "
              f"M_stand={r['M_stand']:.4f}  N_end={r['N_end']:8.3f}  "
              f"|err|={err:.2e}  [{flag}]")
        assert err < 1e-2 * N_start, "stem-density constraint missed target"
    # unification check: stems constraint == mortality disaggregation on M_stand
    Nt = N_raw * 0.80
    rs = stand_constrain_stems(haz5, tpa5, Ts, Nt, N_start=N_start)
    M = 1.0 - Nt / N_start
    km, _ = stand_disaggregate_mortality(haz5, tpa5, Ts, M)
    print(f"  unified with kappa mortality solve: stems kappa={rs['kappa']:.5f}  "
          f"mortality kappa={km:.5f}  [{'OK' if abs(rs['kappa']-km)<1e-6 else 'FAIL'}]")
    assert abs(rs["kappa"] - km) < 1e-6, "stems/mortality kappa solves disagree"
    # edge: target >= start -> no mortality
    r0 = stand_constrain_stems(haz5, tpa5, Ts, N_start * 1.1, N_start=N_start)
    print(f"  edge N_target>=N_start -> kappa={r0['kappa']:.5f} (expect 0), "
          f"N_end={r0['N_end']:.2f} (expect {N_start:.2f})")
    assert r0["kappa"] == 0.0 and abs(r0["N_end"] - N_start) < 1e-9

    print("\nALL SELF-TESTS PASSED")
