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

    print("\nALL SELF-TESTS PASSED")
