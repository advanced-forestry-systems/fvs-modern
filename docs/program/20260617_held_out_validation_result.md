# Held-out spatial-fold validation (red-team issue #1) — result
2026-06-17. Addresses the leakage objection: derive the calibration on a spatial calibration fold (fold A,
counties hashed to A), apply UNCHANGED to a spatially-held-out fold (fold B), report out-of-sample bias.
brms SDImax is plot-level (not fold-derived); ingrowth rate and BAIMULT are derived on fold A only.

## Out-of-sample result (fold B, held out), default -> calibrated bias %

| variant | nA/nB | fold-A rate %/dec | BAIMULT | QMD | BA | TPH |
|---|---|---|---|---|---|---|
| ne | 37/63 | 33.6 | 1.00 | +9.4 -> +3.9 | +7.1 -> +7.9 | -12.2 -> -4.7 |
| sn | 39/46 | 63.1 | 0.70 | +13.2 -> +1.1 | +16.0 -> +16.6 | -6.7 -> +16.1 |
| kt | 13/28 | 43.9 | 0.70 | +32.3 -> +11.8 | +14.6 -> +0.5 | -43.4 -> -30.3 |
| pn | 42/28 | 20.6 | 0.90 | +14.7 -> +10.0 | +13.3 -> +7.6 | -22.1 -> -27.4 |

## What generalizes and what does not

- SIZE LEVERS GENERALIZE. QMD bias is reduced out-of-sample in all four variants (e.g. sn +13.2 -> +1.1,
  kt +32.3 -> +11.8); BA is reduced or held in three of four (kt +14.6 -> +0.5, pn +13.3 -> +7.6). The
  brms SDImax + BAIMODT (diameter-growth) levers are not overfit: they transfer to spatially-held-out
  plots. This rebuts the leakage objection for the size/density-cap levers.
- THE FIXED INGROWTH RATE DOES NOT TRANSFER CLEANLY. TPH improves out-of-sample for ne (-12.2 -> -4.7)
  and kt (partial), but OVER-corrects for sn (-6.7 -> +16.1) and pn (worsens). The fold-A per-variant
  recruitment rate over/under-shoots fold B because recruitment is spatially variable and
  density-dependent. This confirms the red-team's medium issue and the integration roadmap item: the
  recruitment lever needs a DENSITY-DEPENDENT / site-resolved form, not a fixed per-variant %/decade.

## Implications

1. Report the calibration honestly as out-of-sample for the size levers (they hold) and as a prototype
   needing a better functional form for recruitment (it does not transfer).
2. Next: replace the fixed ingrowth rate with a density-dependent recruitment model (function of stand
   density / SDImax headroom and site productivity), then re-run the held-out validation; expect TPH to
   then transfer.
3. The fold sizes are modest (kt nA=13) and folds are heterogeneous (e.g. ne default QMD +4.8 on A vs
   +9.4 on B) — extend to all variants with larger samples and bootstrap CIs.

Data: held_out.csv. Method: held_out_validation.py (county-hash spatial folds, fold-A-derived ingrowth
rate + BAIMULT-min-|QMD| sweep, applied to fold B).
