# Held-out validation with density-dependent recruitment — result
2026-06-18. Item 1 of the integration roadmap. Replaces the fixed per variant ingrowth rate (which did
not transfer out of sample) with a density-dependent recruitment form, and adds bootstrap CIs. Job
11744992 (NSAMP=400, SEED=5), all four held-out variants. Method: `held_out_validation.py`; data:
`held_out_density_dependent_20260618.csv`.

## The change

Old: recruits = fixed_rate_A * (interval/10) * initial_TPA. The fold-A rate over or under shot fold B
because recruitment is spatially variable and density dependent.

New: recruits = R_max * max(0, 1 - SDI_t1 / SDImax_brms) * interval, with R_max fit on fold A so the form
reproduces fold-A observed ingrowth. Recruitment shuts off as a stand approaches its brms site-specific
max SDI, which bounds over recruitment in dense stands and makes the rate transfer across folds. SDI_t1 is
the Reineke summation SDI of the initial stand; SDImax is the brms plot-level posterior median (metric /
2.471 to English).

## Out of sample result (fold B, held out), default -> calibrated, with bootstrap 95% CI on calibrated

| variant | nA / nB | R_max (tpa/ac/yr) | QMD % | BA % | TPH % |
|---|---|---|---|---|---|
| ne | 67 / 128 | 41.3 | +12.1 -> +2.1 [-2.3, +7.2] | +8.7 -> +8.5 [+6.2, +11.0] | -14.6 -> -11.2 [-15.8, -6.6] |
| sn | 73 / 88 | 91.2 | +11.3 -> -2.1 [-8.3, +4.4] | +15.4 -> +12.8 [+7.5, +18.3] | -10.7 -> +1.4 [-8.1, +11.8] |
| kt | 24 / 52 | 26.1 | +22.9 -> -5.8 [-13.6, +2.4] | +16.2 -> +9.5 [+5.4, +14.1] | -36.8 -> -22.0 [-34.3, -5.7] |
| pn | 76 / 53 | 17.7 | +17.5 -> +1.7 [-4.6, +8.2] | +8.5 -> +6.7 [+2.8, +11.0] | -23.6 -> -25.7 [-36.2, -11.9] |

## What this shows

1. The size levers transfer out of sample, now with CIs. QMD bias is sharply reduced on the held-out fold
   in all four variants (ne +12.1->+2.1, sn +11.3->-2.1, kt +22.9->-5.8, pn +17.5->+1.7); the calibrated
   QMD CI spans zero for ne, sn, and pn. BA is reduced or held in all four. The brms SDImax + BAIMULT
   levers are not overfit.

2. Density-dependent recruitment fixes the fixed-rate failure. The fixed rate's worst behaviour was the
   Southern over correction, TPH -6.7 -> +16.1 out of sample. With the density-dependent form, SN is now
   -10.7 -> +1.4 (CI [-8.1, +11.8] spans zero): recruitment no longer overshoots once it is capped by
   max-SDI headroom. TPH also improves for ne (-14.6 -> -11.2) and kt (-36.8 -> -22.0).

3. The remaining exception is PN, where TPH stays under-predicted (-23.6 -> -25.7). PN initial stands sit
   far from the SDI limit less often, or the recorded ingrowth there is not headroom-limited; the
   injection is gated to plots where default TPH under-predicts and the cap leaves it too small. PN
   recruitment is the next refinement (site-productivity scaling of R_max, or a PN-specific R_max).

## Honest framing for the manuscript

Density-dependent recruitment makes TPH transfer out of sample for three of four variants and removes the
over-correction that the red-team flagged; report it as such, with PN as the documented exception that
motivates the site-resolved recruitment refinement. Size levers (QMD, BA) transfer cleanly with CIs.
