# Master results table (red-team pass)
2026-06-18. Item 4 of the integration roadmap. One reconciled set of results with bootstrap CIs, plus the
red-team's outstanding items: brms match rate per variant, the Lake States SDImax note, the volume
definition, and the standing caveats. Reframes the work as a disturbance-aware benchmark plus a prototype
adjustment layer, out-of-sample where claimed.

## 1. The disturbance artifact (established, three ways)

| evidence | result |
|---|---|
| COND stratification, all 22 FIA variants | undisturbed median BA bias +1.8%; harvested +42%; pooled +14% |
| removal-simulation converse test | harvested BA bias collapses to undisturbed when recorded harvest is simulated: ne +51.2->+1.3, sn +110.1->+0.1, pn +41.9->+12.9 |
| fvs-conus projector, same conditions | same pattern; harvest inflation appears in the projector too |

The pooled over-prediction is unsimulated removal, not a growth-equation bias. Validate on undisturbed
plots or simulate the recorded harvest before comparing.

## 2. Engine arms A (default) and B (keyword-calibrated), out-of-sample, with CIs

Spatial held-out fold B (county-hash folds; brms maxSDI + density-dependent recruitment + fold-A BAIMULT
derived on fold A only). FVS engine, COND-undisturbed, one remeasurement cycle. Job 11745221, NSAMP=400.
Values are % bias, default -> calibrated. (ut, ec, wc completing; refresh on job close.)

| variant | nB | BA | TPH | QMD | merch VOL |
|---|---|---|---|---|---|
| ne | 128 | +8.7 -> +8.5 | -14.6 -> -11.2 | +12.1 -> +2.1 | +21.2 -> +21.2 |
| sn | 88 | +15.4 -> +12.8 | -10.7 -> +1.4 | +11.3 -> -2.1 | +9.4 -> +3.5 |
| kt | 51 | +16.6 -> +9.7 | -36.9 -> -22.6 | +23.0 -> -5.6 | +17.0 -> +8.0 |
| pn | 53 | +8.5 -> +6.7 | -23.6 -> -25.7 | +17.5 -> +1.7 | +9.1 -> +8.6 |
| nc | 42 | +3.1 -> +6.1 | -14.5 -> -4.9 | +21.9 -> -0.5 | -9.9 -> -5.6 |
| cr | 15 | +28.0 -> +23.9 | +18.0 -> +27.1 | -1.0 -> -15.4 | +37.1 -> +33.4 |

Reading: QMD bias is sharply reduced out-of-sample in every variant with adequate sample (ne, sn, kt, pn,
nc all collapse toward zero); BA and merch volume improve where BAIMULT engages (sn, kt). cr (n=15) is the
small-sample exception where the size lever over-corrects QMD and TPH; treat western point estimates with
small n as indicative, not final. Full per-arm bootstrap CIs in fourarm_engine_20260618.csv.

## 3. Projector arms A' (default) and C (fvs-conus equations), NE, with CIs

fvs-conus standalone projector, 21,811 COND-undisturbed NE conditions. Values % bias, default -> fvs-conus.

| metric | A' default | C fvs-conus | 95% CI (C) |
|---|---|---|---|
| BA | -12.3 | -7.2 | [-7.3, -7.0] |
| TPA | -7.2 | -6.2 | [-6.5, -6.0] |
| QMD | -5.3 | -3.1 | [-3.2, -3.0] |
| merch VOL (CFNET) | -15.7 | -9.2 | [-9.4, -9.0] |

The engine over-predicts and the projector under-predicts undisturbed (opposite signs, different
machinery), so arms are compared as within-framework improvement, not absolute bias. The fvs-conus
equations reduce BA and volume bias most; the keyword calibration reduces QMD and TPH most. Complementary,
which is the arm-D hypothesis (see 20260618_fourarm_result.md).

## 4. brms SDImax match rate per variant (red-team item)

Fraction of FIA plots in each variant's states carrying a brms site-specific max SDI, with the SDImax
distribution (English, metric / 2.471). Job 11745289.

| variant | FIA plots | brms match | match % | SDImax median [p10-p90] |
|---|---|---|---|---|
| acd | 33,643 | 24,684 | 73.4 | 386 [251-568] |
| ne | 66,732 | 42,156 | 63.2 | 370 [229-543] |
| ec/wc/pn | 64,293 | 38,317 | 59.6 | 371 [180-685] |
| nc | 85,008 | 37,428 | 44.0 | 347 [165-655] |
| sn | 196,282 | 80,671 | 41.1 | 337 [209-495] |
| ca | 43,814 | 12,877 | 29.4 | 347 [158-676] |
| ls | 341,358 | 97,067 | 28.4 | 322 [191-491] |
| ci | 26,876 | 7,580 | 28.2 | 314 [173-548] |
| ie/kt | 73,719 | 16,694 | 22.6 | 307 [167-525] |
| cr | 62,919 | 13,874 | 22.1 | 317 [175-538] |
| cs | 144,607 | 28,763 | 19.9 | 303 [196-446] |
| ut | 67,765 | 12,556 | 18.5 | 265 [147-471] |

Caveat: states are mapped to variants at the state level, so variants sharing states (ie/kt; ec/wc/pn)
report identical rows; a true per-variant rate needs the FVS variant polygon. Match rate is brms coverage
of all FIA plots; the calibration falls back to FIA SICOND / default variant SDIMAX where no brms estimate
exists.

## 5. Lake States SDImax note (red-team item)

Lake States brms coverage is 28.4% with SDImax median 322 [191-491], which sits inside the cross-variant
range (cs 303, sn 337, ne 370), so the LS SDImax values themselves are not out of family. The flagged LS
anomaly therefore is not a wild SDImax value; it should be pinned to the LS-specific calibrated benchmark
(allvar_calibration.csv shows LS among the variants where calibrated RMSE barely moves). Flagged for the
LS calibration rerun; not resolved here, and not papered over.

## 6. Volume definition (red-team item)

Stand merch volume is FVS MCuFt (merchantable cubic feet, summary column) compared to FIA VOLCFNET summed
per plot as VOLCFNET * TPA_UNADJ, both scaled to m3/ha by VOLc = 0.0699055. This is NET merchantable cubic
volume. Gross (VOL_CFGRS) and board-foot (VOL_BFNET) are available in the projector output for a
sensitivity pass; the engine arm currently reports net merch only. Sensitivity to the gross/net and
merch-top definitions is the open volume item.

## 7. Outstanding (carried to next pass)

1. Close ut, ec, wc in the engine four-arm; recompute cross-variant medians on all nine.
2. Extend removal-sim and held-out to more variants and larger n (western undisturbed samples are thin).
3. Volume gross/net + merch-top sensitivity; biomass via FFE.
4. Arm D (combined) once the fvs-conus equations run in the engine (fvs2py) or via per-species multiplier
   emulation, to test the complementarity hypothesis directly.

Sources: fourarm_engine_20260618.csv, fourarm_projector_NE_20260618.csv, brms_match_rate_20260618.csv,
held_out_density_dependent_20260618.csv, 20260617 disturbance benchmark and removal-sim.
