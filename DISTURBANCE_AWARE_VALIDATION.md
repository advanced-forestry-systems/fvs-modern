# Disturbance-aware validation and the four-lever calibration (fvs-modern)
2026-06-17. Companion to CALIBRATION.md. Defines the validation protocol and the deployable adjustment
levers established from the disturbance-stratified FIA benchmark.

## Why this exists

The Bayesian component pipeline (CALIBRATION.md) refits FVS equations. This document covers the
complementary question: how to VALIDATE FVS against FIA correctly, and which adjustment levers actually
move the management-relevant outputs (basal area, density, size, volume, biomass) and at what time scale.

## Validation protocol (mandatory)

Stratify every FIA remeasurement plot by FIA COND history before comparing to an undisturbed FVS run:
- harvested: any TRTCD in {10} (cutting)
- disturbed: any DSTRBCD > 0
- undisturbed: neither
Compare on the undisturbed stratum (the stands the default projection represents), or simulate the
observed removals before comparing. Pooled comparisons that ignore disturbance manufacture a +14-28%
"over-prediction" that is mostly unsimulated harvest (+33-55% on harvested plots). Metrics: BA, TPH, QMD,
and merch volume with matched definitions (FVS MCuFt vs FIA VOLCFNET).

## What is biased on undisturbed plots (19 variants)

| metric | median bias | direction |
|---|---|---|
| basal area | +7% | mild over |
| QMD | +10% (18/19 over) | trees too large |
| TPH | -14% (15/19 under) | too few stems |
| merch volume | +15-24% (East/LS/PNW) | driven by QMD |

Top height is within +/-3% and hides the HT-DBH curve error (below).

## The four deployable levers (by time scale)

1. Maximum SDI -- long-term (multi-decade) density/BA ceiling. Use the plot-level brms site-specific
   SDImax (`brms_SDImax`), converted metric/2.471 -> English, emitted as SDIMAX per species. Replacing
   the FVS default cuts 100-yr BA 5-37% where the default exceeds the FIA self-thinning limit (worst in
   western variants). This is the dominant long-term control.
2. Recruitment / ingrowth -- decadal stem density. FVS adds no background ingrowth in undisturbed runs
   (establishment is disturbance-triggered). Inject a recruitment cohort at the per-variant observed FIA
   rate (17-70%/decade), SIGN-AWARE: only where the variant under-predicts TPH (most). Closes the TPH gap
   and lowers QMD.
3. Diameter growth -- standing size and therefore standing volume. Signed BAIMULT (~0.90 default; size to
   each variant's QMD/volume over-prediction).
4. HT-DBH curve -- height, on which volume/biomass depend. Variant- and size-specific bias (pooled -7% at
   1-3 in to +9% at 19-40 in; IE/KT too tall, SN under-large, PN wrong shape). In DBH-only inventories
   this inflates volume +5-10%. Apply the per-species correction ratios (`htdbh_recalibration_ratios`),
   or refit the curve coefficients where the error is shape (PN, SN). Supply measured heights when
   available (FVS self-calibrates the curve). REGHMULT scales height GROWTH, not the static curve.

## Deployable config

`calibration_config.csv` lists per variant: SDImax source (brms_plot_level), ingrowth rate (%/decade),
inject flag (sign-aware), BAIMULT, and HT-DBH correction. ACD and ADK are NE sub-variants (customR); the
layer applies unchanged because brms SDImax is keyed by plot.

## Validation result (stress test)

Default vs fully-calibrated across all 15 variants on the COND-undisturbed stratum, median |bias|:
QMD 11.1 -> 4.4%, BA 9.9 -> 7.9%, merch volume 17.6 -> 14.3%, TPH 15.8 -> 11.8% (TPH improves further
with sign-aware injection, applied only where the variant under-predicts TPH; the full four-metric
stress table is `calib_final.csv`). QMD, the driver of the volume over-prediction, is more than halved.

## Scripts (calibration/validation/)

multivar_v2.py (COND benchmark), ingrowth_decomp.py (recruitment), maxsdi_fia.py + maxsdi_longterm.py
(SDImax), seed_test.py (ingrowth injection), htdbh_assess.py + vol_height_sens.py (HT-DBH), calib_final.py
(full stress test), plus render_*.py figures. brms_SDImax.csv is the site-specific max-SDI source.
