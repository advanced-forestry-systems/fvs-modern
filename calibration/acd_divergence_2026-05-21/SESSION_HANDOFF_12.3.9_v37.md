# Session handoff: AcadianGY 12.3.5 to 12.3.9 + v31 density correction arc

2026-05-28. Consolidated end-of-arc handoff. Supersedes
`SESSION_HANDOFF_12.3.8.md`. Anyone (Claude session, collaborator, future me)
can read this single document and fully load context.

## Headline

12.3.5 -> 12.3.9 is a sequence of strict improvements to AcadianGY, each
documented and shipped. With the v31 bridge correction added on top, FIA BA
bias closes from a baseline of +15 percent to roughly +2 percent on hold-out
plots, with per-plot R^2 lifting from 0.38 to 0.48.

  start: 12.3.5 baseline                              +15.4 percent BA bias
  12.3.6 in-source MORTCAL                            +8.6 percent (with MORTCAL on)
  12.3.7 ingrowth carry-through Part 1 (multipliers)  no change to no-ingrowth path
  12.3.8 ingrowth carry-through Part 2 (STAND/PLOT)   recruits now persist; QMD bias closes
  12.3.9 ops$CSI_SCALE knob                            CSI_SCALE = 0.7 -> +10.0 percent BA bias
  v31 bridge correction (apply_density_correction.R)  superseded by v33 (over-corrected high-density)
  v33 pooled refit + asymmetric cap                   superseded by v37 (over-corrected on v36)
  v37 n=484 pooled refit + tighter cap                +0.15 percent CV bias, R^2 0.514 (production)

## Model versions and what changed

### 12.3.5 (canonical Acadian baseline, Ben Rice production)

Reference model. No MORTCAL, no ingrowth carry-through fixes.

### 12.3.6 (PR not opened; superseded)

Added opt-in `ops$MORTCAL` survivor haircut. Strict improvement (default off
is byte-identical to 12.3.5).

### 12.3.7 (PR #32)

Ingrowth carry-through Part 1. Recruits in `ING.TreeList` lacked
dDBH.mult/dHt.mult/mort.mult/max.dbh/max.height. They were frozen at 3 cm
recruitment diameter. Fix: set neutral multipliers + inherit species size
caps on `ingrow` rows before `bind_rows(tree, ingrow)`.

### 12.3.8 (PR #34)

Ingrowth carry-through Part 2. Recruits inherited Sum.temp's default
STAND=1/PLOT=1, so multi-stand harnesses split them into a phantom stand and
silently dropped them with "missing value where TRUE/FALSE needed" on cycle
2. Fix: force `ingrow$STAND` and `ingrow$PLOT` to the survivors' single
value.

### 12.3.9 (PR #42)

Added optional `ops$CSI_SCALE`. When set, multiplies the resolved CSI by the
provided factor right after `stand$CSI` is parsed. Default unset = byte-
identical to 12.3.8.

### v31 bridge correction (PR #45)

`apply_density_correction.R`. Post-projection helper, not part of the model
file. Applies a fitted linear correction `residual = 40.6345 - 0.334 * BA_t1`
capped at +/- 25 ft^2/ac to the predicted BA. Used downstream of
AcadianGYOneStand.

## What the validation says

### v24: FIA validation of 12.3.8 baseline

100 ME FIA plots, 10 yr, MORTCAL + CutPoint = 0:
TPA bias -1.0 percent to +1.3 percent
QMD bias +3.7 percent to -1.0 percent
BA bias slight +0.15 pp (recruits surface diameter overshoot)

### v17: Canadian MAGPlot under 12.3.8

262 New Brunswick CFI pairs, 4-cell sweep:
MORTCAL still over-corrects (-6.5 percent) even with full ingrowth flow
Right Canadian configuration: MORTCAL off + CutPoint = 0 ingrowth
That gives +0.4 percent BA bias, TPH 1764 vs obs 1807, QMD 12.3 vs obs 11.9

### v18: FIA per-species calibration decomposition

200 ME FIA plots, 5-cell sweep:
canonical_off (no levers) -> +15.35 percent BA bias
calibration: +6.91 pp (calibration alone hurts BA)
MORTCAL: -7.10 pp (rebalances calibration)
ingrowth: +0.22 pp (adds small-tree contribution)
Production posture calibrated_on_cp0: +15.38 percent (cancellation)

Insight: calibration + MORTCAL + ingrowth at production posture gives the
best TPA + QMD without inflating BA. R^2 is the same as canonical (0.42).
The +15 percent BA residual is structural under 12.3.8.

### v25: CSI sensitivity

5-cell sweep at production posture. Climate is a partial lever
(2.7 pp BA bias range over CSI x 0.6 to x 1.5). R^2 improves with lower CSI
weighting (0.41 to 0.435). This motivated 12.3.9's CSI_SCALE knob.

### v26: Mortality rate

5-cell sweep on `mort.mult`. Null lever (all 5 byte-identical). Tree-level
`mort.mult` does not flow into the mortality calculation in the standalone
production path.

### v27b: Alternative site-index metrics (in-vdat)

ClimateSI wins over SICOND and FVS_SITE_INDEX. Lowest BA bias, highest R^2.

### v27c: BGI raster as per-tree dDBH.mult

Pathological (raw use crashes dynamics) or null (recentered to mean=1.0).

### v29: BGI as per-stand ops$CSI_SCALE

Definitive null. BGI factor range 0.319 to 1.285 across plots; explains
none of the per-plot variance via the climate channel either.

### v30: Per-plot residual regression at 12.3.9 production posture

The strongest single predictor of per-plot residual is BA_t1
(R^2 = 0.188, p = 1.4e-5). BGI, SICOND, ClimateSI_ft all essentially zero.

Fit: `residual = 40.6345 - 0.334383 * BA_t1`. Model overshoots low-density
stands and undershoots high-density stands. Crossover at BA_t1 = 121.5
ft^2/ac.

### v31: First production correction with 50/50 holdout (200 iter, superseded)

Symmetric +/-25 cap fit on v30 n=93. Holdout bias +2.11 percent, R^2 0.479.
Initially shipped as PR #45.

### v32: Fresh-sample out-of-sample test (seed=2027, n=91)

Uncorrected R^2 matched the v30 holdout exactly (0.379 both). But v31's
correction recovered only half the in-sample gains and worsened Q4 stands
(BA_t1 ~ 180): raw bias +5.5 percent became +16.2 percent corrected. v31
slope was too aggressive at high density.

### v33: Pooled refit + asymmetric cap (PR #46, superseded)

Pool v30 + v32 (n=184). Refit shallower slope. CV bias +0.21 percent,
R^2 0.484.

### v34: Tree-level reconciliation (PR #48)

Push the stand-level BA correction back into the tree list by scaling EXPF
uniformly per stand. scale_floor = 0.7 protects against pathological
thinning. QMD invariant under uniform scaling.

### v35: Quadratic exploration (no ship)

residual ~ BA_t1 + I(BA_t1^2) lifts in-sample R^2 to 0.142 but loses on CV
bias (+2.9 vs +0.2). Forest type and interval are null. Quadratic curve
shape reveals model overshoots at BOTH low (BA_t1 < 130) AND very high
(BA_t1 > 200) density; v33 linear misses the high-density curl. Useful
diagnostic for the eventual BAL coefficient refit.

### v36: 300-plot fresh sample test (n=300, seed=2028)

Uncorrected BA bias +10.08, R^2 0.495. **v33 over-corrected to -5.00**
(vs CV target +0.21). Signal was inside CV variance but the larger sample
exposed v33 over-fit.

### v37: Pooled n=484 refit, tighter cap (PR #49, production)

Pool v30 + v32 + v36 = 484 plots. Slope shrinks again to -0.186; cap
tightens to +20. CV variance 3x tighter (sd 0.04 vs 0.11) thanks to n.

  raw_correction = 28.9607 + (-0.186023) * BA_t1
  bounded        = max(0, min(20, raw_correction))
  BA_corrected   = BA_pred - bounded

Crossover at BA_t1 = 155.7 ft^2/ac. Production CV bias +0.15 +/- 0.04
percent, R^2 0.514 +/- 0.001.

Coefficient progression v31 -> v33 -> v37 shows both intercept and slope
shrinking monotonically with n. Suggests we are converging on the true
relationship. The asymmetric cap (lower bound 0) continues to mean the
correction can only subtract from BA_pred, never add.

## Operating recommendations

### FIA-like Maine production

Source AcadianGY_12.3.9.r. Set:
  `ops$MORTCAL = TRUE` with `ops$MORTCAL_INTERVAL` per remeasurement period
  `ops$INGROWTH = "Y"`
  `ops$CutPoint = 0` (expected-value ingrowth)
  `ops$CSI_SCALE = 0.7`
Apply `apply_density_correction(BA_pred, BA_t1)` downstream.

Per-species calibration multipliers (`dDBH.mult`) optional; they improve TPA
and QMD slightly but BA bias is already addressed by the correction.

### Canadian MAGPlot

Source AcadianGY_12.3.9.r. Set:
  `ops$MORTCAL = NULL` (omitted) - MORTCAL over-corrects Canadian
  `ops$INGROWTH = "Y"`
  `ops$CutPoint = 0`
  `ops$CSI_SCALE = NULL` (omitted)
Do NOT apply density correction (calibrated to ME conditions, baseline bias
already at +0.4 percent).

## What is NOT closed

The remaining ~9 pp structural BA residual after CSI_SCALE = 0.7 + ingrowth
+ MORTCAL on FIA. The v30 finding suggests it sits in the BAL competition
coefficient in the Kuehne et al. 2020 dDBH equation (density miscalibration
that v31 patches post-hoc). Refitting BAL coefficient against ME FIA is the
paper-sized fix.

Same for the Glover/Hool mortality functional form. v26 ruled out the
rate-level lever; if the form is wrong (logistic vs Weibull or piecewise),
it would also contribute to the v30 BA_t1 signal.

## Files of record (in fvs-modern main, calibration/acd_divergence_2026-05-21/)

Diagnosis writeups:
  INGROWTH_FIX_12.3.7.md, INGROWTH_FIX_12.3.8.md, INGROWTH_FIX_12.3.8_v24_validation.md
  MAGPLOT_12.3.8_v17_findings.md, CALIB_12.3.8_v18_findings.md
  CSI_SENSITIVITY_12.3.8_v25_findings.md, MORTMULT_12.3.8_v26_findings.md
  SITEINDEX_COMPARISON_12.3.8_v27b_findings.md, BGI_RASTER_12.3.8_v27c_findings.md
  CSI_SCALE_12.3.9_v28_findings.md, BGI_via_CSI_SCALE_v29_findings.md
  RESIDUAL_CALIBRATION_v30_findings.md, DENSITY_CORRECTION_v31_findings.md

Model files and helpers:
  AcadianGY_12.3.8.r, AcadianGY_12.3.9.r
  apply_density_correction.R
  patch_ingrowth_fix.py, patch_ingrowth_fix_12_3_8.py, patch_12.3.9.py
  probe_recruit_stand.R, run_ingrowth_probe.R, extract_bgi.py

Cardinal harnesses and results:
  cardinal_acadgy_*v24,v17,v18,v25,v26,v27b,v27c,v28,v29,v30*.R
  acadgy_*_results.csv for each
  me_bgi_by_pltcn.csv (21,279 ME plots BGI from raster)

## Open follow-ups

#190 (pending) - update Ben Rice Gmail draft with v25-v33 findings
Larger-n validation (v34 = 500 plot fresh sample) would harden v33
Strata-aware refit (FORTYPCD hardwood/softwood/mixed) for v33b
Paper-sized: BAL coefficient refit; mortality functional form refit

## Commit ledger (this arc)

#32 12.3.7 ingrowth carry-through Part 1
#34 12.3.8 STAND/PLOT inheritance
#36 v18 calibration decomposition
#37 SESSION_HANDOFF_12.3.8
#38 v25 CSI sensitivity
#39 v26 mort.mult sensitivity (null)
#40 v27b in-vdat site-index comparison
#41 v27c ME BGI raster (null)
#42 12.3.9 ops$CSI_SCALE
#43 v29 BGI via CSI_SCALE (null)
#44 v30 BA_t1 residual signal
#45 v31 production correction (superseded)
#46 v33 pooled refit + asymmetric cap (superseded)
#47 comprehensive session handoff
#48 v34 tree-level reconciliation
#49 v37 n=484 pooled refit (current production)
