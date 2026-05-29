# Density-dependent bridge correction v33: pooled refit with asymmetric cap

2026-05-28. **Supersedes v31.** v32 fresh-sample test (n=91, seed=2027)
revealed that the v31 coefficients (fit on n=93, seed=42) over-correct
high-density stands. The v32 BA_t1 signal was R^2=0.057 (vs 0.188 in v30) and
the v31 +/-25 correction WORSENED bias from +5.5 to +16.2 percent in
Q4 (BA_t1 mean 180 ft^2/ac) where the model was already well-behaved.

v33 pools both samples (n=184), refits the slope, and applies an asymmetric
cap so the correction can only subtract from BA_pred, never add.

## v33 production formula

    raw_correction = 36.9549 + (-0.235987) * BA_t1
    bounded        = max(0, min(25, raw_correction))   # lower 0, upper 25
    BA_corrected   = BA_pred - bounded

Coefficients fit on the pooled v30 + v32 sample (n=184 ME FIA, 10-year
remeasurement, 12.3.9 production posture with MORTCAL on, CutPoint = 0,
CSI_SCALE = 0.7). Crossover at BA_t1 = 156.6 ft^2/ac. Beyond that, the
correction is bounded to 0 and the raw BA_pred is preserved.

## 5-fold CV results (pooled n=184, 50 random shuffles)

| config              | CV bias        | CV R^2         |
|---------------------|----------------|----------------|
| Uncorrected         | +14.46 percent | 0.380          |
| Symmetric +/-25     | +1.38 +/- 0.12 | 0.472 +/- 0.005 |
| **Asymmetric (0, +25)** | **+0.21 +/- 0.11** | **0.484 +/- 0.003** |
| Asymmetric (0, +20) | +1.65 +/- 0.08 | 0.472 +/- 0.002 |
| Asymmetric (0, +15) | +3.77 +/- 0.08 | 0.457 +/- 0.001 |

The asymmetric (0, +25) configuration wins on both metrics. Tiny CV variance
(SD ~0.003 for R^2) means the result is extremely stable across folds.

## v32 fresh-sample diagnosis (why v31 failed)

| sample | uncorrected bias | v31 +/-25 bias | uncorrected R^2 | v31 R^2 |
|--------|------------------|----------------|------------------|---------|
| v30 (holdout) | +10.86 percent  | +2.11 percent  | 0.379           | 0.479   |
| **v32 (fresh)** | **+17.92** | **+12.95**    | **0.379**       | **0.427** |

Uncorrected R^2 is essentially identical (0.379 vs 0.379) — the model behaves
consistently. But the correction recovered only half the in-sample gains on
fresh data.

Per-quartile decomposition on v32 explained why:

| quartile | BA_t1 mean | v32 raw bias | v32 v31 corrected |
|----------|------------|--------------|---------------------|
| Q1       | 32.8       | +131 percent | +51 percent         |
| Q2       | 84.6       | +9.3 percent | -4.7 percent        |
| Q3       | 122.5      | +11.4        | +11.7               |
| Q4       | 180.0      | **+5.5**     | **+16.2** (worse)   |

v31's slope was too aggressive at high density. With asymmetric cap (0, +25),
Q4 keeps the raw BA_pred (correction bounded to 0), so the high-density
stands are not damaged.

## Comparing v31 and v33 on pooled n=184 (in-sample)

| config | bias | R^2 |
|--------|------|-----|
| Uncorrected           | +14.46 | 0.380 |
| v31 (+/-25)           | +7.25  | 0.462 |
| v31 (asym 0, +25)     | +3.24  | 0.494 |
| v33 (sym +/-25)       | +1.30  | 0.480 |
| **v33 (asym 0, +25)** | **+0.21** | **0.487** |

Both improvements (pooled refit and asymmetric cap) help; the combination is
the production ship.

## Updated bridge helper

`apply_density_correction.R` updated with v33 coefficients. API change:

- `cap` parameter replaced with `upper_cap` (default 25) and `lower_cap`
  (default 0). The defaults give the asymmetric production behavior.
- `ACD_DENSITY_CORRECTION` list now records both caps, CV stats, and notes
  that it supersedes v31.

Smoke test on Cardinal R 4.4.0 verified:
- BA_t1 = 10 (sparse): raw correction +34.6 -> upper cap +25, BA_pred reduced
- BA_t1 = 156.6 (crossover): raw 0, no correction
- BA_t1 >= 180 (dense): raw correction negative -> lower bound 0, BA_pred preserved

## Operating recommendation

Use v33 default behavior `apply_density_correction(BA_pred, BA_t1)`. The
correction now safely:
- closes the BA overshoot at low and mid stand densities,
- leaves high-density predictions alone where the raw model is already close.

Still: do NOT apply to Canadian MAGPlot. The v17 finding stands - Canadian
+0.4 percent baseline bias has no residual to close, and the ME-tuned slope
would damage it.

## Files

  apply_density_correction.R  (updated with v33 coefficients)
  cardinal_acadgy_v31test_v32.R  (v32 harness)
  acadgy_v31test_v32_results.csv (v32 aggregate)
  acdgy_v31test_v32_perplot.csv  (v32 per-plot for n=91)
  DENSITY_CORRECTION_v33_findings.md  (this memo)
  DENSITY_CORRECTION_v31_findings.md  (kept as superseded record)

## What's next

The 14-pp bias closure with asymmetric (0, +25) is now well-validated at
n=184. Three follow-up moves that would harden it further:

1. **Larger n.** Run v34 = 500 ME FIA plots, fresh seed, apply v33 without
   refit. Should match the 50-shuffle CV result (bias ~0.2 percent, R^2
   ~0.484) within +/- 0.5 percent and +/- 0.02 R^2.
2. **Strata-aware refit.** Fit separate slopes for hardwood / softwood /
   mixed forest types (FORTYPCD). The pooled slope may be a weighted
   average of three different per-type slopes.
3. **BAL as a covariate alternative.** v30 already tried BA_t1 + ClimateSI
   (didn't help). BAL (basal area in larger trees) is the actual variable
   in the Kuehne dDBH equation. If BAL data is in vdat, residual ~ BAL +
   BA_t1 might explain more.

For now, v33 is the production ship.
