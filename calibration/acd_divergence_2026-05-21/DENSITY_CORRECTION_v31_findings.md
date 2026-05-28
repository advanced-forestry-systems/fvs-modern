# Density-dependent bridge correction (v31): capped linear, holdout validated

2026-05-28. v30 identified `BA_t1` as the strongest single predictor of the
per-plot residual on 12.3.9 production posture (R^2 = 0.188, p = 1.4e-5). v31
implements the post-projection correction as a ready-to-ship bridge helper
and validates it with 200 random 50/50 holdout splits on the n = 93 sample.

## Production formula

    raw_residual = 40.6345 + (-0.334383) * BA_t1
    capped       = max(-25, min(25, raw_residual))
    BA_corrected = BA_pred - capped

Where BA_pred and BA_t1 are in ft^2/ac. The cap defaults to +/- 25 ft^2/ac
(approximately +/- 26 percent of the mean observed BA), preventing the
correction from introducing larger errors than the original residual standard
deviation of 39 ft^2/ac on edge cases.

## Holdout validation (200 iterations, 50/50 splits)

| config              | TEST raw bias  | TEST corr bias | TEST raw R^2  | TEST corr R^2 |
|---------------------|----------------|----------------|---------------|---------------|
| Linear, no cap      | +10.86 +- 4.45 | **+0.41** +- 7.10 | 0.379 +- 0.146 | 0.481 +- 0.116 |
| **Linear, +/-25**   | **+10.86 +- 4.45** | **+2.11** +- 6.12 | **0.379 +- 0.146** | **0.479 +- 0.115** |
| Linear, +/-20       | +10.86 +- 4.45 | +3.14 +- 5.75 | 0.379 +- 0.146 | 0.471 +- 0.117 |
| Linear, +/-15       | +10.86 +- 4.45 | +4.53 +- 5.35 | 0.379 +- 0.146 | 0.458 +- 0.121 |
| Quadratic, +/-25    | +10.86 +- 4.45 | +3.85 +- 5.90 | 0.379 +- 0.146 | **0.489** +- 0.110 |

Three things to note:

1. **The correction generalizes.** Mean test-set bias drops from +10.86
   percent to +2.11 percent at the +/- 25 cap, with sd 6.1 percent across
   splits. The result is robust, not an artifact of any particular split.

2. **The cap trades small bias closure for variance reduction.** Linear with
   no cap gets test bias to +0.41 percent but only 0.002 better R^2 than
   +/- 25. The cap insurance against bad-actor predictions is essentially
   free.

3. **Quadratic gives slightly higher R^2 (0.489 vs 0.479) but worse mean
   bias (+3.85 vs +2.11).** Linear with +/- 25 cap is the right operating
   choice.

## Bridge helper

`apply_density_correction.R` provides:

  `apply_density_correction(BA_pred, BA_t1, cap = 25)` -> corrected vector
  `apply_density_correction_verbose(BA_pred, BA_t1)` -> diagnostic data.frame
  `ACD_DENSITY_CORRECTION` -> named list of coefficients and provenance

Source the file once at bridge startup; call the function after the
AcadianGYOneStand projection returns. No model code changes required.

Smoke test verified on Cardinal R 4.4.0: cap activates at both BA_t1 extremes
(at 20 ft^2/ac the raw residual is +33.9 but capped to +25; at 220 ft^2/ac
the raw residual is -32.9 but capped to -25). Crossover at BA_t1 = 121.5
produces correction ~ 0 as expected.

## What this is and is not

This is a **post-hoc multiplicative correction**, not a fix to the underlying
equation. Two important caveats:

(a) **Calibration sample is small.** n = 93 ME FIA plots, 10-year
    remeasurement, AcadianGY 12.3.9. Should NOT be applied to:
    - Canadian MAGPlot conditions (v17 showed +0.4 percent baseline; no
      residual to close)
    - Non-ME FIA conditions without revalidation
    - Projection intervals materially different from 5 to 10 years

(b) **The structural problem still exists.** The density-dependent bias is
    the signature of a BAL coefficient miscalibration in the Kuehne et al.
    2020 dDBH equation, or missing density-acceleration in Glover/Hool
    mortality. The right long-term fix is to refit those coefficients
    against ME FIA. The bridge correction buys us a +0.1 R^2 lift and a
    near-zero mean bias for the production pipeline today; the paper-sized
    work is still pending.

## Files

  apply_density_correction.R                              the bridge helper
  acdgy_residualcal_v30_perplot.csv                       the 93-plot calibration data
  DENSITY_CORRECTION_v31_findings.md                      this memo
  RESIDUAL_CALIBRATION_v30_findings.md                    upstream v30 analysis
