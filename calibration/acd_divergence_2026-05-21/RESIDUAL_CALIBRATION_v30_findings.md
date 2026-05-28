# Per-plot residual analysis on 12.3.9 (CSI_SCALE = 0.7): density-dependent calibration lever

2026-05-28. The v27c, v29 input-based tests for BGI were null. v30 asks the
orthogonal question via post-hoc residual regression: does any per-plot
covariate predict the model's per-plot BA error well enough to support a
fitted calibration factor?

## Result

100 ME FIA plots, 12.3.9 production posture (MORTCAL on, CutPoint = 0,
CSI_SCALE = 0.7). 93 plots had complete covariate coverage. Univariate
regressions on per-plot residual = BA_pred - BA_obs:

| predictor                | R^2     | adj R^2 | p-value |
|--------------------------|---------|---------|---------|
| BGI                      | 0.0006  | -0.011  | 0.82    |
| log(BGI)                 | 0.0014  | -0.010  | 0.73    |
| SICOND                   | 0.0015  | -0.010  | 0.72    |
| ClimateSI_ft             | 0.0077  | -0.003  | 0.40    |
| **BA_t1**                | **0.1877** | **0.1788** | **1.4e-5** |
| BA_t1 + I(BA_t1^2)       | 0.2390  | 0.2220  | 4.6e-6  |
| BA_t1 + ClimateSI_ft     | 0.1918  | 0.1739  | 6.9e-5  |
| BA_t1 + BGI              | 0.1891  | 0.1710  | 8.0e-5  |

The clean finding: **BA_t1 (initial stand basal area) is the strongest single
predictor of the model's per-plot residual.** It explains roughly 19 percent
of residual variance with p = 1.4e-5. Adding a quadratic term lifts to 24
percent. BGI, SICOND, ClimateSI_ft add essentially nothing beyond BA_t1.

## Fitted correction

    residual = 40.6345 - 0.3344 * BA_t1     (R^2 = 0.188, n = 93)

The model overshoots low-density stands and undershoots high-density stands,
with the crossover near BA_t1 = 121 ft^2/ac.

| starting BA (ft^2/ac) | predicted overshoot |
|-----------------------|----------------------|
| 50                    | +23.9                |
| 90                    | +10.5                |
| 120                   | +0.5                 |
| 180                   | -19.6                |

## Effect of applying the correction (in-sample)

| metric                | original  | corrected |
|-----------------------|-----------|-----------|
| BA_pred mean          | 107.26    | 96.60     |
| BA_obs mean           | 96.60     | 96.60     |
| BA bias %             | +11.04    | 0.00      |
| R^2 (BA_pred vs obs)  | 0.382     | **0.533** |

BA bias closes to exactly zero (by construction) and per-plot R^2 lifts by
0.15. The lever the entire BGI investigation was looking for is here, just
not in BGI: it's in stand density.

## Interpretation

This is the signature of a density-dependent miscalibration in the dDBH or
mortality equation. The Kuehne et al. 2020 diameter increment uses a BAL
(basal area in larger trees) competition term. If the BAL coefficient is too
weak, the model under-suppresses diameter growth in dense stands, which
shows as overshoot in low-density stands (where BAL is small) and undershoot
in high-density stands (where BAL is large) — exactly the pattern we see.

Alternative explanation: the Glover/Hool mortality is missing density-dependent
acceleration in dense stands, so high-density stands keep too many trees and
their BA undershoots because the model is removing the wrong ones.

Both stories are plausible; both are paper-sized refits. The post-hoc
correction is the practical fix that ships now.

## Why BGI was null even though density isn't

BGI is fundamentally a remote-sensing measure of realized growth: a stand
with high BGI grew faster than the model expected based on what was
visible from the air. The model already sees stand density through BA_t1
explicitly (and through BAL inside the dDBH equation). So BGI is mostly
proxying what the model already knows about density. The residual it
could explain is what is left after density is accounted for, and
within-stand productivity variation conditional on density is small on
this sample.

## What to ship

Three nested options for the customRun bridge in 12.3.9 (no model code
changes; just a post-projection multiplier on BA_pred):

  1. **Linear**: BA_corrected = BA_pred - (40.6 - 0.334 * BA_t1).
     In-sample bias 0 percent, R^2 0.533. Risk: only validated on n = 93,
     could over-fit.

  2. **Capped linear**: same formula but clamp the correction to within
     +/- 25 ft^2/ac so it can't introduce a worse error than the original
     baseline residual sd of 39 ft^2/ac.

  3. **Hold off and refit BAL coefficient**. Re-fit the Kuehne BAL term
     against ME FIA. This is the paper-sized version; closes the bias at
     the equation level instead of post-hoc.

The right autopilot move is option 2 with holdout validation. v31 will
implement it as a bridge function `apply_density_correction(BA_pred, BA_t1)`
and validate on a 50/50 holdout split of the v30 sample.

## Files

`cardinal_acadgy_residualcal_v30.R`, `acdgy_residualcal_v30_perplot.csv`
(per-plot data), `acdgy_residualcal_v30_results.csv` (aggregate). Cardinal
SLURM job 10968650.
