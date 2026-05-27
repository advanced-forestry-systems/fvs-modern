# Climate Site Index sensitivity under AcadianGY 12.3.8: CSI is a partial lever on BA bias

2026-05-27. Probes the structural +11 percent FIA BA residual under 12.3.8 by
scanning CSI (Climate Site Index) at five scale factors holding all other
parameters at the production posture (MORTCAL on, CutPoint = 0 expected-value
ingrowth, no per-species diameter calibration). 100 ME FIA plots, 10 annual
cycles, same sample as v24.

## Result

| CSI scale | BA bias % | TPA  | QMD cm | R^2  |
|-----------|-----------|------|--------|------|
| 0.6       | +9.57     | 1033 | 4.94   | 0.435 |
| 0.8       | +10.38    | 1038 | 4.93   | 0.430 |
| **1.0**   | **+11.05** | **1043** | **4.92** | **0.423** |
| 1.2       | +11.63    | 1048 | 4.91   | 0.418 |
| 1.5       | +12.31    | 1055 | 4.89   | 0.410 |

Observed: BA 94.72 ft^2/ac, TPA 1029, QMD 4.97 in.

The CSI*1.0 row reproduces v24 (BA +11.05 percent, TPA 1043, QMD 4.92) to
four decimals. Sanity check passes.

## Interpretation

Climate sensitivity in the Kuehne et al. 2020 dDBH equations is a real lever
on stand BA, but a partial one. Across the 0.6 to 1.5 scaling range (factor
2.5x in CSI), BA bias moves by 2.74 pp (9.57 to 12.31). The elasticity is
roughly 0.27 pp BA bias per 0.1 CSI scale step in this regime.

Three observations make this more than a curiosity:

(a) The response is monotonic, smooth, and the expected direction. Higher
    CSI gives higher diameter growth gives higher BA. The model is
    climate-sensitive in the predicted direction.

(b) R^2 improves systematically with lower CSI weighting, from 0.41 at
    1.5x to 0.435 at 0.6x. This is not just a mean shift; the per-plot
    fit gets better when CSI's influence is dampened. That is a
    calibration signal, not noise.

(c) CSI*0.6 with MORTCAL on and ingrowth on gives TPA 1033 (obs 1029, off
    by 0.4 percent) and QMD 4.94 (obs 4.97, off by 0.6 percent) - the
    closest the model has come to ground truth on TPA and QMD in any
    configuration tested so far. BA is still +9.6 percent over.

## What this means for the residual

Of the +11 percent BA residual under 12.3.8 production posture, climate
weighting accounts for roughly 1.5 pp (the gap between CSI*0.6 and CSI*1.0).
The remaining ~9.6 percent sits in:

  1. The intrinsic Kuehne et al. dDBH coefficients (the base-rate diameter
     growth, independent of CSI).
  2. The Glover/Hool mortality functional form (logistic in size class) -
     a Weibull or piecewise survival function may capture late-stand
     acceleration that the logistic underweights.
  3. A structural variant-level ceiling specific to the Acadian variant
     on Maine FIA conditions.

The next sensitivity scan to try is mortality functional form: re-run v25
with a Weibull survival function and measure the BA response. That is a
larger commitment than this scan because it requires re-fitting mortality
coefficients, not just scaling an input.

## Action items

Two options consistent with this finding:

**Production tweak.** Set CSI*0.7 as the default scaling factor in the
customRun bridge. Drops BA bias from +11.1 to roughly +10 percent and
improves R^2 by ~0.01 with no other code changes. Low-risk, partial closure.

**Bigger commitment.** Refit the climate term in the Kuehne et al. dDBH
equations against ME FIA data directly. That recovers the full lever but
opens a calibration question (which species, which functional form for
CSI) that is a paper-sized commitment.

The right autopilot move is the production tweak as a 12.3.9 release if Ben
is interested, with the refit reserved for the next research session.

## Files

`cardinal_acadgy_csisensitivity_v25.R` and `acadgy_csisensitivity_v25_results.csv`
in `acd_divergence_2026-05-21/`. Cardinal run as SLURM job 10641397
(approximately 20 minutes on c0039).
