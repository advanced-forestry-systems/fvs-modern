# Calibrated NE vs Calibrated ACD on Acadian-footprint FIA — 2026-05-16

## Goal

With both NE and ACD Bayesian posteriors fit, score the same
Acadian-region (Maine, NH, Vermont; STATECD in {23, 33, 50}) FIA
remeasurement conditions through both calibrations and compare. This
is the apples-to-apples calibrated A/B that goes beyond the raw
runtime A/B in `../ne_vs_acd/`.

## Setup

Cardinal SLURM job 9718410 on branch acd-bridge-fix-2026-05-15.
Engine: calibration/R/19_fia_benchmark_engine.R with environment gate
FVS_ACD_RELABEL=TRUE enabled. The relabel logic at engine line 1490
retags NE-tagged FIA plots whose STATECD falls in the Acadian
footprint as ACD, so the variants loop projects them with ACD
posteriors while the rest of NE stays on NE posteriors.

Plot split after relabel: 13,693 ACD and 16,453 NE conditions
(n = 30,146 ME/NH/VT pairs total).

Posteriors used: from calibration/output/variants/{ne,acd}/. The NE
diameter-growth fit is ADVI (variational); the ACD fit is HMC with
rhat = 1.95 and converged = FALSE (a known limitation noted in the
parameter-level comparison report).

## Result: NE posteriors outperform ACD posteriors on the Acadian footprint

From fia_benchmark_ne_vs_acd_acadian.csv:

| Variant | mean_obs_BA (ft2/ac) | BA %RMSE | BA bias_pct | %RMSE reduction vs default |
| --- | --- | --- | --- | --- |
| ACD_acadian | 100.86 | 28.52 | +51.56 | -3.2 pp (calib worse) |
| NE_acadian  | 97.12  | 23.19 | +27.32 | -1.0 pp (calib worse) |

NE wins on both metrics: 5.3 pp lower %RMSE, about 24 pp less
positive bias. Both calibrations slightly degrade %RMSE compared to
their respective defaults on this footprint, driven by large
positive BA bias under calibration, consistent with the ACD
diameter-growth posterior wide sigma_b0 = 0.656 inflating DDS.

Verdict: For Acadian-region FIA stands, the NE calibrated posterior
is currently the better choice. The ACD posterior is distinct (zero
bit-identical parameters with NE) but not yet better fit. The fix
path is re-running the ACD HMC with more iterations and tighter
step-size control until rhat is under 1.05.

## Parameter-level posterior divergence

Diameter growth (Kuehne):

| Hyperparameter | NE p50 | ACD p50 |
| --- | --- | --- |
| mu_b0 | 0.897 | 0.399 |
| sigma_b0 | 0.147 | 0.656 |

Crown ratio:

| Fixed effect | NE p50 | ACD p50 |
| --- | --- | --- |
| b_Intercept | -0.411 | -0.635 |
| b_DBH_std | +0.054 | +0.127 |
| b_IDBH_stdE2 | +0.005 | -0.012 |

Mortality:

| Fixed effect | NE p50 | ACD p50 |
| --- | --- | --- |
| b_Intercept | +4.727 | +5.316 |
| b_DBH_std | -1.794 | -2.390 |

All 2,261 common posterior variables across the three submodels
differ between NE and ACD (zero bit-identical).

## Artifacts

- submit_benchmark_acd_relabel.sh — SLURM wrapper used to run job 9718410
- fia_benchmark_ne_vs_acd_acadian.csv — side-by-side benchmark output

The full per-variant rich CSV with all metrics is at
calibration/output/comparisons/manuscript_tables/fia_benchmark_results.csv.

## Open follow-ups

1. Re-fit ACD diameter-growth HMC with longer warmup and smaller
   step size until rhat is under 1.05 (currently 1.95). Likely needs
   adaptive reparameterization on the wide between-species spread.
2. Widen the ACD footprint beyond STATECD {23, 33, 50} to include
   New York Adirondacks via county-level filter (FIPS counties
   NY-019, NY-031, NY-033, NY-035, NY-041, NY-049, NY-089, NY-115).
   Currently a TODO in the engine docstring.
3. Repeat the A/B with the 81-row ACD post-pass calibration factors
   applied to see whether the stand-level corrections close the bias
   gap.
