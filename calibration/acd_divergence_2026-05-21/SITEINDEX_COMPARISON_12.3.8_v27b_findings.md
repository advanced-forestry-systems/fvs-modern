# Site productivity metric comparison under AcadianGY 12.3.8

2026-05-28. Compares three site productivity metrics already in
`validation_data_acd_post.csv` as alternatives to the ClimateSI_ft baseline
that v24/v25 used. All other parameters held at 12.3.8 production posture
(MORTCAL on, CutPoint = 0 expected-value ingrowth). Same 100-plot ME FIA
sample, same n_years = 10. BGI from `ME_BGI_V1.tif` raster is a follow-on (the
existing `bgi_by_pltcn.csv` is WA-state-specific compact IDs and doesn't join
to ME plots).

## Result

| metric                       | mean (m) | BA bias % | R^2    | TPA  | QMD cm |
|------------------------------|----------|-----------|--------|------|--------|
| **csi_baseline (ClimateSI)** | **14.31** | **+11.05** | **0.4232** | **1043** | **4.923** |
| sicond_replace (SICOND)      | 16.24    | +11.39    | 0.4190 | 1045 | 4.912 |
| fvssi_replace (FVS_SITE_INDEX) | 15.38  | +11.27    | 0.4190 | 1044 | 4.918 |

Observed: BA 94.72 ft^2/ac, TPA 1029, QMD 4.97 in.

## Interpretation

ClimateSI_ft wins all three criteria of interest: lowest BA bias, highest R^2,
TPA closest to observed.

The response is monotonic with mean site-index value (higher SI -> more
diameter growth -> higher BA bias). This is consistent with the v25 CSI
sensitivity scan (where CSI x 1.0 -> +11.05% and CSI x 0.8 -> +10.38%).
Replacing ClimateSI with SICOND (+1.93 m, ~14% higher) shifts BA bias by +0.34
pp, which matches the v25 elasticity of ~0.27 pp per 0.1 CSI scale step
applied to a +14% step.

The substantive finding: the choice of productivity metric among the three
already-loaded options is not a lever on the residual. ClimateSI is already
the best of these three. Two adjacent questions remain:

(a) Does the ME-specific BGI from `ME_BGI_V1.tif` give a different signal?
    BGI is a regression-residual-style index measuring how much actual growth
    exceeds or falls short of expected. As a plot-level dDBH.mult it could
    explain plot-level variance even if its mean equals 1.0. Task #183
    tracks the raster extraction and v27c follow-up.

(b) Does scaling ClimateSI per the v25 finding (CSI x 0.7) and then layering
    per-species calibration on top get below +9 percent BA bias? The
    decomposition in v18 + v25 + v26 + v27b suggests this should give
    approximately +9.6 pp - some calibration interaction. A v28 test
    combining the two levers would close that question.

## Status of the residual under 12.3.8

After this scan, the picture is:

  ClimateSI is the right choice among in-vdat metrics      ~0.3 pp loss if swapped
  CSI scaling (v25, x 0.6 to x 0.7)                        ~1.5 pp partial closure
  Mortality (v26 mort.mult)                                null lever
  Per-species calibration (v18, in production)             cancelled by MORTCAL
  Ingrowth (v24, in production)                            ~0.2 pp
  Productivity metric choice (v27b)                        no closure available

Remaining structural ~9 pp BA bias is in the intrinsic Kuehne et al. dDBH
coefficients or the Glover/Hool mortality functional form. Both are
paper-sized commitments and outside autopilot scope.

## Files

`cardinal_acadgy_siteindex_v27b.R` and `acdgy_siteindex_v27b_results.csv` in
`acd_divergence_2026-05-21/`. Cardinal run as SLURM job 10873744
(approximately 10 minutes on c0320).
