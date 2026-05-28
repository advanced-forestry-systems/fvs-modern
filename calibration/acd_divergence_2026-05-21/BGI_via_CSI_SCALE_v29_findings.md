# BGI routed through ops$CSI_SCALE per stand (v29): null lever, definitive

2026-05-28. v27c found BGI as a per-tree dDBH.mult was null. 12.3.9 exposes
`ops$CSI_SCALE`, which routes climate weighting through the proper channel
that feeds dDBH, dHt, and ingrowth equations. v29 tests whether BGI's
plot-level signal, routed through CSI_SCALE instead of dDBH.mult, captures
variance the per-tree approach missed.

## Result

Same 100 ME FIA plots, 10 yr, 12.3.9 production posture (MORTCAL on,
CutPoint = 0 expected-value ingrowth):

| config              | CSI_SCALE used per stand                | BA bias % | R^2    | TPA  | QMD cm |
|---------------------|-----------------------------------------|-----------|--------|------|--------|
| csi_scale_1.0       | unset (matches v24)                     | +11.05    | 0.4232 | 1043 | 4.923  |
| **csi_scale_0.7**   | **0.7 uniform**                         | **+9.96** | **0.4325** | **1035** | **4.940** |
| csi_bgi_recenter    | BGI / mean(BGI)                         | +11.04    | 0.4230 | 1043 | 4.915  |
| csi_bgi_0.7x        | 0.7 * BGI / mean(BGI)                   | +9.99     | 0.4308 | 1035 | 4.933  |

BGI factor range across the 100 plots: 0.319 to 1.285. So some plots had
climate weighted as low as 32 percent, others as high as 129 percent. Plenty
of variance to exploit, if BGI were carrying useful signal.

## Three findings

(a) **csi_bgi_recenter is null.** BGI variance routed through the climate
    channel produces no improvement (+11.04 vs +11.05; R^2 0.4230 vs 0.4232).
    Identical conclusion to v27c per-tree dDBH.mult, now via the cleaner CSI
    path. The signal is genuinely not there.

(b) **csi_bgi_0.7x slightly degrades the uniform 0.7 result.** Stacking BGI
    variance on top of the global 0.7 shrink gives R^2 0.4308 vs uniform
    0.4325 (-0.0017). The variance term is not just neutral, it adds a small
    amount of noise to the fit.

(c) **The lever is in the mean, not the variance.** Uniform `CSI_SCALE = 0.7`
    beats both BGI variants. The +9 percent residual at CSI_SCALE = 0.7 is
    structural in the model's diameter and mortality equations, not a missing
    productivity signal across plots.

## What this rules out for closing the residual

After this scan, all candidates for "BGI carries useful information for the
+11 percent FIA BA residual" are dead:

  - Per-tree dDBH.mult (v27c bgi_recenter): null
  - Per-stand CSI replacement (v27c bgi_as_dmult): pathological (R^2 -10.44)
  - Per-stand ops$CSI_SCALE (v29 csi_bgi_recenter): null
  - Stacked with uniform 0.7 (v29 csi_bgi_0.7x): slightly worse

The remaining productivity question worth asking on this data is whether BGI
predicts the residual when applied as a post-hoc correction (as it does in
the fia_cem_projections pipeline). That would be the multiplicative
correction `corrected_BA = predicted_BA * (BGI / ref)^str` applied AFTER
projection. v30 could test that, but the per-plot signal demonstrated by
v27c/v29 suggests it would also be marginal.

The cleaner takeaway: **stop layering BGI into AcadianGY**. The lever is the
intrinsic diameter coefficients, the mortality functional form, or a
variant-level ceiling, none of which are external productivity inputs.

## Files

`cardinal_acadgy_bgicsi_v29.R`, `acdgy_bgicsi_v29_results.csv` and
`patch_v29.py` in `acd_divergence_2026-05-21/`. Cardinal run as SLURM job
10886783 (about 18 minutes on c0319).
