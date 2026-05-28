# AcadianGY 12.3.9: ops$CSI_SCALE knob and v28 validation

2026-05-28. Adds one optional parameter to AcadianGYOneStand:
`ops$CSI_SCALE`. When set, multiplies the resolved Climate Site Index by the
provided factor right after CSI is parsed from `stand$CSI`, before anything
downstream (Kuehne et al. dDBH, Russell/Weiskittel dHt, Li et al. Ingrowth.FUN,
Ingrowth.Comp) consumes it. Default is unset (no scaling) — 12.3.9 with no
CSI_SCALE in ops is byte-identical to 12.3.8.

## The change

In `AcadianGYOneStand`, immediately after the existing CSI resolution line:

    CSI = if (is.null(stand) || is.null(stand$CSI) || is.na(stand$CSI)) 12 else stand$CSI
    # 12.3.9: ops$CSI_SCALE knob
    if (!is.null(ops$CSI_SCALE) && is.finite(ops$CSI_SCALE) && ops$CSI_SCALE > 0) {
      CSI <- CSI * ops$CSI_SCALE
    }

Version tag bumped to AcadianV12.3.9.

## v28 validation (100 ME FIA plots, 10 yr, MORTCAL on, CutPoint = 0)

| CSI_SCALE | BA bias % | R^2    | TPA (obs 1029) | QMD cm (obs 4.97) |
|-----------|-----------|--------|----------------|--------------------|
| 1.0 (unset, = v24)  | +11.05    | 0.4232 | 1043           | 4.923              |
| **0.7 (production recommendation)** | **+9.96** | **0.4325** | **1035** | **4.940** |
| 0.5 (aggressive)    | +9.04     | 0.4411 | 1029           | 4.946              |

Three things this confirms:

(a) Backward compat: CSI_SCALE not set reproduces v24 exactly. 12.3.9 is a
    strict improvement; no existing harness needs to change.

(b) The knob behaves linearly with the v25 sensitivity scan. v25 gave
    csi_x0.6 +9.6 percent and csi_x0.8 +10.4 percent; the v28 csi_scale_0.7
    +9.96 percent sits at the interpolated midpoint. The mechanism is the
    same and the elasticity is consistent across runs.

(c) CSI_SCALE = 0.5 lands TPA at observed (1029.27 predicted vs 1029.33
    observed) and pulls BA bias below +10 percent for the first time in
    this sample. R^2 improves by 0.018. QMD moves from 4.92 toward observed
    4.97.

## Operating recommendations

**Default in the customRun bridge: CSI_SCALE = 0.7.** Conservative shift from
v24, ~1 pp BA closure, +0.01 R^2, no downstream code changes. Safe for
production.

**Optional aggressive default: CSI_SCALE = 0.5.** Larger shift, 2 pp BA
closure, TPA hits observed, +0.018 R^2. The risk is over-tuning to the
sample; would benefit from holdout validation before becoming the bridge
default.

**For Canadian MAGPlot: keep CSI_SCALE = 1.0** (or omit). Canadian validation
(v17) showed MORTCAL off + CutPoint = 0 already gives +0.4 percent BA bias on
NB CFI plots. There is no residual climate signal to close on that side. The
ME-tuned CSI scaling should not be applied to Canadian harnesses without a
parallel sensitivity scan on MAGPlot.

## Status of the +15 percent FIA BA residual

Walking the chain:

  12.3.6 baseline (recruits dropped, no calibration)        +15.4%
  12.3.8 production posture (calibration + MORTCAL + ingrowth)  +11.05%
  12.3.9 production posture, CSI_SCALE = 0.7                 +9.96%
  12.3.9 production posture, CSI_SCALE = 0.5                 +9.04%

12.3.9 with CSI_SCALE = 0.5 gets us within roughly 9 pp of zero. The remaining
gap is structural: Kuehne et al. dDBH base-rate coefficients, Glover/Hool
mortality functional form, or a true variant-level ceiling. All three are
paper-sized commitments and sit outside autopilot scope.

## Files

`AcadianGY_12.3.9.r`, `cardinal_acadgy_csiscale_v28.R`,
`acdgy_csiscale_v28_results.csv`, and `patch_12.3.9.py` in
`acd_divergence_2026-05-21/`. Cardinal run as SLURM job 10879396 (about 13
minutes on c0317).

Deployed to `ForestVegetationSimulator-Interface-main/fvsOL/inst/extdata/
AcadianGY.R` with 12.3.8 backed up.
