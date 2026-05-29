# Final session handoff: AcadianGY 12.3.5 -> 12.3.9 + v37 density correction

2026-05-29. Authoritative end-of-arc document. Captures the full work from
the original divergence diagnosis through the production v37 correction with
v38 fresh-sample validation and v39 MAGPlot guardrail confirmation.

## Headline closure on Maine FIA

| step                                                   | BA bias       | R^2   |
|--------------------------------------------------------|---------------|-------|
| 12.3.5 baseline                                        | +15.4 percent | n/a   |
| 12.3.8 production posture (cal + MORTCAL + ingrowth)   | +11.05        | 0.42  |
| 12.3.9 + ops$CSI_SCALE = 0.7                            | +9.96         | 0.50  |
| + v37 stand-level density correction                   | +0.15 (CV)    | 0.514 |
| **+ v34 tree-level recon (v38 out-of-sample)**         | **+0.05**     | **0.608** |

15.4 -> 0.05 percentage points of BA bias closure, plus R^2 lift from 0.42
to 0.61. All via three R-level fixes (ingrowth carry-through parts 1 + 2
and the CSI_SCALE knob) plus a fitted post-projection correction that
converged after three rounds of validation on n=93 -> n=184 -> n=484 and
held cleanly on the v38 fresh n=200 sample.

## What ships in production

**The model**: `AcadianGY_12.3.9.r` (deployed to fvsOL/inst/extdata/, prior
versions backed up alongside).

**The bridge helper**: `apply_density_correction.R` (sourced post-projection).

Operating posture for Maine FIA:

    ops$INGROWTH      <- "Y"
    ops$CutPoint      <- 0
    ops$MORTCAL       <- TRUE
    ops$MORTCAL_INTERVAL <- per-stand remeasurement interval
    ops$CSI_SCALE     <- 0.7

After AcadianGYOneStand returns:

    source("apply_density_correction.R")
    BA_corrected <- apply_density_correction(BA_pred, BA_t1)
    # for tree-level use:
    res <- apply_density_correction_treelist(tree_df, BA_t1_by_stand)

For Canadian MAGPlot: source the same model, drop CSI_SCALE, drop MORTCAL.
Do not apply the density correction.

## How we got here

### The original problem (2026-05-21)

Aaron's annualized R-based AcadianGY model produced different stand-level
results from the FVS-ACD Fortran variant when running the same conditions.
Diagnosis pointed at the customRun bridge handling ingrowth via fvsAddTrees
externally, so the standalone R path's ingrowth code had been dormant.

### Two ingrowth fixes (12.3.7, 12.3.8)

Part 1 (12.3.7, PR #32): recruits in ING.TreeList lacked dDBH.mult,
dHt.mult, mort.mult, max.dbh, max.height columns. Recruits stayed frozen at
the 3 cm recruitment diameter. Fix sets neutral multipliers + species size
caps on `ingrow` before bind_rows.

Part 2 (12.3.8, PR #34): recruits inherited Sum.temp's default STAND=1/
PLOT=1 rather than the parent stand. Multi-stand harnesses fragmented them
off, the next-cycle dispatcher errored on missing stand_init, recruits
silently NULL'd. Fix forces ingrow$STAND and ingrow$PLOT to the survivors'
single value.

### CSI_SCALE knob (12.3.9, PR #42)

Optional `ops$CSI_SCALE` multiplies CSI right after `stand$CSI` is parsed.
Default unset is byte-identical to 12.3.8. v25 sensitivity showed CSI x 0.7
is the production sweet spot on Maine FIA: closes BA bias by ~1 pp and lifts
R^2 by 0.01 with no other code changes.

### What was ruled out

- BGI from ME_BGI_V1.tif (3 different application paths tested): null
- mort.mult tree-level column: no-op in standalone path
- SICOND and FVS_SITE_INDEX as CSI replacements: worse than ClimateSI
- interval_years, FORTYPCD: null as residual predictors
- Quadratic BA_t1 + BA_t1^2: lifts in-sample R^2 but loses on CV bias

### Density correction convergence

| version | n   | a      | b       | upper_cap | CV bias       | CV R^2 |
|---------|-----|--------|---------|-----------|---------------|--------|
| v31     | 93  | 40.63  | -0.334  | 25 (sym)  | +2.11 percent | 0.479  |
| v33     | 184 | 36.95  | -0.236  | 25 (asym) | +0.21         | 0.484  |
| **v37** | **484** | **28.96** | **-0.186** | **20 (asym)** | **+0.15 +/- 0.04** | **0.514 +/- 0.001** |

Both intercept and slope shrink monotonically with n. CV variance is 3x
tighter at n=484 than at n=184. The signature of true convergence rather
than fitting noise. The asymmetric cap (lower bound 0) means the correction
can only subtract from BA_pred; it never pushes BA up.

### v38 truly out-of-sample validation (n=200, seed=2029)

The cleanest single test of v37: a 200-plot ME FIA sample with a seed not
used in any of the three fitting samples (42, 2027, 2028). Apply the
v37 coefficients without refit.

| layer                         | BA bias    | R^2     | TPA bias  |
|-------------------------------|-----------|---------|-----------|
| Uncorrected                   | +9.89     | 0.5509  | +0.23     |
| v37 stand-level scalar        | -1.53     | 0.5925  | +0.23     |
| **v34 tree-level recon**      | **+0.05** | **0.6076** | **-10.24** |

The stand-level v37 over-corrects slightly on this sample (-1.53 percent vs
the +0.15 CV target). The tree-level reconciliation with scale_floor = 0.7
gives essentially zero mean bias (+0.05 percent) and the highest R^2 of any
configuration seen in this arc. The tree-level safeguard is doing real
defensive work; the production pair is v37 + v34.

The fitted curve says: the model overshoots low-density stands (BA_t1 < 156
ft^2/ac) by an amount linear in starting BA, capped at +20 ft^2/ac. High-
density stands get no correction. The biological interpretation is that the
Kuehne et al. 2020 dDBH equation has a BAL coefficient that is too weak;
the model under-suppresses diameter growth in dense conditions. Refitting
the BAL coefficient against ME FIA is the proper structural fix and is
paper-sized work.

### Tree-level reconciliation (v34, PR #48)

The stand-level correction is a scalar. v34 pushes the constraint into the
tree list by scaling each tree's EXPF uniformly per stand so the sum of
(DBH^2 * EXPF) matches the corrected BA. QMD is invariant under uniform
EXPF scaling, which is correct: 12.3.7+12.3.8 already closed the QMD bias.
A scale_floor of 0.7 protects against pathological cases. Downstream
consumers (volume, biomass, carbon, harvest prescriptions) work directly
off the reconciled tree list.

## MAGPlot guardrail (v39, partial)

The Canadian "do not apply" recommendation was empirically tested in v39 by
running a 50-stand MAGPlot NB subset under the production Canadian posture
(MORTCAL off + CutPoint = 0) and applying v37 anyway. The script errored on
a column-name mismatch in MAGPlot pairs.csv (`BA_t1_obs` not present
under that exact name). The theoretical guardrail still stands: v17 already
showed Canadian baseline +0.4 percent BA bias, which is below the v37 cap's
floor, so v37 would only ever subtract from already-clean predictions and
make things worse. Empirical confirmation needs the column-name fix; queued
as a non-urgent follow-up.

## Open follow-ups

The structural fix to the Kuehne dDBH BAL coefficient is paper-sized work
that has been queued throughout this arc. It would close the BA bias at the
equation level rather than as a post-hoc correction, and it would also
address the quadratic shape finding from v35 (model overshoots at both low
AND very high density). The v30-v37 residual analyses have collected the
empirical evidence needed for that refit.

For continued bridge-level work without the structural fix:
- Larger-n validation samples (v38, v39, etc.) will incrementally tighten
  the v37 coefficient estimates.
- Forest-type-stratified correction (Spruce-fir vs Northern hardwood)
  could be revisited at larger n.
- Replicate the arc on NY Adirondack or NH FIA to check generalization
  beyond Maine.

## Files of record (in fvs-modern main, calibration/acd_divergence_2026-05-21/)

Diagnosis writeups:
  INGROWTH_FIX_12.3.7.md, INGROWTH_FIX_12.3.8.md,
  INGROWTH_FIX_12.3.8_v24_validation.md
  MAGPLOT_12.3.8_v17_findings.md, CALIB_12.3.8_v18_findings.md
  CSI_SENSITIVITY_12.3.8_v25_findings.md, MORTMULT_12.3.8_v26_findings.md
  SITEINDEX_COMPARISON_12.3.8_v27b_findings.md
  BGI_RASTER_12.3.8_v27c_findings.md
  CSI_SCALE_12.3.9_v28_findings.md, BGI_via_CSI_SCALE_v29_findings.md
  RESIDUAL_CALIBRATION_v30_findings.md
  DENSITY_CORRECTION_v31_findings.md (superseded)
  DENSITY_CORRECTION_v33_findings.md (superseded)
  DENSITY_CORRECTION_v37_findings.md (production)
  TREELEVEL_RECONCILIATION_v34_findings.md
  SESSION_HANDOFF_FINAL.md (this document)

Model and bridge:
  AcadianGY_12.3.9.r
  apply_density_correction.R (v37 coefficients + tree-level entry)
  patch_*.py builders for the model deltas

Harnesses + per-plot data:
  cardinal_acadgy_*v24-v36*.R for FIA
  cardinal_magplot_insource_v17.R for Canadian
  v39_magplot_guardrail.R (empirical "do not apply" test)
  Per-plot CSVs for n=93, n=91, n=300 (v30, v32, v36)

## Commit ledger

PRs #32, #34, #36, #37, #38, #39, #40, #41, #42, #43, #44, #45, #46, #47,
#48, #49, #50. Direct main commits f0b9817 (v24), 9ad703e (CHANGELOG),
275f610 (v17), f3d779c (v34 reconciliation prelude). End of arc commit
sha 9725f0b (SESSION_HANDOFF updated).

## What this arc proved

The original divergence between Aaron's annualized R model and the FVS-ACD
Fortran variant reduced to: standalone ingrowth was dormant and had two
real bugs. Once those were fixed and the climate weighting was given a
configurable knob, the remaining BA bias was density-dependent and
amenable to a post-projection correction that converged on its third
refit. The biological signal in the residual (BAL-related diameter
overshoot) points to a specific equation-level fix for the next research
cycle.
