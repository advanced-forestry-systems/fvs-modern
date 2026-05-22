# Addressing the ingrowth residual: root cause localized

2026-05-22. The residual FIA basal area over-projection after MORTCAL (+8.6 to
+10.9%) is a QMD / small-tree effect (#127): the model QMD runs about 3.7% above
observed (5.15 vs 4.97), because observed stands gain small-tree ingrowth that
the model does not produce. This note pins down why the model recruits nothing.

## The recruitment path is inert across every documented lever

FIA validation (100 plots, 12.3.6), varying the ingrowth controls:

| config | BA bias | TPA | QMD (obs 4.97) |
|---|---|---|---|
| INGROWTH off (v18) | +10.9% | 1019.32 | 5.150 |
| INGROWTH on, CutPoint 0.95 | +10.9% | 1019.32 | 5.150 |
| INGROWTH on, CutPoint 0.50 | +10.9% | 1019.32 | 5.150 |
| INGROWTH on, CutPoint 0 | +10.9% | 1019.32 | 5.150 |

All four are byte-identical (TPA 1019.32299652593 in every case). Toggling
INGROWTH and lowering CutPoint from 0.95 to 0.50 to 0 changes nothing: the model
adds zero recruits.

## Why (code level)

- `ING.TreeList` (line 1565) recruits only on plots with `IPH > 0`.
- `Ingrowth.FUN` (Li et al. 2011 GNLS) computes the ingrowth probability
  `PI = logistic(a0 + a1*BA + a2*PHW + ... )`. For a typical Acadian FIA stand
  (BA ~30 m2/ha, PHW ~0.5, CSI ~14, MinDBH 3) `PI` is about 0.42.
- The gate (line 1503) is a HARD threshold:
  `IPH = if (CutPoint == 0) IPH*PI else ifelse(PI >= CutPoint, IPH, 0)`.
  With CutPoint 0.95 or 0.50, `PI` (0.42) is below both, so `IPH = 0` and no
  plot recruits. This explains why 0.95 and 0.50 are identical.
- `CutPoint == 0` switches to the expected value `IPH*PI` (> 0), which SHOULD
  recruit on every plot. It does not change the result, so the recruited trees
  are not being retained into the carried-forward tree state in the standalone
  AcadianGYOneStand path. CutPoint propagation is correct (line 1914 reads
  `ops$CutPoint`, line 2501 passes it to `Ingrowth.FUN`), so the zero is
  structural in the ingrowth-to-tree-list-to-next-cycle handoff, not a
  mis-set parameter.

## The fix (model-code work, #127)

Two steps, in order:

1. Make recruitment actually carry through. Trace why `ING.TreeList` output
   (with `CutPoint = 0`, expected-value mode) does not appear in the projected
   stand. Confirm the `ingrow` rows are bound (line 2586) AND survive the
   `rtnVars` subset and the next annual cycle. This is a standalone-path gap; in
   the customRun bridge ingrowth is added separately via `fvsAddTrees`, which is
   why it was never exercised here.
2. Once recruitment carries through, recalibrate the Li et al. probability /
   count for FIA-like conditions. A hard `PI >= CutPoint` gate with `PI ~ 0.42`
   is too restrictive; the expected-value formulation (CutPoint = 0) is the
   correct way to apply this hurdle model for prediction, but the coefficients
   should be checked against observed FIA recruitment so the produced ingrowth
   matches the observed small-tree gain (which would pull QMD from 5.15 toward
   4.97 and remove the residual BA over-projection).

This is the genuine stand-level fix. Diameter and mortality multipliers cannot
close the residual (shown earlier: diameter calibration makes BA worse, MORTCAL
plateaus at +8.6%). The lever is recruitment.

## Files

`cardinal_acadgy_cutpoint_v19.R` / `acadgy_cutpoint_v19_results.csv` (CutPoint
0.95 vs 0.50) and `cardinal_acadgy_cp0_v20.R` / `acadgy_cp0_v20_results.csv`
(CutPoint 0 expected-value vs hard threshold).
