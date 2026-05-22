# Wiring the annual calibration into the model: tree-level vs stand-level

2026-05-22. Experiment: feed the fitted per-species annual diameter multipliers
(`acd_annual_calibration.csv`, dDBH.mult / dHt.mult) into the live model and
re-validate against FIA, to see whether the calibration that improves per-tree
diameter growth also improves stand-level basal area. It does not. This is the
key finding and it redirects how the calibration should be used.

## FIA result (200 Maine plots, 12.3.6, MORTCAL handled in-source)

| config | BA bias | BA R2 | TPA (obs 1060) | QMD in (obs 5.02) |
|---|---|---|---|---|
| canonical (no calibration) | +15.4% | 0.417 | 1084 | 5.25 |
| calibrated diameter only | +20.2% | 0.356 | 1098 | 5.32 |
| calibrated + MORTCAL | +13.3% | 0.449 | 1015 | 5.34 |
| MORTCAL alone (v16) | +8.6% | 0.483 | 1002 | 5.26 |

## Interpretation

The fitted diameter multipliers are mostly above 1 because the raw Acadian
`dDBH_fun` under-predicts per-tree diameter growth (that fit is correct at the
tree level: it cut diameter-growth RMSE). But the FIA stand error is a basal area
OVER-projection driven by too little mortality, not by too little growth. So
applying the diameter multipliers boosts growth on top of an already too-high BA
and makes the stand bias worse (+15.4 to +20.2) while lowering R2.

In other words the model's stand BA is "right-ish" through compensating errors:
under-predicted diameter growth roughly offsets under-predicted mortality.
Calibrating the diameter component alone breaks that balance. Adding MORTCAL back
(calibrated + MORTCAL) pulls the bias to +13.3% with the best R2 of the three
(0.449), but it is still worse on mean bias than MORTCAL alone (+8.6%): the
diameter boost partially undoes the mortality correction.

## Implications for calibration strategy

1. The annualized diameter calibration serves a TREE-LEVEL objective (diameter
   growth accuracy). It should not be applied to improve STAND-LEVEL basal area;
   for stand BA it is counterproductive.
2. For stand BA on FIA-like Maine conditions, MORTCAL alone is best. On Canadian
   CFI (MAGPlot, already unbiased) the same diameter calibration would push BA
   positive (over-project), the mirror of the FIA result, so it is not a fix
   there either.
3. The genuine stand-level fix is the ingrowth / recruitment submodel (#127),
   not diameter or mortality multipliers. The residual after MORTCAL is a QMD /
   small-tree composition effect.
4. If both tree-level accuracy and stand-level BA are wanted simultaneously, the
   calibration must be done jointly at the stand level (constrained so diameter,
   mortality and ingrowth corrections are mutually consistent), not as an
   independent per-component multiplier.

## Files

`cardinal_acadgy_calib_v17.R` (the calibration-aware FIA harness) and
`acadgy_calib_v17_results.csv` (the table above).
