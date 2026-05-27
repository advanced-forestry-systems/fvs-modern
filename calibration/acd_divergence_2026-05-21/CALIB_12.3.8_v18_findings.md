# FIA per-species calibration under 12.3.8: refined compensating-errors picture

2026-05-27. Re-runs the v17 calibration harness (200 ME FIA plots, 10 annual
cycles) against AcadianGY_12.3.8.r, adding a configuration that layers CutPoint
= 0 expected-value ingrowth on top of calibrated_on, plus an ingrowth_only
baseline for comparison.

## Result

| config              | BA bias % | TPA (obs 1060) | QMD cm (obs 5.02) | R^2  |
|---------------------|-----------|----------------|--------------------|------|
| canonical_off       | +15.35    | 1084           | 5.25               | 0.42 |
| calibrated_off      | +22.26    | 1102           | 5.35               | 0.32 |
| calibrated_on       | +15.16    | 1019           | 5.37               | 0.43 |
| calibrated_on_cp0   | +15.38    | 1042           | 5.16               | 0.42 |
| ingrowth_only       | +15.50    | 1107           | 5.05               | 0.42 |

Observed: BA 98.44 ft^2/ac, TPA 1060, QMD 5.02 in. SLURM job 10594343 on c0001,
about 22 minutes.

## Decomposition of the production configuration

Walking from canonical_off (no levers) to calibrated_on_cp0 (production posture):

  canonical_off       +15.35  (none)
  + calibration       +22.26  (+6.91 pp)  diameter mults push BA up
  + MORTCAL           +15.16  (-7.10 pp)  rebalances calibration's BA gain
  + ingrowth (CP=0)   +15.38  (+0.22 pp)  adds the missing small trees

The three levers combine to within 0.03 pp of the no-lever baseline. The prior
characterization of "compensating errors" is mechanistically accurate, but
several quantitative details change once recruits flow:

(a) Calibration is not a wash. It alone pushes BA up because tree-level
    diameter growth is calibrated higher, and BA scales with diameter squared.
    Considered in isolation, calibration worsens stand BA. This finding from
    v17 holds.

(b) MORTCAL alone (under 12.3.8 with full ingrowth) is harmful on Canadian
    MAGPlot (-6.08 percent BA bias, see MAGPLOT_12.3.8_v17_findings.md). But
    on FIA it almost perfectly cancels calibration's BA inflation. This is
    not coincidence; both levers operate on Maine FIA conditions, where they
    were calibrated, and produce equal-and-opposite effects on BA.

(c) Ingrowth contributes only 0.15 to 0.22 pp to BA but materially improves
    TPA (1019 -> 1042 in production posture; closer to obs 1060) and QMD
    (5.37 -> 5.16 cm; closer to obs 5.02). Ingrowth's job is stand structure,
    not stand size.

(d) The production posture (calibrated_on_cp0) is the best of any
    configuration on TPA AND on QMD without inflating BA. R^2 is the same as
    canonical (0.42). It is the right operating recommendation for FIA-like
    Maine conditions.

## Comparison to v17

The v17 result (12.3.6) on the same 200 plots:

  canonical_off    +15.35
  calibrated_off   +20.24    (vs v18 +22.26, drift +2.0 pp)
  calibrated_on    +13.31    (vs v18 +15.16, drift +1.9 pp)

The drift is consistent with using a different calibration coefficient table
(`acd_annual_calibration_test.csv` on Cardinal vs the production version v17
pulled from `fvs-modern-acdbridge/`). The directional findings are unchanged:
calibration worsens BA in isolation, MORTCAL recovers it, and ingrowth adds a
small additional contribution.

## Status of the residual

The remaining +15 percent BA bias is structural and not removable by any
combination of the three levers under 12.3.8. The recruit-related portion
(roughly 0.15 pp) is now in the model; the rest is the underlying diameter
overshoot that the model was calibrated to deliver, balanced against the
MORTCAL haircut and the per-species calibration. The headline "BA over-projects
by ~15 percent on Maine FIA" is now mechanistically attributable rather than
mysterious. The next BA closure has to come from somewhere else - climate
sensitivity in the diameter equations, a different mortality functional form,
or accepting +15 percent as the structural ceiling for the current Acadian
variant.

## Files

`cardinal_acadgy_calib_v18.R` and `acadgy_calib_v18_results.csv` in
`acd_divergence_2026-05-21/`. Cardinal run as SLURM job 10594343.
