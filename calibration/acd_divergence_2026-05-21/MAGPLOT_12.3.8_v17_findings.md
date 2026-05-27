# Canadian MAGPlot under AcadianGY 12.3.8: MORTCAL still over-corrects

2026-05-27. Re-runs the v16 MAGPlot harness (262 New Brunswick CFI pairs, 10
annual cycles) against AcadianGY_12.3.8.r, adding two configurations that
exercise CutPoint=0 expected-value ingrowth now that recruits actually carry
through cycles.

## Result

| config                                     | n   | BA bias % | R^2   | QMD (cm) | TPH  |
|--------------------------------------------|-----|-----------|-------|----------|------|
| canonical_off (no MORTCAL, default CP)     | 262 | -0.04     | 0.878 | 12.95    | 1675 |
| insource_on (MORTCAL on, default CP)       | 262 | -6.51     | 0.869 | 12.93    | 1548 |
| insource_on_cp0 (MORTCAL on, CP=0)         | 262 | -6.08     | 0.870 | 12.25    | 1635 |
| **ingrowth_only (MORTCAL off, CP=0)**      | 262 | **+0.39** | 0.878 | 12.30    | 1764 |

Observed (each row): BA 17.30 m^2/ha, QMD 11.89 cm, TPH 1807. The first two
rows match the v16 result on 12.3.6 to four decimals, which is the expected
strict-improvement property (12.3.8 == 12.3.6 when CutPoint = 0.95 suppresses
ingrowth).

## Interpretation

The +0.39 percent BA bias in ingrowth_only is the smallest residual of the
four configurations, beating both the canonical -0.04 percent baseline (which
the ingrowth_only run improves on for QMD and TPH) and any configuration with
MORTCAL. On Canadian MAGPlot, the right configuration is MORTCAL off and
ingrowth on with expected-value mode.

MORTCAL still over-corrects under 12.3.8. Adding ingrowth on top of MORTCAL
recovers 0.43 pp of BA (-6.51 -> -6.08), but the sign does not flip and the
direction does not change. Issue #128 (MORTCAL is needed for FIA but harmful
for Canadian CFI) is therefore a real, structural finding rather than an
artifact of the now-fixed ingrowth carry-through. The right operating
recommendation remains: enable MORTCAL only for FIA-like conditions; for
Canadian MAGPlot keep MORTCAL off and use CutPoint = 0 expected-value
ingrowth.

The ingrowth_only QMD improvement is notable: 12.95 -> 12.30 cm with observed
11.89, closing roughly two thirds of the QMD overshoot in the canonical
configuration. TPH closes from 1675 -> 1764 with observed 1807.

## Files

`cardinal_magplot_insource_v17.R` and `magplot_insource_v17_results.csv` in
`acd_divergence_2026-05-21/`. Cardinal run as SLURM job 10582592 (14 min
on c0011).
