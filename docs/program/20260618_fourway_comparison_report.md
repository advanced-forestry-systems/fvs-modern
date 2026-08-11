# FVS-CONUS Four-Way Calibration Comparison

Default versus keyword-calibrated FVS across CONUS variants on disturbance-clean FIA remeasurement

Date: 2026-06-18
Job: SLURM 11759257 (calib_4way), OSC Cardinal, COMPLETED in 19 minutes, exit 0
Source script: diagnostics_2026-06-16/calib_4way.py
Aggregation: diagnostics_2026-06-16/calib_4way_aggregate.R
Repository: holoros/fvs-modern, branch conus-sf-integration-2026-05-21

## Summary

This benchmark projects every CONUS variant forward on COND-undisturbed FIA remeasurement pairs and compares two parameter sets against observed growth: the stock FVS defaults and the keyword-calibrated parameters. Bias is reported as percent signed bias, computed as 100 times the sum of predicted minus observed over the sum of observed, so a positive value means overprediction. Four response variables are evaluated: basal area (BA), trees per hectare (TPH), quadratic mean diameter (QMD), and volume (VOL).

Calibration reduces median absolute bias on every metric across the 951 evaluated conditions. The gain is largest for QMD and TPH and smallest for VOL.

| Metric | Median \|bias\| default | Median \|bias\| calibrated | Reduction |
|---|---|---|---|
| BA  | 11.4 % | 7.7 %  | 32 % |
| TPH | 16.7 % | 7.6 %  | 54 % |
| QMD | 9.8 %  | 4.1 %  | 58 % |
| VOL | 14.1 % | 10.9 % | 23 % |

The variant comparison rests on the full 951 conditions and is reliable. The species, ecoregion, and landowner breakdowns do not, because the per-condition CN key joined the projection output to the metadata keys at only 12.9 percent. Those strata are reported with that limitation stated plainly, and the root cause and fix are documented below.

## CN join shortfall (read before using the stratified results)

The aggregation joins the per-condition predictions in `calib_4way_20260618_percond.csv` to the four metadata keys (dominant species group, EPA Level II ecoregion, OWNGRPCD landowner) held in `conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds`. The join was attempted on the t1 plot key against several candidate columns:

| Candidate RDS column | Match rate |
|---|---|
| PLT_CN_cond1 | 12.9 % |
| CN_cond1 | 0 % |
| plot_key | 0 % |

Only PLT_CN_cond1 matched, and only for 123 of 951 conditions. The aggregation script proceeds when the best match clears 10 percent, so it ran rather than aborting, but the stratified output is built on that 123-condition subset. The species dimension produced no rows at all, because no species group reached the 20-condition minimum after the join. The ecoregion and landowner dimensions each retained only two strata.

Root cause: the per-condition CN emitted by `calib_4way.py` does not align with the CN keys stored in the remeasurement RDS for most conditions. The fix is to capture the species, ecoregion, and landowner keys directly inside `calib_4way.py` at projection time, carrying them through to the per-condition CSV, rather than recovering them by a post hoc CN join in the aggregation step. Until that change is made, treat everything below the variant section as indicative only.

## Bias by variant (full 951-condition basis)

Signed percent bias, default value then calibrated value, on the disturbance-clean basis. Variants are ordered by evaluated condition count.

| Variant | n | inj | BA def | BA cal | TPH def | TPH cal | QMD def | QMD cal | VOL def | VOL cal |
|---|---|---|---|---|---|---|---|---|---|---|
| acd | 145 | 1 | 6.3 | 5.0 | -18.0 | -12.3 | 9.9 | 5.0 | 13.1 | 10.9 |
| ne  | 112 | 1 | 9.9 | 9.7 | -11.4 | -4.4 | 8.7 | 2.9 | 24.0 | 22.9 |
| sn  | 97  | 1 | 15.2 | 14.0 | -13.8 | -0.5 | 14.9 | 5.7 | 14.1 | 9.6 |
| ec  | 77  | 1 | -13.2 | -1.5 | -30.8 | -6.1 | 14.7 | -2.2 | -16.9 | -6.4 |
| wc  | 77  | 1 | 4.5 | 2.3 | -16.7 | -14.6 | 6.4 | 0.7 | -2.2 | -3.6 |
| pn  | 77  | 1 | 11.4 | 7.7 | -16.9 | -17.1 | 9.8 | 4.1 | 12.3 | 10.1 |
| ls  | 64  | 0 | 16.8 | 13.2 | -4.1 | -9.9 | 7.2 | 7.7 | 21.1 | 18.4 |
| nc  | 54  | 1 | 6.0 | 8.5 | -10.7 | 4.4 | 18.4 | 9.4 | -8.0 | -5.2 |
| ie  | 47  | 1 | -0.9 | -2.9 | -24.0 | -13.0 | 10.3 | 2.1 | -11.6 | -14.6 |
| kt  | 45  | 1 | 12.4 | 5.3 | -40.2 | -16.5 | 26.3 | 3.9 | 18.6 | 7.1 |
| ci  | 44  | 1 | -0.3 | -4.1 | -22.0 | -8.0 | 8.7 | -1.8 | -11.4 | -15.5 |
| cs  | 40  | 0 | 20.4 | 18.3 | 6.3 | 4.0 | 3.8 | 3.6 | 20.5 | 18.5 |
| ca  | 40  | 0 | 7.6 | 4.6 | 3.9 | -1.4 | 4.8 | 4.9 | 8.4 | 7.0 |
| cr  | 19  | 0 | 25.9 | 20.3 | 7.6 | 1.1 | 4.4 | 4.2 | 33.6 | 28.8 |
| ut  | 13  | 1 | 13.1 | 11.4 | -20.6 | -7.6 | 13.8 | 4.6 | 54.5 | 49.9 |

Where calibration helps most. The largest corrections occur where the default carried the largest error. The kt variant is the clearest case: TPH moves from -40.2 to -16.5, QMD from +26.3 to +3.9, and VOL from +18.6 to +7.1. The ec variant improves on all four metrics, with BA moving from -13.2 to -1.5, TPH from -30.8 to -6.1, and QMD from +14.7 to -2.2. The sn variant nearly removes its TPH bias (-13.8 to -0.5) and more than halves QMD (+14.9 to +5.7). The ie and ci variants both cut QMD bias by roughly three quarters.

Where calibration does not help or overshoots. Calibration can degrade variants whose default bias was already near zero, because it pulls them past the target. The ci variant moves on BA from -0.3 to -4.1 and on VOL from -11.4 to -15.5; the ie variant moves on BA from -0.9 to -2.9 and on VOL from -11.6 to -14.6. The ls TPH bias worsens from -4.1 to -9.9, and the nc TPH bias flips sign from -10.7 to +4.4, an overcorrection rather than a reduction. VOL is the stubborn metric: ne, cs, ls, and cr all retain large positive volume bias after calibration (22.9, 18.5, 18.4, and 28.8 respectively), and ut remains at +49.9 on a sample of only 13.

Note on the inject flag. The per-variant output records an inject status (inj). Four variants ran with inj=0 (ls, cs, ca, cr) yet still show a default versus calibrated difference. That should be reconciled before these calibrated columns are quoted elsewhere, since inj=0 would normally imply no calibrated parameters were applied.

## Bias by EPA Level II ecoregion (123-condition subset, indicative only)

Two ecoregions cleared the 20-condition minimum.

| Ecoregion L2 | n | Metric | Default | Calibrated |
|---|---|---|---|---|
| 6.2 | 65 | BA  | 9.2 | 5.6 |
| 8.1 | 23 | BA  | 4.4 | 2.4 |
| 6.2 | 65 | TPH | -27.5 | -17.9 |
| 8.1 | 23 | TPH | -24.2 | -21.1 |
| 6.2 | 65 | QMD | 19.8 | 10.4 |
| 8.1 | 23 | QMD | 11.5 | 6.8 |
| 6.2 | 65 | VOL | 9.4 | 4.7 |
| 8.1 | 23 | VOL | 19.8 | 17.8 |

Ecoregion 6.2 (Western Cordillera) improves strongly on every metric, most of all QMD and VOL. Ecoregion 8.1 (Mixed Wood Plains) improves modestly, and its volume bias barely moves (19.8 to 17.8). All other EPA Level II regions fell below the 20-condition floor in the joined subset and cannot be assessed here.

## Bias by landowner (123-condition subset, indicative only)

OWNGRPCD groups National Forest (10), Other federal (20), State and local (30), and Private (40). Only National Forest and Private cleared the 20-condition minimum.

| Landowner | n | Metric | Default | Calibrated |
|---|---|---|---|---|
| Private | 55 | BA  | 8.2 | 6.5 |
| National Forest | 56 | BA  | 10.4 | 7.0 |
| Private | 55 | TPH | -17.1 | -10.0 |
| National Forest | 56 | TPH | -25.4 | -17.6 |
| Private | 55 | QMD | 12.1 | 5.0 |
| National Forest | 56 | QMD | 19.5 | 10.8 |
| Private | 55 | VOL | 11.7 | 9.0 |
| National Forest | 56 | VOL | 12.4 | 8.1 |

Calibration improves both classes on all four metrics. The gains are larger for National Forest, particularly on QMD (19.5 to 10.8) and TPH (-25.4 to -17.6), which also carried the larger default bias. Other federal and State and local land had too few joined conditions to report, so this comparison covers only the two best-sampled ownership classes.

## Most populated ecoregion by landowner cross

The cross is degenerate in this subset. Only two cells cleared the 20-condition minimum: National Forest within ecoregion 6.2 (calibrated QMD bias about +11) and Private within ecoregion 8.1 (calibrated QMD bias about +6). The two keys are therefore nearly collinear in the joined data, with National Forest concentrated in the western Cordillera and Private concentrated in the Mixed Wood Plains. Their separate contributions cannot be disentangled until the join coverage improves.

## Where calibration helps and where it does not

Calibration delivers its clearest benefit on QMD and TPH and on the variants that started with the largest default error (kt, ec, sn, ie, ci). It is least effective on volume, which retains substantial positive bias in several eastern and Rocky Mountain variants after calibration (ne, cs, ls, cr, ut). It can mildly worsen variants whose default was already accurate (ie and ci on BA and VOL, ls on TPH) and can overcorrect to a sign flip (nc on TPH). The stratified picture across ecoregion and landowner is consistent with the variant story, calibration helping most where bias was largest, but it rests on 123 of 951 conditions and should not be cited until the CN join is fixed by carrying the keys through `calib_4way.py` directly.

## Files

| File | Location (Cardinal) |
|---|---|
| Per-condition predictions | diagnostics_2026-06-16/calib_4way_20260618_percond.csv |
| Margin summary | diagnostics_2026-06-16/calib_4way_margins.csv |
| Landowner figure | diagnostics_2026-06-16/fig_4way_landowner.png |
| Ecoregion cross figure | diagnostics_2026-06-16/fig_4way_ecoregion.png |
| This report | docs/program/20260618_fourway_comparison_report.md |
