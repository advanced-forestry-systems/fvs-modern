# FVS four-arm result, arm D complementarity test (all variants)

Date 2026-06-18. Branch conus-sf-integration-2026-05-21, holoros/fvs-modern. Source diagnostics_2026-06-16/fourarm_abcd_20260618.csv, SLURM jobs 11745408 (gen_baimult) and 11745409 (fourarm_abcd), both COMPLETED.

## Design

Four arms run inside the FVS engine, evaluated out of sample on spatial fold B with bootstrap confidence intervals. Arm A is the default engine. Arm B adds the brms maxSDI ceiling, density-dependent recruitment, and the fold-A BAIMULT. Arm C applies per-species BAIMULT against observed and default diameter growth, emulating the fvs-conus growth equations. Arm D combines C with the brms maxSDI and density recruitment layer. Metrics are percent bias in basal area (BA), trees per hectare (TPH), quadratic mean diameter (QMD), and merchantable volume (VOL). Twenty one variants entered the analysis. The em variant dropped out because its fold B held too few plots (A=8, B=17).

## Median |bias| across variants

| Metric | A default | B density | C BAIMULT | D combined | Winner |
|---|---|---|---|---|---|
| Basal area | 6.1 | 4.9 | 6.7 | 5.6 | B |
| Trees per hectare | 16.3 | 9.3 | 17.2 | 9.3 | B |
| QMD | 12.5 | 5.9 | 15.1 | 5.7 | D |
| Merch volume | 12.0 | 10.7 | 13.0 | 11.9 | B |

## Complementarity verdict

The combined arm D does not stack the two layers. D beats both B and C only on QMD, and there by a margin (5.7 versus B 5.9) that sits inside the noise. On BA, TPH, and VOL the density layer alone (B) is the best or tied-best arm, and adding the C growth signal moves the median the wrong way. The per-species BAIMULT signal (C) does not help at the all-variant median: C runs worse than the default A on BA, QMD, and VOL, so D mostly tracks B and inherits little benefit from C. The hypothesis that the fvs-conus growth signal complements the density layer is not supported by these medians. The workhorse is the density layer; the growth-equation emulation adds cost without a median gain.

Arm by arm, B wins BA (4.9), B and D tie on TPH (9.3), D edges B on QMD (5.7 versus 5.9), and B wins VOL (10.7).

## Basal area (BA) by variant, with bootstrap CIs for B and D

| Variant | A | B | C | D | B 95% CI | D 95% CI |
|---|---|---|---|---|---|---|
| ne | +9.2 | +8.1 | +8.9 | +8.0 | [+6.2,+10.1] | [+6.1,+9.9] |
| acd | +7.5 | +6.8 | +6.7 | +7.6 | [+5.3,+8.3] | [+6.1,+9.1] |
| sn | +12.2 | +9.4 | +19.7 | +20.4 | [+4.6,+14.1] | [+15.4,+25.5] |
| ls | +12.6 | +7.2 | +12.7 | +9.7 | [+3.2,+11.5] | [+5.6,+14.2] |
| cs | +19.0 | +15.1 | +20.3 | +19.2 | [+9.0,+23.0] | [+12.9,+27.3] |
| ie | -5.4 | -4.9 | -2.7 | -2.8 | [-8.1,-1.1] | [-5.9,+0.7] |
| kt | +14.0 | +7.6 | +6.8 | +2.6 | [+3.3,+12.1] | [-1.2,+6.5] |
| ci | -0.5 | -2.1 | -1.4 | -2.4 | [-6.0,+2.3] | [-6.1,+1.8] |
| bm | +1.5 | +3.9 | -1.2 | +1.9 | [-0.5,+8.8] | [-1.9,+6.2] |
| cr | +23.2 | +20.3 | +18.9 | +16.9 | [+11.8,+32.7] | [+8.7,+28.9] |
| tt | -3.0 | -0.5 | -2.5 | +0.0 | [-4.6,+4.0] | [-4.1,+4.6] |
| ut | +39.9 | +35.0 | +40.2 | +35.2 | [+13.8,+70.8] | [+14.1,+71.1] |
| ca | +5.4 | +3.3 | +6.6 | +4.2 | [-0.8,+7.9] | [+0.1,+9.0] |
| ws | +3.2 | +4.0 | +4.2 | +5.0 | [-0.2,+9.0] | [+0.6,+10.0] |
| nc | +6.1 | +7.9 | +5.4 | +7.5 | [+2.8,+14.4] | [+2.5,+14.0] |
| so | -4.6 | -1.4 | +0.0 | +3.0 | [-6.1,+4.7] | [-2.0,+9.7] |
| ec | -12.4 | -0.7 | -11.2 | +1.6 | [-4.5,+3.4] | [-2.1,+5.6] |
| wc | +3.8 | +3.1 | +6.5 | +5.6 | [-0.1,+6.7] | [+2.5,+9.1] |
| oc | -2.0 | -2.2 | +4.1 | +3.6 | [-5.9,+1.8] | [+0.2,+7.3] |
| op | +5.9 | +4.7 | +13.1 | +12.3 | [-2.0,+11.6] | [+5.8,+19.1] |
| pn | +10.4 | +9.0 | +7.1 | +6.1 | [+5.8,+12.8] | [+2.9,+9.7] |
| **median |bias|** | 6.1 | 4.9 | 6.7 | 5.6 | | |

## Trees per hectare (TPH) by variant, with bootstrap CIs for B and D

| Variant | A | B | C | D | B 95% CI | D 95% CI |
|---|---|---|---|---|---|---|
| ne | -13.0 | -11.2 | -12.7 | -10.7 | [-15.7,-6.4] | [-15.3,-6.0] |
| acd | -19.1 | -8.8 | -18.2 | -9.2 | [-13.0,-4.1] | [-13.4,-4.6] |
| sn | -13.6 | -3.7 | -16.0 | -8.1 | [-10.7,+3.8] | [-14.8,-0.9] |
| ls | -0.5 | -6.7 | -0.5 | -8.0 | [-15.5,+2.2] | [-16.5,+0.7] |
| cs | +2.5 | +1.3 | +2.1 | -0.0 | [-8.5,+12.8] | [-9.7,+11.3] |
| ie | -28.4 | -15.4 | -28.2 | -15.5 | [-24.5,-3.9] | [-24.6,-4.1] |
| kt | -28.7 | -15.5 | -29.3 | -15.1 | [-26.8,-3.1] | [-26.4,-2.6] |
| ci | -29.7 | -21.0 | -29.8 | -21.1 | [-29.8,-11.7] | [-30.0,-11.8] |
| bm | -25.9 | -16.5 | -24.3 | -14.3 | [-25.4,-6.6] | [-22.9,-4.4] |
| cr | +13.0 | +17.1 | +14.3 | +20.1 | [-1.4,+49.3] | [+2.0,+51.5] |
| tt | -17.1 | +4.9 | -17.2 | +4.8 | [-5.2,+16.5] | [-5.3,+16.4] |
| ut | -8.2 | -11.8 | -8.2 | -11.9 | [-27.9,+21.5] | [-27.9,+21.3] |
| ca | -2.2 | -2.0 | -2.2 | -2.1 | [-14.4,+14.3] | [-14.7,+13.9] |
| ws | -2.2 | +3.9 | -2.5 | +3.5 | [-9.8,+19.8] | [-10.0,+19.3] |
| nc | -16.3 | -3.8 | -16.3 | -3.6 | [-19.0,+13.8] | [-18.7,+14.0] |
| so | -15.9 | +1.9 | -19.6 | -2.4 | [-14.1,+20.5] | [-17.9,+15.1] |
| ec | -32.1 | -11.7 | -36.1 | -14.7 | [-20.1,-1.9] | [-23.6,-4.6] |
| wc | -22.5 | -20.6 | -22.5 | -21.4 | [-28.8,-10.9] | [-29.7,-11.8] |
| oc | -18.3 | -9.3 | -18.4 | -9.3 | [-17.9,+0.3] | [-17.9,+0.2] |
| op | -13.1 | -6.4 | -13.1 | -4.4 | [-22.4,+10.9] | [-20.3,+12.0] |
| pn | -22.6 | -22.7 | -22.6 | -21.7 | [-30.9,-13.1] | [-30.0,-12.2] |
| **median |bias|** | 16.3 | 9.3 | 17.2 | 9.3 | | |

## QMD (QMD) by variant, with bootstrap CIs for B and D

| Variant | A | B | C | D | B 95% CI | D 95% CI |
|---|---|---|---|---|---|---|
| ne | +11.0 | +1.1 | +11.0 | +1.0 | [-2.4,+4.8] | [-2.5,+4.7] |
| acd | +10.1 | -4.7 | +9.4 | -3.9 | [-7.9,-1.7] | [-7.1,-0.8] |
| sn | +11.2 | -1.8 | +16.6 | +5.7 | [-7.1,+3.6] | [-0.2,+11.9] |
| ls | +6.9 | +5.9 | +6.9 | +7.9 | [+0.9,+12.4] | [+2.8,+14.4] |
| cs | +11.0 | +10.1 | +11.9 | +12.8 | [+0.1,+21.4] | [+2.4,+24.2] |
| ie | +15.8 | -5.3 | +17.1 | -4.3 | [-11.1,+0.6] | [-10.1,+1.4] |
| kt | +21.5 | -6.1 | +17.8 | -8.8 | [-12.0,+0.7] | [-14.5,-2.1] |
| ci | +13.1 | -5.3 | +12.7 | -5.7 | [-10.4,+0.2] | [-10.6,-0.2] |
| bm | +23.7 | +4.0 | +21.0 | +1.8 | [-4.4,+14.6] | [-6.2,+12.2] |
| cr | +2.8 | -13.0 | +0.6 | -15.2 | [-20.4,-4.3] | [-22.4,-6.8] |
| tt | +12.5 | -11.1 | +12.9 | -10.8 | [-15.0,-6.9] | [-14.7,-6.6] |
| ut | +18.7 | -0.8 | +18.9 | -0.7 | [-14.5,+20.9] | [-14.3,+21.3] |
| ca | +4.8 | -10.9 | +5.4 | -10.2 | [-17.1,-3.6] | [-16.6,-2.8] |
| ws | +6.3 | -11.1 | +6.5 | -10.8 | [-17.2,-3.8] | [-17.1,-3.6] |
| nc | +20.3 | -5.9 | +19.8 | -6.4 | [-16.7,+4.5] | [-17.0,+4.0] |
| so | +13.3 | -12.9 | +17.6 | -9.2 | [-22.7,-3.3] | [-19.3,+1.1] |
| ec | +23.0 | -4.7 | +26.3 | -2.4 | [-9.5,+0.0] | [-7.3,+2.4] |
| wc | +15.5 | -1.0 | +16.6 | +0.3 | [-6.2,+4.2] | [-5.0,+5.5] |
| oc | +10.7 | -7.7 | +13.9 | -4.9 | [-12.5,-2.9] | [-9.8,-0.1] |
| op | +11.4 | -10.1 | +15.1 | -7.1 | [-17.8,-1.2] | [-15.2,+2.2] |
| pn | +19.0 | +2.7 | +17.2 | +0.8 | [-2.8,+8.2] | [-4.5,+6.1] |
| **median |bias|** | 12.5 | 5.9 | 15.1 | 5.7 | | |

## Merch volume (VOL) by variant, with bootstrap CIs for B and D

| Variant | A | B | C | D | B 95% CI | D 95% CI |
|---|---|---|---|---|---|---|
| ne | +19.9 | +19.4 | +19.5 | +19.1 | [+16.5,+22.2] | [+16.2,+21.9] |
| acd | +14.5 | +10.9 | +13.0 | +12.6 | [+8.7,+13.1] | [+10.4,+14.9] |
| sn | +7.9 | +2.3 | +17.8 | +17.3 | [-4.1,+9.4] | [+10.6,+25.2] |
| ls | +17.8 | +12.3 | +18.0 | +17.1 | [+6.4,+19.5] | [+11.1,+24.8] |
| cs | +24.4 | +19.9 | +25.9 | +25.3 | [+9.4,+32.1] | [+14.5,+38.0] |
| ie | -7.2 | -7.9 | -4.3 | -5.7 | [-14.8,-0.8] | [-12.0,+1.1] |
| kt | +8.5 | +0.8 | +1.7 | -4.2 | [-4.9,+6.7] | [-9.5,+1.1] |
| ci | -12.0 | -13.6 | -11.9 | -12.8 | [-19.4,-7.4] | [-18.4,-7.0] |
| bm | -10.7 | -6.8 | -11.6 | -7.4 | [-11.6,-1.5] | [-12.0,-2.6] |
| cr | +29.5 | +27.0 | +25.0 | +23.2 | [+15.2,+43.6] | [+11.8,+39.0] |
| tt | -13.8 | -11.4 | -13.4 | -11.0 | [-15.9,-7.2] | [-15.5,-6.7] |
| ut | +51.3 | +50.3 | +52.7 | +51.4 | [+9.0,+147.2] | [+9.2,+150.0] |
| ca | +11.6 | +10.7 | +12.9 | +11.9 | [+1.2,+20.8] | [+2.0,+22.2] |
| ws | -3.8 | -2.6 | -2.9 | -1.5 | [-8.2,+4.2] | [-7.1,+5.2] |
| nc | -7.8 | -4.4 | -8.2 | -4.7 | [-9.1,+2.1] | [-9.4,+1.8] |
| so | -11.2 | -8.7 | -7.3 | -4.6 | [-14.0,-2.1] | [-9.1,+1.6] |
| ec | -19.6 | -9.9 | -18.8 | -8.3 | [-15.3,-3.9] | [-13.7,-2.3] |
| wc | -2.2 | -2.3 | -0.8 | -0.9 | [-7.8,+3.9] | [-6.2,+5.0] |
| oc | +20.5 | +19.8 | +27.9 | +26.4 | [+13.1,+25.8] | [+20.1,+32.4] |
| op | +17.9 | +17.5 | +25.2 | +24.9 | [+4.3,+29.3] | [+12.3,+36.5] |
| pn | +11.1 | +10.7 | +8.4 | +8.2 | [+3.6,+18.9] | [+1.4,+16.0] |
| **median |bias|** | 12.0 | 10.7 | 13.0 | 11.9 | | |

## Variants where arm D regresses

Per-variant cases where D carries a larger absolute bias than the default A on a metric (threshold 0.5 points):

- sn: BA 12 to 20; VOL 8 to 17
- ls: TPH 0 to 8; QMD 7 to 8
- cs: QMD 11 to 13; VOL 24 to 25
- ci: BA 0 to 2; VOL 12 to 13
- cr: TPH 13 to 20; QMD 3 to 15
- ut: TPH 8 to 12
- ca: QMD 5 to 10
- ws: BA 3 to 5; TPH 2 to 4; QMD 6 to 11
- nc: BA 6 to 8
- wc: BA 4 to 6
- oc: BA 2 to 4; VOL 20 to 26
- op: BA 6 to 12; VOL 18 to 25

The recurring pattern is volume and basal area drift in southern and western variants (sn, cs, oc, op) once the density ceiling pulls stems out, and QMD overshoot in a few small-fold variants (cr, ws, ca) where the recruitment layer is poorly constrained.

## Figure

![Four-arm out-of-sample bias by metric](fourarm_abcd_20260618.png)
