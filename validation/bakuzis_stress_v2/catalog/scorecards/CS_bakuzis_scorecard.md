# CS: Bakuzis scorecard (catalog v3, October 29, 2026)

Engine arms: 233 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 41.8 and p99.9 64.0 m2 ha-1, SDI p99.9 1130 trees ha-1 from 302 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **RECOMMENDED**; recommended arm default_native; admissible: default_native; default_gompit; calibrated_native; calibrated_gompit.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 233 / **4** | 233 / **1** | 233 / **46** | 233 / **19** | 1500 / **9** | 1500 / **17** | 1500 / **10** | 1500 / **17** |
| T1 monotone HT_dom | 233 / **29** | 233 / **53** | 233 / **28** | 233 / **86** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 233 / **0** | 233 / **1** | 233 / **5** | 233 / **5** | 964 / **110** | 829 / **650** | 1000 / **89** | 877 / **563** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 201 / **5** | 177 / **15** | 70 / **1** | 73 / **4** | 1435 / **1238** | 1348 / **620** | 1380 / **1137** | 1285 / **531** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 233 / **0** | 233 / **1** | 233 / **0** | 233 / **1** | 1500 / **1202** | 1500 / **552** | 1500 / **1049** | 1500 / **355** |
| T4 max BA vs observed p99 | 233 / **72** | 233 / **70** | 233 / **5** | 233 / **5** | 1500 / **1469** | 1500 / **1406** | 1500 / **1462** | 1500 / **1311** |
| T4 max BA vs observed p99.9 | 233 / **1** | 233 / **2** | 233 / **1** | 233 / **1** | 1500 / **1340** | 1500 / **750** | 1500 / **1196** | 1500 / **434** |
| T5 Eichhorn VOL~HT site independence | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 233 / **2** | 233 / **2** | 233 / **24** | 233 / **29** | 1500 / **39** | 1500 / **162** | 1500 / **51** | 1500 / **196** |
| T9 VOL inversions (>10% step drop) | 233 / **2** | 233 / **2** | 233 / **17** | 233 / **19** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 233 | 0.4 | 2.4 | 0.4 | 0.02 (201) | -1.71 | 38.1 vs 18.5 | 0.51 | 0.78 | yes | yes |
| default, GOMPIT | 233 | 0.4 | 2.4 | 2.1 | 0.08 (177) | -1.66 | 36.0 vs 18.5 | 0.54 | 0.78 | yes | yes |
| calibrated, native | 233 | 0.4 | 2.4 | 0.4 | 0.01 (70) | -1.74 | 29.2 vs 18.5 | 0.28 | 0.78 | yes | yes |
| calibrated, GOMPIT | 233 | 0.4 | 2.4 | 0.9 | 0.05 (73) | -1.80 | 27.8 vs 18.5 | 0.37 | 0.78 | yes | yes |
| CONUS b1, native | 1500 | 86.1 | 87.7 | 79.5 | 0.86 (1435) | -0.37 | 86.5 vs 18.5 | 0.78 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 48.8 | 51.3 | 43.4 | 0.46 (1348) | -1.26 | 59.8 vs 18.5 | 0.44 |  | no | no |
| CONUS b2, native | 1500 | 75.9 | 78.0 | 71.7 | 0.82 (1380) | -0.39 | 75.5 vs 18.5 | 0.78 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 28.7 | 31.1 | 29.9 | 0.41 (1285) | -1.31 | 53.3 vs 18.5 | 0.46 |  | no | no |

## Reading and component attribution (arm default, native)

- **Monotone height and QMD** (T1 monotone HT_dom): 29 of 233 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 72 of 233 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 4 of 233; T3 Reineke slope (trajectory OLS) 5 of 201; T4 max BA vs observed p99.9 1 of 233; T9 BA inversions (>10% step drop) 2 of 233; T9 VOL inversions (>10% step drop) 2 of 233.

Neither strong relation (Eichhorn, Reineke) is flagged on this arm among evaluated trajectories.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

