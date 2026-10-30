# NE: Bakuzis scorecard (catalog v6, fifth version (calibrated arms rerun after the overlay fix), October 30, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 51.3 and p99.9 62.7 m2 ha-1, SDI p99.9 1498 trees ha-1 from 3711 fitting plots. Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **RECOMMENDED**; recommended arm default_gompit; admissible: default_native; default_gompit; calibrated_native; calibrated_gompit.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **8** | 300 / **0** | 300 / **19** | 300 / **3** | 1500 / **1** | 1500 / **0** | 1500 / **9** |
| T1 monotone HT_dom | 300 / **9** | 300 / **53** | 300 / **11** | 300 / **55** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | 1080 / **51** | 1118 / **58** | 1072 / **629** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 159 / **10** | 161 / **12** | 74 / **4** | 85 / **7** | 1423 / **1286** | 1368 / **1167** | 1245 / **271** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **1** | 300 / **6** | 300 / **1** | 300 / **4** | 1500 / **1106** | 1500 / **838** | 1500 / **188** |
| T4 max BA vs observed p99 | 300 / **54** | 300 / **44** | 300 / **34** | 300 / **25** | 1500 / **1470** | 1500 / **1462** | 1500 / **1257** |
| T4 max BA vs observed p99.9 | 300 / **8** | 300 / **5** | 300 / **5** | 300 / **2** | 1500 / **1444** | 1500 / **1424** | 1500 / **835** |
| T5 Eichhorn VOL~HT site independence | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **0** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **1** | 1500 / **20** | 1500 / **12** | 1500 / **209** |
| T9 VOL inversions (>10% step drop) | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **2** | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.7 | 2.4 | 0.3 | 0.06 (159) | -1.71 | 42.6 vs 24.9 | 0.46 | 1.00 | yes | yes |
| default, GOMPIT | 300 | 1.3 | 3.4 | 3.3 | 0.07 (161) | -1.74 | 42.0 vs 24.9 | 0.54 | 1.00 | yes | yes |
| calibrated, native | 300 | 0.3 | 1.9 | 0.3 | 0.05 (74) | -1.68 | 36.5 vs 24.9 | 0.36 | 1.00 | yes | yes |
| calibrated, GOMPIT | 300 | 0.7 | 2.4 | 2.0 | 0.08 (85) | -1.78 | 36.3 vs 24.9 | 0.47 | 1.00 | yes | yes |
| CONUS b1, native | 1500 | 94.1 | 95.2 | 75.7 | 0.90 (1423) | -0.14 | 105.0 vs 24.9 | 0.87 |  | no | no |
| CONUS b2, native | 1500 | 91.4 | 92.7 | 61.6 | 0.85 (1368) | -0.17 | 89.5 vs 24.9 | 0.86 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 55.3 | 57.8 | 16.2 | 0.22 (1245) | -1.57 | 60.3 vs 24.9 | 0.47 |  | no | no |

## Reading and component attribution (arm default, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 53 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 44 of 300 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T3 Reineke slope (trajectory OLS) 12 of 161; T4 peak SDI / SDImax 6 of 300; T4 max BA vs observed p99.9 5 of 300; T9 BA inversions (>10% step drop) 1 of 300; T9 VOL inversions (>10% step drop) 1 of 300.

Neither strong relation (Eichhorn, Reineke) is flagged on this arm among evaluated trajectories.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

