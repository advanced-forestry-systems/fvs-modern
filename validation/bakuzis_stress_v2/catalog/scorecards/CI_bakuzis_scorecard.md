# CI: Bakuzis scorecard (catalog v3, October 29, 2026)

Engine arms: 297 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 56.5 and p99.9 61.3 m2 ha-1, SDI p99.9 1313 trees ha-1 from 425 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm calibrated_native; admissible: none; rule-sensitive (alternative Reineke gate admits a different set).

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 297 / **0** | 297 / **0** | 297 / **0** | 297 / **0** | 1500 / **83** | 1500 / **14** | 1500 / **63** | 1500 / **14** |
| T1 monotone HT_dom | 297 / **18** | 297 / **74** | 297 / **22** | 297 / **74** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 297 / **6** | 297 / **12** | 297 / **6** | 297 / **12** | 739 / **347** | 565 / **470** | 721 / **401** | 475 / **405** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 79 / **8** | 74 / **24** | 58 / **10** | 74 / **23** | 488 / **225** | 336 / **55** | 418 / **149** | 319 / **57** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 297 / **0** | 297 / **1** | 297 / **0** | 297 / **1** | 1500 / **71** | 1500 / **10** | 1500 / **34** | 1500 / **10** |
| T4 max BA vs observed p99 | 297 / **6** | 297 / **17** | 297 / **5** | 297 / **17** | 1500 / **276** | 1500 / **73** | 1500 / **162** | 1500 / **55** |
| T4 max BA vs observed p99.9 | 297 / **1** | 297 / **11** | 297 / **1** | 297 / **11** | 1500 / **228** | 1500 / **45** | 1500 / **102** | 1500 / **35** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 297 / **0** | 297 / **12** | 297 / **0** | 297 / **12** | 1500 / **564** | 1500 / **977** | 1500 / **659** | 1500 / **1032** |
| T9 VOL inversions (>10% step drop) | 297 / **0** | 297 / **4** | 297 / **0** | 297 / **4** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 297 | 0.3 | 1.9 | 0.3 | 0.10 (79) | -1.67 | 31.9 vs 20.9 | 0.54 | 0.99 | no | yes |
| default, GOMPIT | 297 | 1.7 | 3.9 | 1.0 | 0.32 (74) | -1.58 | 25.9 vs 20.9 | 0.48 | 0.99 | no | no |
| calibrated, native | 297 | 0.3 | 1.9 | 0.3 | 0.17 (58) | -1.82 | 30.3 vs 20.9 | 0.51 | 0.99 | no | no |
| calibrated, GOMPIT | 297 | 1.7 | 3.9 | 1.0 | 0.31 (74) | -1.59 | 25.9 vs 20.9 | 0.48 | 0.99 | no | no |
| CONUS b1, native | 1500 | 12.8 | 14.6 | 5.5 | 0.46 (488) | -1.41 | 23.6 vs 20.9 | 0.64 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 3.0 | 4.0 | 1.1 | 0.16 (336) | -1.69 | 20.2 vs 20.9 | 0.52 |  | no | no |
| CONUS b2, native | 1500 | 5.9 | 7.2 | 2.7 | 0.36 (418) | -1.61 | 22.6 vs 20.9 | 0.64 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 2.3 | 3.2 | 1.1 | 0.18 (319) | -1.70 | 20.0 vs 20.9 | 0.52 |  | no | no |

## Reading and component attribution (arm calibrated, native)

- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 10 of 58 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T1 monotone HT_dom 22 of 297; T1 monotone QMD 6 of 297; T4 max BA vs observed p99 5 of 297; T4 max BA vs observed p99.9 1 of 297.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

