# BM: Bakuzis scorecard (catalog v4 (null-TPA fix), October 29, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 55.3 and p99.9 68.0 m2 ha-1, SDI p99.9 1506 trees ha-1 from 474 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **RECOMMENDED**; recommended arm default_native; admissible: default_native; default_gompit; calibrated_native; calibrated_gompit; rule-sensitive (alternative Reineke gate admits a different set).

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **1** | 1500 / **18** | 1500 / **17** | 1500 / **16** | 1500 / **18** |
| T1 monotone HT_dom | 300 / **104** | 300 / **103** | 300 / **109** | 300 / **110** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **0** | 300 / **4** | 300 / **0** | 300 / **3** | 918 / **351** | 704 / **566** | 892 / **358** | 605 / **513** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 67 / **0** | 62 / **0** | 67 / **0** | 63 / **0** | 714 / **337** | 385 / **41** | 561 / **228** | 328 / **46** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | 1500 / **101** | 1500 / **12** | 1500 / **75** | 1500 / **11** |
| T4 max BA vs observed p99 | 300 / **22** | 300 / **18** | 300 / **28** | 300 / **22** | 1500 / **699** | 1500 / **220** | 1500 / **505** | 1500 / **166** |
| T4 max BA vs observed p99.9 | 300 / **2** | 300 / **2** | 300 / **2** | 300 / **2** | 1500 / **415** | 1500 / **57** | 1500 / **258** | 1500 / **48** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **1** | 1500 / **294** | 1500 / **965** | 1500 / **353** | 1500 / **1061** |
| T9 VOL inversions (>10% step drop) | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.7 | 2.4 | 0.3 | 0.00 (67) | -1.71 | 42.0 vs 17.9 | 0.62 | 1.00 | yes | yes |
| default, GOMPIT | 300 | 0.7 | 2.4 | 0.3 | 0.00 (62) | -1.71 | 42.0 vs 17.9 | 0.57 | 1.00 | yes | yes |
| calibrated, native | 300 | 0.7 | 2.4 | 0.3 | 0.00 (67) | -1.71 | 42.0 vs 17.9 | 0.63 | 1.00 | yes | yes |
| calibrated, GOMPIT | 300 | 0.7 | 2.4 | 0.3 | 0.00 (63) | -1.72 | 42.0 vs 17.9 | 0.57 | 1.00 | yes | yes |
| CONUS b1, native | 1500 | 25.5 | 27.8 | 8.9 | 0.47 (714) | -1.29 | 47.6 vs 17.9 | 0.78 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 3.8 | 4.9 | 1.5 | 0.11 (385) | -1.71 | 30.6 vs 17.9 | 0.52 |  | no | yes |
| CONUS b2, native | 1500 | 16.2 | 18.2 | 6.2 | 0.41 (561) | -1.52 | 40.1 vs 17.9 | 0.75 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 3.2 | 4.2 | 1.4 | 0.14 (328) | -1.77 | 26.7 vs 17.9 | 0.52 |  | no | yes |

## Reading and component attribution (arm default, native)

- **Monotone height and QMD** (T1 monotone HT_dom): 104 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 1 of 300; T4 max BA vs observed p99 22 of 300; T4 max BA vs observed p99.9 2 of 300; T9 BA inversions (>10% step drop) 1 of 300.

Strong relations (Eichhorn, Reineke) flagged: T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

