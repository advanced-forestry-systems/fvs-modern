# CA: Bakuzis scorecard (catalog v3, October 29, 2026)

Engine arms: 142 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 81.6 and p99.9 92.3 m2 ha-1, SDI p99.9 1838 trees ha-1 from 142 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm calibrated_gompit; admissible: none; rule-sensitive (alternative Reineke gate admits a different set).

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 142 / **0** | 142 / **0** | 142 / **2** | 142 / **0** | 1500 / **42** | 1500 / **75** | 1500 / **43** | 1500 / **79** |
| T1 monotone HT_dom | 142 / **1** | 142 / **17** | 142 / **13** | 142 / **25** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 142 / **0** | 142 / **3** | 142 / **1** | 142 / **4** | 1194 / **257** | 1071 / **539** | 1175 / **296** | 1061 / **595** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 89 / **69** | 34 / **19** | 46 / **20** | 22 / **9** | 730 / **206** | 523 / **84** | 552 / **192** | 422 / **98** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 142 / **1** | 142 / **1** | 142 / **2** | 142 / **1** | 1500 / **73** | 1500 / **20** | 1500 / **21** | 1500 / **8** |
| T4 max BA vs observed p99 | 142 / **52** | 142 / **10** | 142 / **12** | 142 / **3** | 1500 / **435** | 1500 / **166** | 1500 / **218** | 1500 / **62** |
| T4 max BA vs observed p99.9 | 142 / **1** | 142 / **1** | 142 / **7** | 142 / **2** | 1500 / **301** | 1500 / **64** | 1500 / **72** | 1500 / **24** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 142 / **0** | 142 / **0** | 142 / **2** | 142 / **3** | 1500 / **545** | 1500 / **988** | 1500 / **700** | 1500 / **1125** |
| T9 VOL inversions (>10% step drop) | 142 / **0** | 142 / **0** | 142 / **1** | 142 / **2** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 142 | 0.7 | 3.9 | 0.7 | 0.78 (89) | -0.83 | 54.4 vs 27.4 | 0.86 | 1.00 | no | no |
| default, GOMPIT | 142 | 0.7 | 3.9 | 0.7 | 0.56 (34) | -1.45 | 38.9 vs 27.4 | 0.64 | 1.00 | no | no |
| calibrated, native | 142 | 0.7 | 3.9 | 1.4 | 0.43 (46) | -1.56 | 47.3 vs 27.4 | 0.76 | 1.00 | no | no |
| calibrated, GOMPIT | 142 | 0.7 | 3.9 | 0.7 | 0.41 (22) | -1.83 | 37.8 vs 27.4 | 0.58 | 1.00 | no | no |
| CONUS b1, native | 1500 | 19.1 | 21.2 | 6.7 | 0.28 (730) | -1.62 | 56.4 vs 27.4 | 0.50 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 4.3 | 5.4 | 2.1 | 0.16 (523) | -1.78 | 42.2 vs 27.4 | 0.33 |  | no | yes |
| CONUS b2, native | 1500 | 4.7 | 5.9 | 2.4 | 0.35 (552) | -1.67 | 41.3 vs 27.4 | 0.47 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 1.6 | 2.4 | 0.9 | 0.23 (422) | -1.83 | 35.4 vs 27.4 | 0.34 |  | no | no |

## Reading and component attribution (arm calibrated, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 25 of 142 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 9 of 22 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T1 monotone QMD 4 of 142; T4 peak SDI / SDImax 1 of 142; T4 max BA vs observed p99 3 of 142; T4 max BA vs observed p99.9 2 of 142; T9 BA inversions (>10% step drop) 3 of 142; T9 VOL inversions (>10% step drop) 2 of 142.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

