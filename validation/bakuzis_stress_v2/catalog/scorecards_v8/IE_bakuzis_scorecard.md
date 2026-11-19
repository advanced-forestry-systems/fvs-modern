# IE: Bakuzis scorecard (catalog v8, sixth version (D5 plot-size fix, all arms rerun), October 30, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 60.7 and p99.9 71.4 m2 ha-1, SDI p99.9 1623 trees ha-1 from 1051 fitting plots. Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm calibrated_gompit; admissible: none.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **0** | 300 / **0** | 300 / **19** | 300 / **0** | 1500 / **57** | 1500 / **19** | 1500 / **48** | 1500 / **19** |
| T1 monotone HT_dom | 300 / **50** | 300 / **96** | 300 / **77** | 300 / **117** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **8** | 300 / **26** | 300 / **13** | 300 / **27** | 962 / **290** | 853 / **581** | 924 / **271** | 805 / **541** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 44 / **14** | 41 / **23** | 47 / **22** | 38 / **21** | 837 / **445** | 631 / **75** | 744 / **346** | 555 / **80** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | 1500 / **204** | 1500 / **57** | 1500 / **118** | 1500 / **34** |
| T4 max BA vs observed p99 | 300 / **12** | 300 / **17** | 300 / **17** | 300 / **17** | 1500 / **774** | 1500 / **408** | 1500 / **611** | 1500 / **290** |
| T4 max BA vs observed p99.9 | 300 / **1** | 300 / **4** | 300 / **4** | 300 / **6** | 1500 / **591** | 1500 / **217** | 1500 / **407** | 1500 / **132** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **0** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **2** | 300 / **22** | 300 / **19** | 300 / **25** | 1500 / **297** | 1500 / **620** | 1500 / **302** | 1500 / **764** |
| T9 VOL inversions (>10% step drop) | 300 / **0** | 300 / **8** | 300 / **3** | 300 / **8** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.3 | 1.9 | 0.0 | 0.32 (44) | -1.91 | 28.8 vs 21.7 | 0.52 | 1.00 | no | no |
| default, GOMPIT | 300 | 0.3 | 1.9 | 0.0 | 0.56 (41) | -1.60 | 22.4 vs 21.7 | 0.49 | 1.00 | no | no |
| calibrated, native | 300 | 0.3 | 1.9 | 0.0 | 0.47 (47) | -1.77 | 26.9 vs 21.7 | 0.51 | 1.00 | no | no |
| calibrated, GOMPIT | 300 | 0.3 | 1.9 | 0.0 | 0.55 (38) | -1.65 | 22.2 vs 21.7 | 0.48 | 1.00 | no | no |
| CONUS b1, native | 1500 | 36.9 | 39.4 | 15.9 | 0.53 (837) | -1.18 | 57.1 vs 21.7 | 0.65 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 14.5 | 16.3 | 4.9 | 0.12 (629) | -1.72 | 38.4 vs 21.7 | 0.43 |  | no | no |
| CONUS b2, native | 1500 | 26.3 | 28.6 | 10.1 | 0.47 (743) | -1.44 | 48.5 vs 21.7 | 0.65 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 8.7 | 10.3 | 3.3 | 0.14 (555) | -1.76 | 35.1 vs 21.7 | 0.43 |  | no | no |

## Reading and component attribution (arm calibrated, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 117 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 21 of 38 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T1 monotone QMD 27 of 300; T4 max BA vs observed p99 17 of 300; T4 max BA vs observed p99.9 6 of 300; T9 BA inversions (>10% step drop) 25 of 300; T9 VOL inversions (>10% step drop) 8 of 300.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

