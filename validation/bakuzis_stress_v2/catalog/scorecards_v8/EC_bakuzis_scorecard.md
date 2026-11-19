# EC: Bakuzis scorecard (catalog v8, sixth version (D5 plot-size fix, all arms rerun), October 30, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 83.4 and p99.9 96.4 m2 ha-1, SDI p99.9 1789 trees ha-1 from 369 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm calibrated_gompit; admissible: none; rule-sensitive (alternative Reineke gate admits a different set).

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **30** | 300 / **5** | 300 / **1** | 300 / **1** | 1500 / **47** | 1500 / **18** | 1500 / **39** | 1500 / **18** |
| T1 monotone HT_dom | 300 / **123** | 300 / **149** | 300 / **39** | 300 / **105** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **1** | 300 / **18** | 300 / **1** | 300 / **18** | 1044 / **392** | 905 / **700** | 1017 / **416** | 813 / **633** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 29 / **5** | 24 / **8** | 172 / **100** | 61 / **16** | 718 / **234** | 498 / **85** | 577 / **174** | 411 / **103** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **1** | 1500 / **80** | 1500 / **27** | 1500 / **26** | 1500 / **18** |
| T4 max BA vs observed p99 | 300 / **4** | 300 / **4** | 300 / **21** | 300 / **14** | 1500 / **371** | 1500 / **97** | 1500 / **152** | 1500 / **47** |
| T4 max BA vs observed p99.9 | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **1** | 1500 / **176** | 1500 / **35** | 1500 / **42** | 1500 / **13** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **19** | 300 / **26** | 300 / **1** | 300 / **9** | 1500 / **427** | 1500 / **960** | 1500 / **484** | 1500 / **1073** |
| T9 VOL inversions (>10% step drop) | 300 / **11** | 300 / **23** | 300 / **1** | 300 / **5** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.3 | 1.9 | 0.3 | 0.17 (29) | -1.86 | 24.1 vs 21.2 | 0.39 | 1.00 | no | no |
| default, GOMPIT | 300 | 0.3 | 1.9 | 0.3 | 0.33 (24) | -2.04 | 22.5 vs 21.2 | 0.43 | 1.00 | no | no |
| calibrated, native | 300 | 0.3 | 1.9 | 0.3 | 0.58 (172) | -1.09 | 49.4 vs 21.2 | 0.87 | 1.00 | no | no |
| calibrated, GOMPIT | 300 | 0.3 | 1.9 | 0.3 | 0.26 (61) | -1.56 | 31.7 vs 21.2 | 0.61 | 1.00 | no | no |
| CONUS b1, native | 1500 | 10.1 | 11.7 | 6.8 | 0.33 (718) | -1.59 | 55.3 vs 21.2 | 0.59 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 2.3 | 3.2 | 2.3 | 0.17 (498) | -1.77 | 38.5 vs 21.2 | 0.42 |  | no | yes |
| CONUS b2, native | 1500 | 2.7 | 3.6 | 2.8 | 0.30 (577) | -1.68 | 45.1 vs 21.2 | 0.57 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 0.9 | 1.5 | 1.7 | 0.25 (411) | -1.80 | 34.4 vs 21.2 | 0.42 |  | no | no |

## Reading and component attribution (arm calibrated, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 105 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 16 of 61 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 1 of 300; T1 monotone QMD 18 of 300; T4 peak SDI / SDImax 1 of 300; T4 max BA vs observed p99 14 of 300; T4 max BA vs observed p99.9 1 of 300; T9 BA inversions (>10% step drop) 9 of 300; T9 VOL inversions (>10% step drop) 5 of 300.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

