# UT: Bakuzis scorecard (catalog v4 (null-TPA fix), October 29, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 64.6 and p99.9 76.9 m2 ha-1, SDI p99.9 1818 trees ha-1 from 594 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm calibrated_gompit; admissible: none.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **0** | 300 / **0** | 300 / **1** | 300 / **1** | 1500 / **51** | 1500 / **16** | 1500 / **44** | 1500 / **16** |
| T1 monotone HT_dom | 300 / **37** | 300 / **92** | 300 / **38** | 300 / **89** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **0** | 300 / **5** | 300 / **0** | 300 / **5** | 791 / **281** | 782 / **211** | 836 / **336** | 730 / **217** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 67 / **38** | 37 / **11** | 59 / **17** | 38 / **10** | 147 / **47** | 119 / **62** | 127 / **35** | 111 / **46** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | 1500 / **4** | 1500 / **0** | 1500 / **0** | 1500 / **0** |
| T4 max BA vs observed p99 | 300 / **13** | 300 / **13** | 300 / **5** | 300 / **5** | 1500 / **84** | 1500 / **66** | 1500 / **56** | 1500 / **56** |
| T4 max BA vs observed p99.9 | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **1** | 1500 / **23** | 1500 / **23** | 1500 / **12** | 1500 / **10** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **0** | 300 / **1** | 300 / **0** | 300 / **2** | 1500 / **395** | 1500 / **450** | 1500 / **427** | 1500 / **499** |
| T9 VOL inversions (>10% step drop) | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.3 | 1.9 | 0.3 | 0.57 (67) | -1.12 | 38.0 vs 22.9 | 0.88 | 1.00 | no | no |
| default, GOMPIT | 300 | 0.3 | 1.9 | 0.3 | 0.30 (37) | -1.55 | 25.1 vs 22.9 | 0.47 | 1.00 | no | no |
| calibrated, native | 300 | 0.3 | 1.9 | 0.3 | 0.29 (59) | -1.55 | 35.7 vs 22.9 | 0.79 | 1.00 | no | no |
| calibrated, GOMPIT | 300 | 0.3 | 1.9 | 0.3 | 0.26 (38) | -1.59 | 27.4 vs 22.9 | 0.54 | 1.00 | no | no |
| CONUS b1, native | 1500 | 1.1 | 1.8 | 0.3 | 0.32 (147) | -1.59 | 23.8 vs 22.9 | 0.71 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 0.7 | 1.2 | 0.1 | 0.52 (119) | -1.54 | 25.6 vs 22.9 | 0.79 |  | no | no |
| CONUS b2, native | 1500 | 0.6 | 1.1 | 0.1 | 0.28 (127) | -1.62 | 24.9 vs 22.9 | 0.72 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 0.4 | 0.9 | 0.1 | 0.41 (111) | -1.63 | 25.8 vs 22.9 | 0.80 |  | no | no |

## Reading and component attribution (arm calibrated, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 89 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 10 of 38 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 1 of 300; T1 monotone QMD 5 of 300; T4 max BA vs observed p99 5 of 300; T4 max BA vs observed p99.9 1 of 300; T9 BA inversions (>10% step drop) 2 of 300.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

