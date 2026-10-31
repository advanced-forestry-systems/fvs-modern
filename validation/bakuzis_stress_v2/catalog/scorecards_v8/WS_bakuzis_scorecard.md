# WS: Bakuzis scorecard (catalog v8, sixth version (D5 plot-size fix, all arms rerun), October 30, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 86.2 and p99.9 107.2 m2 ha-1, SDI p99.9 1991 trees ha-1 from 605 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **RECOMMENDED**; recommended arm calibrated_native; admissible: calibrated_native; rule-sensitive (alternative Reineke gate admits a different set).

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **9** | 300 / **5** | 300 / **19** | 300 / **7** | 1500 / **25** | 1500 / **19** | 1500 / **21** | 1500 / **28** |
| T1 monotone HT_dom | 300 / **19** | 300 / **59** | 300 / **34** | 300 / **90** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **1** | 300 / **25** | 300 / **2** | 300 / **29** | 885 / **369** | 780 / **502** | 835 / **361** | 730 / **515** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 99 / **16** | 38 / **18** | 80 / **5** | 38 / **11** | 382 / **130** | 304 / **80** | 303 / **110** | 258 / **94** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **2** | 300 / **1** | 300 / **2** | 300 / **1** | 1500 / **22** | 1500 / **11** | 1500 / **9** | 1500 / **6** |
| T4 max BA vs observed p99 | 300 / **43** | 300 / **5** | 300 / **44** | 300 / **6** | 1500 / **186** | 1500 / **105** | 1500 / **89** | 1500 / **62** |
| T4 max BA vs observed p99.9 | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **1** | 1500 / **50** | 1500 / **13** | 1500 / **20** | 1500 / **18** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **8** | 300 / **17** | 300 / **18** | 300 / **26** | 1500 / **740** | 1500 / **951** | 1500 / **837** | 1500 / **1095** |
| T9 VOL inversions (>10% step drop) | 300 / **5** | 300 / **10** | 300 / **12** | 300 / **17** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.3 | 1.9 | 1.0 | 0.16 (99) | -1.59 | 46.5 vs 30.0 | 0.86 | 1.00 | no | yes |
| default, GOMPIT | 300 | 0.3 | 1.9 | 0.3 | 0.47 (38) | -2.15 | 37.9 vs 30.0 | 0.54 | 1.00 | no | no |
| calibrated, native | 300 | 0.3 | 1.9 | 0.7 | 0.06 (80) | -1.65 | 44.3 vs 30.0 | 0.74 | 1.00 | yes | yes |
| calibrated, GOMPIT | 300 | 0.3 | 1.9 | 0.3 | 0.29 (38) | -1.92 | 35.9 vs 30.0 | 0.51 | 1.00 | no | no |
| CONUS b1, native | 1500 | 3.3 | 4.4 | 1.9 | 0.34 (382) | -1.69 | 31.9 vs 30.0 | 0.52 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 0.9 | 1.5 | 0.9 | 0.26 (304) | -1.82 | 30.7 vs 30.0 | 0.45 |  | no | no |
| CONUS b2, native | 1500 | 1.3 | 2.1 | 0.9 | 0.36 (303) | -1.81 | 26.3 vs 30.0 | 0.49 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 1.2 | 1.9 | 0.7 | 0.36 (258) | -1.85 | 28.9 vs 30.0 | 0.45 |  | no | no |

## Reading and component attribution (arm calibrated, native)

- **Monotone height and QMD** (T1 monotone HT_dom): 34 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 44 of 300 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 19 of 300; T1 monotone QMD 2 of 300; T3 Reineke slope (trajectory OLS) 5 of 80; T4 peak SDI / SDImax 2 of 300; T4 max BA vs observed p99.9 1 of 300; T9 BA inversions (>10% step drop) 18 of 300; T9 VOL inversions (>10% step drop) 12 of 300.

Strong relations (Eichhorn, Reineke) flagged: T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

