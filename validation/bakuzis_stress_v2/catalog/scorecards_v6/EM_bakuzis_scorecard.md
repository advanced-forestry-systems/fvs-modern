# EM: Bakuzis scorecard (catalog v6, fifth version (calibrated arms rerun after the overlay fix), October 30, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 58.9 and p99.9 67.8 m2 ha-1, SDI p99.9 1687 trees ha-1 from 865 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm conus_b2_gompit; admissible: none.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | 1500 / **96** | 1500 / **24** | 1500 / **84** | 1500 / **22** |
| T1 monotone HT_dom | 300 / **59** | 300 / **110** | 300 / **57** | 300 / **103** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **9** | 300 / **9** | 300 / **12** | 300 / **8** | 800 / **379** | 633 / **451** | 798 / **421** | 562 / **384** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 202 / **130** | 82 / **11** | 191 / **117** | 80 / **11** | 356 / **109** | 239 / **63** | 303 / **93** | 235 / **70** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | 1500 / **17** | 1500 / **0** | 1500 / **4** | 1500 / **1** |
| T4 max BA vs observed p99 | 300 / **206** | 300 / **47** | 300 / **190** | 300 / **47** | 1500 / **273** | 1500 / **99** | 1500 / **163** | 1500 / **72** |
| T4 max BA vs observed p99.9 | 300 / **168** | 300 / **33** | 300 / **146** | 300 / **36** | 1500 / **170** | 1500 / **34** | 1500 / **78** | 1500 / **24** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **0** | 1 / **1** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **2** | 300 / **16** | 300 / **3** | 300 / **16** | 1500 / **543** | 1500 / **837** | 1500 / **587** | 1500 / **860** |
| T9 VOL inversions (>10% step drop) | 300 / **5** | 300 / **9** | 300 / **7** | 300 / **10** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 22.3 | 27.4 | 0.0 | 0.64 (202) | -1.04 | 53.1 vs 23.7 | 0.84 | 1.00 | no | no |
| default, GOMPIT | 300 | 7.7 | 11.2 | 0.0 | 0.13 (82) | -1.65 | 30.9 vs 23.7 | 0.50 | 1.00 | no | no |
| calibrated, native | 300 | 22.3 | 27.4 | 0.0 | 0.61 (191) | -1.03 | 52.0 vs 23.7 | 0.81 | 1.00 | no | no |
| calibrated, GOMPIT | 300 | 9.7 | 13.5 | 0.0 | 0.14 (80) | -1.68 | 31.2 vs 23.7 | 0.47 | 1.00 | no | no |
| CONUS b1, native | 1500 | 9.8 | 11.4 | 1.9 | 0.31 (356) | -1.59 | 25.5 vs 23.7 | 0.57 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 2.1 | 3.0 | 0.3 | 0.26 (239) | -1.73 | 21.1 vs 23.7 | 0.50 |  | no | no |
| CONUS b2, native | 1500 | 4.3 | 5.5 | 0.3 | 0.31 (303) | -1.69 | 27.1 vs 23.7 | 0.57 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 1.3 | 2.1 | 0.1 | 0.30 (235) | -1.74 | 23.3 vs 23.7 | 0.49 |  | no | no |

## Reading and component attribution (arm default, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 110 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 11 of 82 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 47 of 300 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99.9): 33 of 300 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.

Minor (under 10 percent of evaluated trajectories): T1 monotone QMD 9 of 300; T9 BA inversions (>10% step drop) 16 of 300; T9 VOL inversions (>10% step drop) 9 of 300.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

