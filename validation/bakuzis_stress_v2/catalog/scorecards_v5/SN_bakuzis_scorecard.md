# SN: Bakuzis scorecard (catalog v4 (null-TPA fix), October 29, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 48.3 and p99.9 61.8 m2 ha-1, SDI p99.9 1443 trees ha-1 from 4682 fitting plots. Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **RECOMMENDED**; recommended arm default_native; admissible: default_native; default_gompit; calibrated_native; calibrated_gompit.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **7** | 300 / **2** | 300 / **22** | 300 / **8** | 1500 / **0** | 1500 / **0** | 1500 / **35** |
| T1 monotone HT_dom | 300 / **76** | 300 / **84** | 300 / **126** | 300 / **119** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **0** | 300 / **0** | 300 / **1** | 300 / **1** | 631 / **46** | 771 / **8** | 592 / **272** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 218 / **4** | 169 / **4** | 120 / **4** | 130 / **3** | 1446 / **1366** | 1334 / **1194** | 1168 / **517** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **0** | 300 / **1** | 300 / **0** | 300 / **1** | 1500 / **1363** | 1500 / **1185** | 1500 / **239** |
| T4 max BA vs observed p99 | 300 / **77** | 300 / **71** | 300 / **4** | 300 / **4** | 1500 / **1481** | 1500 / **1487** | 1500 / **1235** |
| T4 max BA vs observed p99.9 | 300 / **1** | 300 / **1** | 300 / **0** | 300 / **1** | 1500 / **1472** | 1500 / **1446** | 1500 / **757** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **0** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **1** | 300 / **1** | 300 / **1** | 300 / **1** | 1500 / **39** | 1500 / **11** | 1500 / **262** |
| T9 VOL inversions (>10% step drop) | 300 / **0** | 300 / **0** | 300 / **1** | 300 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.3 | 1.9 | 0.0 | 0.02 (218) | -1.73 | 43.2 vs 21.6 | 0.56 | 1.00 | yes | yes |
| default, GOMPIT | 300 | 0.3 | 1.9 | 0.3 | 0.02 (169) | -1.69 | 41.3 vs 21.6 | 0.51 | 1.00 | yes | yes |
| calibrated, native | 300 | 0.0 | 1.3 | 0.0 | 0.03 (120) | -1.73 | 37.2 vs 21.6 | 0.41 | 1.00 | yes | yes |
| calibrated, GOMPIT | 300 | 0.3 | 1.9 | 0.7 | 0.02 (130) | -1.79 | 37.2 vs 21.6 | 0.43 | 1.00 | yes | yes |
| CONUS b1, native | 1500 | 97.9 | 98.5 | 91.2 | 0.94 (1446) | 0.12 | 118.1 vs 21.6 | 1.07 |  | no | no |
| CONUS b2, native | 1500 | 95.5 | 96.4 | 81.2 | 0.90 (1334) | 0.11 | 94.9 vs 21.6 | 1.06 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 49.1 | 51.6 | 20.2 | 0.44 (1168) | -1.41 | 53.8 vs 21.6 | 0.48 |  | no | no |

## Reading and component attribution (arm default, native)

- **Monotone height and QMD** (T1 monotone HT_dom): 76 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 77 of 300 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 7 of 300; T3 Reineke slope (trajectory OLS) 4 of 218; T4 max BA vs observed p99.9 1 of 300; T9 BA inversions (>10% step drop) 1 of 300.

Strong relations (Eichhorn, Reineke) flagged: T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

