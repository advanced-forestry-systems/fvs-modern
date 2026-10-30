# LS: Bakuzis scorecard (catalog v4 (null-TPA fix), October 29, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 43.9 and p99.9 55.7 m2 ha-1, SDI p99.9 1425 trees ha-1 from 1066 fitting plots. Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **RECOMMENDED**; recommended arm calibrated_native; admissible: calibrated_native; calibrated_gompit.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **1** | 300 / **0** | 300 / **10** | 300 / **1** | 1500 / **0** | 1500 / **24** | 1500 / **0** | 1500 / **27** |
| T1 monotone HT_dom | 300 / **14** | 300 / **39** | 300 / **18** | 300 / **48** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **0** | 300 / **1** | 300 / **3** | 300 / **3** | 998 / **32** | 816 / **454** | 995 / **25** | 820 / **385** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 216 / **47** | 156 / **21** | 68 / **5** | 84 / **4** | 1157 / **1018** | 971 / **210** | 1120 / **932** | 889 / **150** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **12** | 300 / **12** | 300 / **0** | 300 / **3** | 1500 / **934** | 1500 / **235** | 1500 / **746** | 1500 / **149** |
| T4 max BA vs observed p99 | 300 / **197** | 300 / **133** | 300 / **51** | 300 / **37** | 1500 / **1206** | 1500 / **1063** | 1500 / **1301** | 1500 / **1065** |
| T4 max BA vs observed p99.9 | 300 / **73** | 300 / **41** | 300 / **16** | 300 / **4** | 1500 / **1183** | 1500 / **817** | 1500 / **1177** | 1500 / **652** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **0** | 300 / **4** | 300 / **6** | 300 / **11** | 1500 / **10** | 1500 / **289** | 1500 / **10** | 1500 / **353** |
| T9 VOL inversions (>10% step drop) | 300 / **2** | 300 / **5** | 300 / **4** | 300 / **10** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 14.7 | 19.1 | 3.0 | 0.22 (216) | -1.61 | 42.2 vs 16.8 | 0.69 | 1.00 | no | no |
| default, GOMPIT | 300 | 9.7 | 13.5 | 3.7 | 0.13 (156) | -1.64 | 37.3 vs 16.8 | 0.54 | 1.00 | no | no |
| calibrated, native | 300 | 3.3 | 6.0 | 0.0 | 0.07 (68) | -1.68 | 33.3 vs 16.8 | 0.45 | 1.00 | yes | yes |
| calibrated, GOMPIT | 300 | 1.3 | 3.4 | 1.0 | 0.05 (84) | -1.78 | 32.0 vs 16.8 | 0.48 | 1.00 | yes | yes |
| CONUS b1, native | 1500 | 77.1 | 79.2 | 61.9 | 0.88 (1157) | -0.14 | 92.9 vs 16.8 | 0.97 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 53.3 | 55.8 | 19.5 | 0.22 (971) | -1.56 | 51.8 vs 16.8 | 0.45 |  | no | no |
| CONUS b2, native | 1500 | 75.7 | 77.8 | 52.7 | 0.83 (1120) | -0.14 | 78.3 vs 16.8 | 0.97 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 42.2 | 44.7 | 12.5 | 0.17 (889) | -1.63 | 47.3 vs 16.8 | 0.47 |  | no | no |

## Reading and component attribution (arm calibrated, native)

- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 51 of 300 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 10 of 300; T1 monotone HT_dom 18 of 300; T1 monotone QMD 3 of 300; T3 Reineke slope (trajectory OLS) 5 of 68; T4 max BA vs observed p99.9 16 of 300; T9 BA inversions (>10% step drop) 6 of 300; T9 VOL inversions (>10% step drop) 4 of 300.

Strong relations (Eichhorn, Reineke) flagged: T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

