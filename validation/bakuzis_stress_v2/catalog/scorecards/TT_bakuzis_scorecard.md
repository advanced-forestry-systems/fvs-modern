# TT: Bakuzis scorecard (catalog v3, October 29, 2026)

Engine arms: 142 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 57.2 and p99.9 72.6 m2 ha-1, SDI p99.9 1598 trees ha-1 from 142 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm conus_b2_gompit; admissible: none.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 142 / **1** | 142 / **1** | 142 / **1** | 142 / **1** | 1500 / **60** | 1500 / **16** | 1500 / **40** | 1500 / **11** |
| T1 monotone HT_dom | 142 / **9** | 142 / **39** | 142 / **16** | 142 / **37** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 142 / **0** | 142 / **2** | 142 / **0** | 142 / **2** | 683 / **261** | 481 / **258** | 705 / **355** | 350 / **222** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 31 / **26** | 7 / **3** | 19 / **5** | 12 / **2** | 296 / **101** | 175 / **50** | 294 / **59** | 161 / **36** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 142 / **0** | 142 / **0** | 142 / **0** | 142 / **0** | 1500 / **15** | 1500 / **2** | 1500 / **3** | 1500 / **2** |
| T4 max BA vs observed p99 | 142 / **29** | 142 / **11** | 142 / **4** | 142 / **3** | 1500 / **198** | 1500 / **80** | 1500 / **138** | 1500 / **44** |
| T4 max BA vs observed p99.9 | 142 / **1** | 142 / **1** | 142 / **1** | 142 / **1** | 1500 / **69** | 1500 / **26** | 1500 / **21** | 1500 / **7** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 142 / **1** | 142 / **2** | 142 / **2** | 142 / **3** | 1500 / **662** | 1500 / **889** | 1500 / **680** | 1500 / **979** |
| T9 VOL inversions (>10% step drop) | 142 / **0** | 142 / **0** | 142 / **2** | 142 / **2** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 142 | 0.7 | 3.9 | 0.7 | 0.84 (31) | -0.84 | 34.0 vs 16.8 | 0.91 | 1.00 | no | no |
| default, GOMPIT | 142 | 0.7 | 3.9 | 0.7 | not testable (7) | -1.65 | 21.5 vs 16.8 | 0.44 | 1.00 | no | no |
| calibrated, native | 142 | 0.7 | 3.9 | 0.7 | 0.26 (19) | -1.68 | 32.3 vs 16.8 | 0.85 | 1.00 | no | no |
| calibrated, GOMPIT | 142 | 0.7 | 3.9 | 0.7 | 0.17 (12) | -1.74 | 22.3 vs 16.8 | 0.56 | 1.00 | no | no |
| CONUS b1, native | 1500 | 3.7 | 4.8 | 1.4 | 0.34 (296) | -1.62 | 17.7 vs 16.8 | 0.47 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 1.1 | 1.7 | 0.3 | 0.29 (175) | -1.68 | 13.7 vs 16.8 | 0.37 |  | no | no |
| CONUS b2, native | 1500 | 1.3 | 2.0 | 0.7 | 0.20 (294) | -1.68 | 23.0 vs 16.8 | 0.50 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 0.3 | 0.8 | 0.3 | 0.22 (161) | -1.75 | 16.9 vs 16.8 | 0.37 |  | no | no |

## Reading and component attribution (arm default, native)

- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 26 of 31 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 29 of 142 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 1 of 142; T1 monotone HT_dom 9 of 142; T4 max BA vs observed p99.9 1 of 142; T9 BA inversions (>10% step drop) 1 of 142.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

