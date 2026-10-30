# CR: Bakuzis scorecard (catalog v6, fifth version (calibrated arms rerun after the overlay fix), October 30, 2026)

Engine arms: 300 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 64.8 and p99.9 77.2 m2 ha-1, SDI p99.9 1714 trees ha-1 from 2370 fitting plots. Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm calibrated_gompit; admissible: none; rule-sensitive (alternative Reineke gate admits a different set).

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 300 / **1** | 300 / **0** | 300 / **6** | 300 / **1** | 1500 / **67** | 1500 / **10** | 1500 / **58** | 1500 / **10** |
| T1 monotone HT_dom | 300 / **79** | 300 / **218** | 300 / **98** | 300 / **224** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 300 / **1** | 300 / **2** | 300 / **1** | 300 / **9** | 774 / **277** | 875 / **487** | 770 / **270** | 867 / **430** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 129 / **35** | 33 / **10** | 72 / **9** | 24 / **5** | 277 / **49** | 165 / **25** | 199 / **42** | 137 / **34** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 300 / **0** | 300 / **2** | 300 / **0** | 300 / **0** | 1500 / **13** | 1500 / **1** | 1500 / **4** | 1500 / **1** |
| T4 max BA vs observed p99 | 300 / **36** | 300 / **10** | 300 / **20** | 300 / **2** | 1500 / **185** | 1500 / **59** | 1500 / **108** | 1500 / **32** |
| T4 max BA vs observed p99.9 | 300 / **0** | 300 / **0** | 300 / **0** | 300 / **0** | 1500 / **84** | 1500 / **16** | 1500 / **32** | 1500 / **7** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 300 / **2** | 300 / **9** | 300 / **8** | 300 / **15** | 1500 / **740** | 1500 / **690** | 1500 / **749** | 1500 / **663** |
| T9 VOL inversions (>10% step drop) | 300 / **2** | 300 / **10** | 300 / **6** | 300 / **18** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 300 | 0.0 | 1.3 | 0.0 | 0.27 (129) | -1.46 | 42.7 vs 24.0 | 0.89 | 1.00 | no | no |
| default, GOMPIT | 300 | 0.0 | 1.3 | 0.3 | 0.30 (33) | -1.75 | 26.4 vs 24.0 | 0.56 | 1.00 | no | no |
| calibrated, native | 300 | 0.0 | 1.3 | 0.0 | 0.12 (72) | -1.58 | 37.9 vs 24.0 | 0.81 | 1.00 | no | yes |
| calibrated, GOMPIT | 300 | 0.0 | 1.3 | 0.0 | 0.21 (24) | -1.89 | 26.9 vs 24.0 | 0.55 | 1.00 | no | no |
| CONUS b1, native | 1500 | 5.3 | 6.6 | 0.9 | 0.18 (277) | -1.74 | 20.5 vs 24.0 | 0.51 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 1.1 | 1.7 | 0.2 | 0.15 (165) | -1.70 | 26.3 vs 24.0 | 0.52 |  | no | yes |
| CONUS b2, native | 1500 | 1.9 | 2.8 | 0.3 | 0.21 (199) | -1.75 | 19.6 vs 24.0 | 0.52 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 0.4 | 0.9 | 0.1 | 0.25 (137) | -1.72 | 25.5 vs 24.0 | 0.53 |  | no | no |

## Reading and component attribution (arm calibrated, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 224 of 300 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 5 of 24 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 1 of 300; T1 monotone QMD 9 of 300; T4 max BA vs observed p99 2 of 300; T9 BA inversions (>10% step drop) 15 of 300; T9 VOL inversions (>10% step drop) 18 of 300.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

