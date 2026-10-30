# PN: Bakuzis scorecard (catalog v4 (null-TPA fix), October 29, 2026)

Engine arms: 38 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 104.6 and p99.9 105.6 m2 ha-1, SDI p99.9 2230 trees ha-1 from 38 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm calibrated_gompit; admissible: none.

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 38 / **2** | 38 / **0** | 38 / **1** | 38 / **0** | 1500 / **116** | 1500 / **7** | 1500 / **99** | 1500 / **5** |
| T1 monotone HT_dom | 38 / **9** | 38 / **22** | 38 / **12** | 38 / **23** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 38 / **0** | 38 / **0** | 38 / **0** | 38 / **0** | 895 / **278** | 1055 / **757** | 916 / **302** | 1001 / **597** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 31 / **9** | 23 / **14** | 31 / **8** | 23 / **14** | 987 / **628** | 824 / **110** | 928 / **509** | 738 / **143** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 38 / **0** | 38 / **0** | 38 / **0** | 38 / **0** | 1500 / **209** | 1500 / **46** | 1500 / **69** | 1500 / **74** |
| T4 max BA vs observed p99 | 38 / **26** | 38 / **3** | 38 / **23** | 38 / **3** | 1500 / **882** | 1500 / **437** | 1500 / **723** | 1500 / **290** |
| T4 max BA vs observed p99.9 | 38 / **26** | 38 / **2** | 38 / **23** | 38 / **2** | 1500 / **873** | 1500 / **412** | 1500 / **706** | 1500 / **273** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **0** | 1 / **0** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 38 / **0** | 38 / **0** | 38 / **0** | 38 / **0** | 1500 / **451** | 1500 / **628** | 1500 / **489** | 1500 / **651** |
| T9 VOL inversions (>10% step drop) | 38 / **0** | 38 / **0** | 38 / **0** | 38 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 38 | 34.2 | 50.1 | 0.0 | 0.29 (31) | -1.38 | 100.2 vs 65.3 | 0.63 | 1.00 | no | no |
| default, GOMPIT | 38 | 2.6 | 13.5 | 0.0 | 0.61 (23) | -2.12 | 79.4 vs 65.3 | 0.57 | 1.00 | no | no |
| calibrated, native | 38 | 31.6 | 47.5 | 0.0 | 0.26 (31) | -1.41 | 100.9 vs 65.3 | 0.59 | 1.00 | no | no |
| calibrated, GOMPIT | 38 | 2.6 | 13.5 | 0.0 | 0.61 (23) | -2.10 | 79.9 vs 65.3 | 0.57 | 1.00 | no | no |
| CONUS b1, native | 1500 | 56.1 | 58.6 | 18.1 | 0.64 (987) | -0.86 | 109.8 vs 65.3 | 0.75 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 27.5 | 29.8 | 4.9 | 0.13 (824) | -1.78 | 73.1 vs 65.3 | 0.43 |  | no | no |
| CONUS b2, native | 1500 | 46.9 | 49.4 | 6.6 | 0.55 (928) | -1.26 | 95.7 vs 65.3 | 0.74 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 18.1 | 20.2 | 6.1 | 0.19 (738) | -1.80 | 66.3 vs 65.3 | 0.43 |  | no | no |

## Reading and component attribution (arm calibrated, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 23 of 38 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 14 of 23 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T4 max BA vs observed p99 3 of 38; T4 max BA vs observed p99.9 2 of 38.

Strong relations (Eichhorn, Reineke) flagged: T3. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

