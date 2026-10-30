# WC: Bakuzis scorecard (catalog v6, fifth version (calibrated arms rerun after the overlay fix), October 30, 2026)

Engine arms: 226 seed-passing FIA stands, 100 yr at a 5 yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 98.5 and p99.9 109.4 m2 ha-1, SDI p99.9 2087 trees ha-1 from 226 fitting plots (fewer than 1,000 plots, so the p99.9 is the sample maximum). Reineke screen qmd_rate 0.10 cm yr-1, relative-density gate 0.55.

Catalog verdict (primary rule): **NOT RECOMMENDED for long horizons (least-bad arm named)**; recommended arm default_gompit; admissible: none; rule-sensitive (alternative Reineke gate admits a different set).

How to read the site-based relations. The engine arms run real FIA stands, and the site and density classes are terciles of the condition site index (SICOND, which is species-specific) and live basal area, so the classes differ in species, age and structure as well as in site and density. Site ordering (T2), Eichhorn (T5), Sukachev (T6) and yield-density (T7) are therefore confounded by the stand selection and flag in nearly every variant; they are reported, not gated, and the synthetic 36-cell matrix (memo Section 7) is the clean test of those relations. The engine verdicts rest on the Reineke axis (T3, T3b), the magnitude envelope (T4), the plateau (T8), inversions (T9) and the initial-step shock (T0). Top height (T1 HT_dom) can fall when the tallest trees die and is a minor flag on GOMPIT arms, not a growth defect. A pass on any relation is a failure to falsify, not a confirmation of accuracy.

## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)

| Relation | default, native | default, GOMPIT | calibrated, native | calibrated, GOMPIT | CONUS b1, native | CONUS b1, GOMPIT | CONUS b2, native | CONUS b2, GOMPIT |
|---|---|---|---|---|---|---|---|---|
| T0 initial-step shock | 226 / **4** | 226 / **1** | 226 / **5** | 226 / **1** | 1500 / **19** | 1500 / **25** | 1500 / **16** | 1500 / **21** |
| T1 monotone HT_dom | 226 / **62** | 226 / **99** | 226 / **54** | 226 / **99** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T1 monotone QMD | 226 / **1** | 226 / **6** | 226 / **2** | 226 / **7** | 1152 / **522** | 1062 / **685** | 1107 / **571** | 1002 / **678** |
| T2 site ordering HT | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T2 site ordering VOL | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T3 Reineke slope (trajectory OLS) | 165 / **100** | 117 / **20** | 168 / **110** | 118 / **20** | 935 / **342** | 892 / **234** | 802 / **332** | 732 / **264** |
| T3b frontier slope (rq tau 0.95, pooled limb) | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T4 peak SDI / SDImax | 226 / **1** | 226 / **0** | 226 / **4** | 226 / **0** | 1500 / **145** | 1500 / **67** | 1500 / **63** | 1500 / **38** |
| T4 max BA vs observed p99 | 226 / **117** | 226 / **35** | 226 / **123** | 226 / **39** | 1500 / **552** | 1500 / **313** | 1500 / **285** | 1500 / **163** |
| T4 max BA vs observed p99.9 | 226 / **88** | 226 / **17** | 226 / **94** | 226 / **21** | 1500 / **352** | 1500 / **146** | 1500 / **142** | 1500 / **73** |
| T5 Eichhorn VOL~HT site independence | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T6 Sukachev per-tree VOL falls with density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T7 yield-density final BA vs initial density | 1 / **1** | 1 / **1** | 1 / **1** | 1 / **1** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |
| T8 pinned BA plateau | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** | 1 / **0** |
| T9 BA inversions (>10% step drop) | 226 / **1** | 226 / **2** | 226 / **1** | 226 / **2** | 1500 / **637** | 1500 / **929** | 1500 / **754** | 1500 / **1012** |
| T9 VOL inversions (>10% step drop) | 226 / **1** | 226 / **1** | 226 / **0** | 226 / **0** | NE (0 NA) | NE (0 NA) | NE (0 NA) | NE (0 NA) |

## Magnitude and gates (engine and projector arms)

| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| default, native | 226 | 10.6 | 15.3 | 0.4 | 0.61 (165) | -1.11 | 77.0 vs 50.4 | 0.65 | 1.00 | no | no |
| default, GOMPIT | 226 | 3.1 | 6.3 | 0.4 | 0.17 (117) | -1.58 | 61.8 vs 50.4 | 0.60 | 1.00 | no | yes |
| calibrated, native | 226 | 18.1 | 23.7 | 4.9 | 0.65 (168) | -1.09 | 78.6 vs 50.4 | 0.69 | 1.00 | no | no |
| calibrated, GOMPIT | 226 | 4.0 | 7.4 | 2.2 | 0.17 (118) | -1.57 | 61.5 vs 50.4 | 0.61 | 1.00 | no | yes |
| CONUS b1, native | 1500 | 23.0 | 25.2 | 12.3 | 0.37 (935) | -1.61 | 75.8 vs 50.4 | 0.54 |  | no | no |
| CONUS b1, GOMPIT | 1500 | 9.7 | 11.3 | 6.3 | 0.26 (892) | -1.82 | 63.8 vs 50.4 | 0.44 |  | no | no |
| CONUS b2, native | 1500 | 9.5 | 11.1 | 5.9 | 0.41 (802) | -1.72 | 60.0 vs 50.4 | 0.50 |  | no | no |
| CONUS b2, GOMPIT | 1500 | 4.9 | 6.1 | 3.8 | 0.36 (732) | -1.92 | 53.7 vs 50.4 | 0.44 |  | no | no |

## Reading and component attribution (arm default, GOMPIT)

- **Monotone height and QMD** (T1 monotone HT_dom): 99 of 226 evaluated flagged. diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question.
- **Site ordering (height, volume)** (T2 site ordering HT): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Site ordering (height, volume)** (T2 site ordering VOL): 1 of 1 evaluated flagged. height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one.
- **Reineke slope, within stand (RD gated)** (T3 Reineke slope (trajectory OLS)): 20 of 117 evaluated flagged. mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify.
- **Magnitude vs FIA envelope (BA p99, p99.9; SDI)** (T4 max BA vs observed p99): 35 of 226 evaluated flagged. the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it.
- **Eichhorn (volume at matched height across sites)** (T5 Eichhorn VOL~HT site independence): 1 of 1 evaluated flagged. the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site.
- **Sukachev (per-tree volume vs density)** (T6 Sukachev per-tree VOL falls with density): 1 of 1 evaluated flagged. mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7.
- **Yield-density (final BA vs initial density)** (T7 yield-density final BA vs initial density): 1 of 1 evaluated flagged. self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier.

Minor (under 10 percent of evaluated trajectories): T0 initial-step shock 1 of 226; T1 monotone QMD 6 of 226; T4 max BA vs observed p99.9 17 of 226; T9 BA inversions (>10% step drop) 2 of 226; T9 VOL inversions (>10% step drop) 1 of 226.

Strong relations (Eichhorn, Reineke) flagged: T3, T5. These outrank the others; start there.

Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.

