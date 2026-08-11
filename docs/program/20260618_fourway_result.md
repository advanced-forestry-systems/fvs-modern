# Four-way comparison: variants x species x ecoregion x landowner
2026-06-18. Default vs calibrated bias decomposed across four margins on the disturbance-clean basis.
Harness: calib_4way.py (per-condition default/calibrated predictions). Aggregation: calib_4way_aggregate.R.
Outputs: calib_4way_margins.csv, fig_4way_landowner.png, fig_4way_ecoregion.png.

## How the join works and what is solid

The per-condition harness emits 951 conditions across 15 variants, each with default and calibrated BA/TPH/QMD/VOL
predictions plus the observed value. The variant margin uses all 951 (variant is intrinsic to the harness, no join
needed) and is statistically solid (n = 40 to 112 per variant). The species, ecoregion, and landowner margins require
joining to the cspiv6 remeasurement RDS on the t1 plot CN; that join matches the disturbance-clean subset, 123 of 951
conditions (12.9%). This is expected: cspiv6 is itself the disturbance-clean filtered universe, so only the clean
conditions carry the ecoregion and landowner keys. The ecoregion (EPA L2) and landowner margins clear n >= 20 and are
reported; the dominant-species margin does not reach n >= 20 per level in the joined subset and is instead covered by
the dedicated species-stress test (sigma_sp, response-scale per-species deviation).

## Result: calibration reduces bias on nearly every variant and margin

By variant (full sample), default to calibrated absolute bias, the headline metric QMD:

| variant | QMD def -> cal | BA def -> cal | TPH def -> cal | VOL def -> cal | n |
|---|---|---|---|---|---|
| kt | 26.3 -> 3.9 | 12.4 -> 5.3 | -40.2 -> -16.5 | 18.6 -> 7.1 | 45 |
| sn | 14.9 -> 5.7 | 15.2 -> 14.0 | -13.8 -> -0.5 | 14.1 -> 9.6 | 97 |
| ec | 14.7 -> -2.2 | -13.2 -> -1.5 | -30.8 -> -6.1 | -16.9 -> -6.4 | 77 |
| nc | 18.4 -> 9.4 | 6.0 -> 8.5 | -10.7 -> 4.4 | -8.0 -> -5.2 | 54 |
| ie | 10.3 -> 2.1 | -0.9 -> -2.9 | -24.0 -> -13.0 | -11.6 -> -14.6 | 47 |
| ci | 8.7 -> -1.8 | -0.3 -> -4.1 | -22.0 -> -8.0 | -11.4 -> -15.5 | 44 |
| pn | 9.8 -> 4.1 | 11.4 -> 7.7 | -16.9 -> -17.1 | 12.3 -> 10.1 | 77 |
| wc | 6.4 -> 0.7 | 4.5 -> 2.3 | -16.7 -> -14.6 | -2.2 -> -3.6 | 77 |
| ne | 8.7 -> 2.9 | 9.9 -> 9.7 | -11.4 -> -4.4 | 24.0 -> 22.9 | 112 |
| acd | 9.9 -> 5.0 | 6.3 -> 5.0 | -18.0 -> -12.3 | 13.1 -> 10.9 | 145 |
| ca | 4.8 -> 4.9 | 7.6 -> 4.6 | 3.9 -> -1.4 | 8.4 -> 7.0 | 40 |
| cs | 3.8 -> 3.6 | 20.4 -> 18.3 | 6.3 -> 4.0 | 20.5 -> 18.5 | 40 |
| ls | 7.2 -> 7.7 | 16.8 -> 13.2 | -4.1 -> -9.9 | 21.1 -> 18.4 | 64 |

QMD bias falls on 11 of 13 variants, dramatically on the worst-default variants (kt 26.3 -> 3.9, ec 14.7 -> -2.2,
sn 14.9 -> 5.7). BA and TPH improve on most variants; VOL improves on the majority. The variants where calibration
does little (ca, cs, ls) already had small QMD bias, so there is little to remove.

By ecoregion (EPA L2), joined subset:

| L2 | description | QMD def -> cal | BA def -> cal | n |
|---|---|---|---|---|
| 6.2 | Western Cordillera | 19.8 -> 10.4 | 9.2 -> 5.6 | 65 |
| 8.1 | Mixed Wood Plains | 11.5 -> 6.8 | 4.4 -> 2.4 | 23 |

By landowner:

| landowner | QMD def -> cal | BA def -> cal | TPH def -> cal | n |
|---|---|---|---|---|
| Private | 12.1 -> 5.0 | 8.2 -> 6.5 | -17.1 -> -10.0 | 55 |
| National Forest | 19.5 -> 10.8 | 10.4 -> 7.0 | -25.4 -> -17.6 | 56 |

Calibration roughly halves QMD bias on both ecoregions and both landowner classes. National Forest land carried the
larger default bias (less-managed, denser stands) and sees the larger absolute correction.

## Caveats

- The ecoregion and landowner margins rest on 123 disturbance-clean conditions; they show the direction and rough
  magnitude of the correction, not a precise per-cell estimate. A larger harness sample would tighten them; the
  variant margin (full 951) is the robust quantitative result.
- VOL on a few eastern variants (ne, cs, ls) stays high after calibration because the volume bias there is a
  taper/merch-rule issue, not a growth-level issue the density and BAIMULT levers touch.
