# Post-refit calibrated NE vs ACD comparison v3

Source: calibrated_ne_vs_acd_v3/fia_benchmark_results_refit_only.csv
Origin: SLURM job 9914785 pass 1 (round 16, before round-17 ACD fallback patch)
Note: 10022214 in flight with round-17 patches; v3 will be updated when complete

## Per-variant BA metrics

| VARIANT | n | calib RMSE | default RMSE | calib R^2 | default R^2 | bias_pct calib | bias_pct default |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AK | 56 | 52.03 | 68.37 | 0.825 | 0.697 | +267.30 | +221.81 |
| BM | 3369 | 24.49 | 27.29 | 0.793 | 0.743 | +26.06 | +13.17 |
| CA | 603 | 44.18 | 48.43 | 0.791 | 0.749 | +140.58 | +123.43 |
| EC | 2792 | 34.46 | 37.66 | 0.801 | 0.763 | +54.75 | +50.69 |
| IE | 650 | 31.73 | 35.66 | 0.809 | 0.759 | +69.79 | +57.31 |
| LS | 26746 | 19.28 | 19.09 | 0.823 | 0.826 | +55.80 | +47.67 |
| NC | 418 | 52.38 | 62.83 | 0.795 | 0.706 | +107.71 | +83.69 |
| NE | 14717 | 21.15 | 21.65 | 0.842 | 0.834 | +20.06 | +11.33 |
| PN | 2822 | 61.02 | 63.59 | 0.647 | 0.617 | +311.58 | +233.47 |
| SN | 36945 | 32.37 | 33.25 | 0.588 | 0.566 | +140.67 | +159.37 |
| SO | 2783 | 27.09 | 30.03 | 0.807 | 0.762 | +28.66 | +15.47 |
| WC | 4447 | 46.93 | 60.13 | 0.827 | 0.717 | +83.44 | +57.33 |
| OVERALL | 96348 | 29.89 | 31.91 | 0.788 | 0.758 | +90.77 | +89.59 |

**OVERALL improvement:** calibrated RMSE 29.89 vs default 31.91 = **6.3% reduction**

## Caveat: ACD missing from this run

This CSV is from job 9914785 which produced the 12-variant table
but without an ACD row (project_condition_default returned NA for
variant_code=ACD). The round-17 patch routes ACD default-path to NE.
Job 10022214 currently running with that patch will produce a 13-
variant table when it completes. v3 will be updated at that point.