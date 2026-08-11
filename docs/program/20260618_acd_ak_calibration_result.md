# ACD and AK calibration against Canadian NFI (MAGPlot): result
2026-06-18. Closes the ACD calibration; precisely localizes the one remaining AK engine bug.

## ACD (Acadian, New Brunswick): calibrated and validated

Default FVS-Acadian projected against 262 protocol-consistent NB remeasurement stands (5-year intervals),
standalone-binary path (validated machinery), observed compiled from MAGPlot stem-expansion:

| metric | default bias | proj vs obs | calibration factor (obs/proj) |
|---|---|---|---|
| BA  | -0.0% | 17.3 vs 17.3 m2/ha | 1.000 (no adjustment needed) |
| QMD | +8.9% | 12.9 vs 11.9 cm | 0.918 |
| TPH | -7.3% | 1675 vs 1807 /ha | 1.079 |
| BA R2 | 0.878 | | |

Reading: default Acadian basal-area growth on Canadian maritime forest is essentially unbiased and well
correlated (R2 0.88), confirming the engine transfers across the international border with no BA-level
adjustment. The residual is a modest diameter-vs-density split: QMD over-predicted ~9% and TPH under-predicted
~7% (the model thins slightly too fast and grows survivors slightly too large), the same signature the CONUS
recalibration corrects. The QMD 0.918 / TPH 1.079 factors are the ACD MAGPlot calibration; they match the
direction of the keyword-calibration QMD correction already applied CONUS-wide. ACD calibration is done.

## AK (Alaska / coastal BC): data and crosswalk ready; one engine bug localized

Everything up to the engine projection is complete and verified:

- 2,451 protocol-consistent BC remeasurement pairs (matched subplots and DBH tag limit across visits,
  plausible BA change), mean interval 24.7 y, BA 51.4 -> 60.6 m2/ha, median annual BA increment 0.34 m2/ha/yr.
- Species crosswalk 100% (BC genus.epithet -> FIA codes the AK variant accepts: western hemlock, redcedar,
  Sitka spruce, Douglas-fir, Pacific silver fir, lodgepole, etc.).
- FVS-AK binary works: it reads and projects the validated NB Acadian tree lists correctly (BA 34.9 m2/ha,
  TPA 1548 for both ACD and AK on the same NB stand).

The bug, localized precisely by controlled tests:

- FVS reads DBH correctly for BC lists (returned QMD always equals the input QMD to the cm).
- FVS applies TREE_COUNT (the stand expansion) only up to a stand-specific multiplicative factor that is 1.0
  for NB lists but about 0.14 for BC lists. Scaling tree_count by 10 scales the reported TPA by exactly 10
  (60 -> 600), so the engine does use the column, but renormalizes it by a wrong constant for these lists.
  The constant shifts with the record structure (collapsing 86 records to 2 changed it from 0.14 to 0.08),
  which points at the AK variant's per-record sampling/expansion renormalization (tree tripling or the
  stockable-area normalization), not at the input data, the SQLite schema, the column case, the dtypes, or
  the species mapping, all of which were ruled out with side-by-side NB-vs-BC tests.
- Stand setup is identical and correct in both cases (NUMBER OF PLOTS = 1, sampling weight 1.0,
  stockable 1.0), so the divisor is internal to the expansion code path.

This is a maintainer-level FVS-source fix in the AK variant's inventory expansion (same class as the Route A
tree-init blocker), not a data or harness problem. The BC observed-growth dataset is ready to validate against
the moment that constant is corrected; no recompute of the pairs is needed. Diagnostic scripts:
diagnostics_2026-06-16/magplot/ (ak_runner.py, ak_scale.py, bc_clean.py).

## Net

ACD: calibrated and validated against Canadian NFI (the maritime cross-border test passes). AK: fully staged
and the engine bug pinpointed to one expansion-renormalization constant in the AK variant; one focused source
fix from a number.
