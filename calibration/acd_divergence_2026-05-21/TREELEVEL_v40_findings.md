# Tree-level reconciliation v40: TPA-preserving, size-weighted

2026-05-29. v34's uniform EXPF scaling crashes TPA by ~10 percent on v36/v38
samples because BA and TPA scale together when EXPF is multiplied uniformly.
v40 fixes this by scaling EXPF as a function of DBH so BA reaches the target
while sum(EXPF) stays at the raw value. Biological story: the model is
keeping too many big trees alive; reduce big-tree EXPF and shift weight
toward smaller trees.

## The math

Per stand, parameterize the per-tree scale factor as

  s_i = 1 + beta * (DBH_mean - DBH_i) / DBH_mean

where `DBH_mean` is the EXPF-weighted mean DBH. By construction:

  sum(s_i * EXPF_i) = sum(EXPF_i)     # TPA preserved exactly

Solve for beta from the BA target:

  beta = (BA_target - BA_raw) * DBH_mean / (DBH_mean * M2 - M3)

where M2 = sum(EXPF * DBH^2), M3 = sum(EXPF * DBH^3). For typical stands
M3 > DBH_mean * M2 (because moments grow super-linearly with DBH), so the
denominator is negative and beta > 0 when reducing BA. Small trees end up
with s_i > 1 (more weight), big trees with s_i < 1 (less weight).

Safety bounds: s_min = 0.1, s_max = 2.0 prevent pathological cases.

## v36 comparison (n=300, real data)

| metric    | raw      | v34 uniform | v40 TPA-preserving |
|-----------|----------|-------------|---------------------|
| BA bias   | +10.08   | -0.30       | **-1.78**          |
| TPA bias  | -0.66    | -11.44      | **-0.66 (preserved)** |
| BA R^2    | 0.495    | 0.549       | 0.535              |
| TPA R^2   | 0.745    | 0.699       | **0.745**          |

v40 preserves TPA exactly while closing BA to within 2 percent. The cost on
BA R^2 is small (0.535 vs v34's 0.549).

## Per-quartile (BA_t1)

| q | BA_t1 mean | raw BA | v34 BA | v40 BA | raw TPA | v34 TPA | v40 TPA |
|---|------------|--------|--------|--------|---------|---------|---------|
| Q1 | 27 | +60.7 | +21.5 | **+6.3** | +13.4 | -16.6 | +13.4 |
| Q2 | 70 | +9.4  | -10.7 | -10.8 | -7.0  | -23.9 | -7.0  |
| Q3 | 110 | +7.1 | -0.9  | -0.9  | -2.5  | -9.5  | -2.5  |
| Q4 | 167 | +0.8 | +0.2  | +0.2  | -1.1  | -1.5  | -1.1  |

v40 is substantially better in Q1 (+6.3 vs +21.5) because v34's scale_floor =
0.7 caps how aggressively low-density stands can be thinned uniformly. v40
has no scale floor (s_i bounded only at [0.1, 2.0] for safety), so it can
close more of the Q1 overshoot.

## QMD behavior

v34: QMD invariant (uniform EXPF scaling preserves QMD).
v40: QMD decreases by ~5 percent because big trees are thinned preferentially.

This is the main trade-off. Volume tables and biomass equations that are
QMD-sensitive may prefer v34 for that reason. For TPA-driven downstream
(counts, carbon per acre, harvest density decisions), v40 is the right
choice.

## Recommendation

v40 is the new production default. The API exposes both:

    # Default v40 (TPA-preserving)
    res <- apply_density_correction_treelist_tpa(tree_df, BA_t1_by_stand)

    # v34 (uniform, QMD-preserving)
    res <- apply_density_correction_treelist(tree_df, BA_t1_by_stand,
                                              scale_floor = 0.7)

Both functions return `tree_corrected` and `stand_summary` with the same
structure.

## Caveats unchanged

- Apply only at the production posture (12.3.9 + CSI_SCALE = 0.7).
- Do NOT apply to Canadian MAGPlot (v17 baseline already close to zero;
  the correction would mostly degrade fit, though the v39 50-stand subset
  showed sample-dependent behavior).
- Not a substitute for refitting the BAL coefficient in the Kuehne dDBH
  equation. Paper-sized structural fix remains the long-term goal.

## Files

  apply_density_correction.R   updated to expose both v34 and v40
  TREELEVEL_v40_findings.md    this memo
  v36_perplot.csv              the data used for the comparison
