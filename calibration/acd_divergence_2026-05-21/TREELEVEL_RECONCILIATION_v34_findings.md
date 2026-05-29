# Tree-level reconciliation with v33 stand BA correction (v34)

2026-05-29. Extends the v33 stand-level scalar correction to a tree-level
function that scales each tree's EXPF uniformly so the sum of (DBH^2 * EXPF)
matches the corrected BA per stand. Result: downstream consumers (volume
tables, carbon accounting, FIA-format tree lists) get a tree list that is
internally consistent with the BA prediction the bridge ships.

## Two functions

`apply_density_correction(BA_pred, BA_t1, ...)` returns the scalar corrected
stand BA. Same as before.

`apply_density_correction_treelist(tree, BA_t1_by_stand, ...)` takes a tree
data.frame with STAND, DBH, EXPF (and any other columns) plus a per-stand
BA_t1 lookup, and returns:

  tree_corrected   the data.frame with EXPF scaled in place (by stand)
  stand_summary    per-stand BA_raw, BA_corrected, scale_factor, TPA_raw,
                   TPA_corrected diagnostics

Per stand, the scale factor is

    s = max(scale_floor, min(1.0, BA_corrected / BA_raw))

with `scale_floor = 0.7` by default. The upper bound of 1.0 is mathematical;
the v33 correction can only subtract, so s is never above 1.

## Why scale_floor matters

The v33 correction is linear in BA_t1 and clipped to [0, +25] ft^2/ac. For a
realistically calibrated stand the implied scale factor is mild. But for an
edge-case stand where the raw model under-projected BA (BA_pred small)
while BA_t1 was moderate, the implied scale factor can be unrealistically
low (the n=91 v32 sample had min raw scale -1.486 on one pathological
case, q05 = 0.415).

Empirically on v32 (n=91) treating each plot's BA_pred as the raw and
applying the tree-level scaling:

| safeguard      | BA bias  | TPA bias | R^2(BA) |
|----------------|----------|----------|---------|
| raw (no corr)  | +17.92   | +7.29    | (n/a)   |
| no floor       | +5.05    | -5.22    | 0.486   |
| **floor=0.7**  | **+6.98**| **-3.22**| **0.496** |
| floor=0.8      | +8.23    | -1.91    | 0.491   |

floor = 0.7 is the production default: meaningful BA closure (+17.92 ->
+6.98), TPA stays close to observed (sd -3.22), and the per-plot fit R^2
is the best of the three.

## QMD is invariant under uniform EXPF scaling

A useful property: scaling every tree's EXPF in a stand by the same factor
preserves QMD. QMD = sqrt(sum(DBH^2 * EXPF) / sum(EXPF)). Scaling EXPF by s
multiplies numerator and denominator by s, so QMD is unchanged. The
correction therefore moves BA and TPA together, leaves the QMD invariant.

Implications:
- If the model's QMD bias was small (and after the 12.3.7+12.3.8 fix it is),
  this is correct behavior.
- If QMD has a separate bias, the tree-level correction does not address it;
  a separate dDBH-shrink layer would be needed.

## Smoke test (representative, on Cardinal R 4.4.0)

Stand S1, BA_t1 = 50 ft^2/ac (low density, ~4 trees, raw BA = 24.8):
  scale_factor 0.7 (floor active), TPA 180 -> 126 (a 30 percent thin)

Stand S2, BA_t1 = 180 ft^2/ac (high density, ~4 trees, raw BA = 36.9):
  scale_factor 1.0 (no correction), TPA 110 -> 110 unchanged

Both behaviors are as designed.

## What this is and is not

This is a **tree-level reconciliation that preserves the v33 BA constraint
by construction**. The tree list it returns will compute the v33-corrected
BA exactly (subject to scale_floor activation). For downstream uses that
operate on the tree list directly (volume, biomass, carbon, harvest
prescriptions), the resulting numbers are derived from a tree population
that is internally consistent with the corrected stand BA.

It is **NOT a refit of the per-tree mortality model**. We are interpreting
the v33 correction as additional uniform mortality and reducing EXPF
accordingly. The biological story behind that interpretation is plausible
(the bias signal lives at the stand-density level; mortality
under-predicts in dense stands) but it is post-hoc, not a fitted mortality
function.

For the proper structural fix, the BAL coefficient in the Kuehne dDBH
equation or the Glover/Hool mortality functional form would need to be
refit against ME FIA mortality data directly. The bridge correction is
the pragmatic interim.

## API

    source("apply_density_correction.R")

    # Stand-level scalar (unchanged from v33)
    ba_corr <- apply_density_correction(BA_pred, BA_t1)

    # Tree-level reconciliation (new in v34)
    res <- apply_density_correction_treelist(
      tree_df,
      BA_t1_by_stand = c(stand_id = ba_t1_value, ...),
      dbh_units      = "cm",   # or "in"
      scale_floor    = 0.7,    # production default
      upper_cap      = 25,
      lower_cap      = 0
    )
    res$tree_corrected   # the updated tree list
    res$stand_summary    # diagnostics

## Files

  apply_density_correction.R   updated (now exposes both stand- and
                               tree-level entry points)
  TREELEVEL_RECONCILIATION_v34_findings.md  this memo
