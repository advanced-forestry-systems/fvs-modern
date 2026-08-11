# Dashboard output stress test

`fvs_dashboard_stress.py` scans every state series JSON for the three national
FVS engines and checks the invariants the projections must satisfy:

* finite / non-negative values
* uncertainty band contains its value (`lo <= value <= hi`)
* monotone reserve (no-harvest grows; gompit exempt -- it caps by design)
* scenario ordering: reserve >= extensive >= harvest >= intensive each year
* every managed bucket <= reserve
* engine ordering: gompit <= default at the final year
* harvest flux >= 0
* FIA-anchored states reconcile to fia.json at the anchor year (+-2%)

## Result: clean (49 states x 3 engines x 4 scenarios)

First pass found two things, one real and one false positive:

* **False positive (3):** gompit reserve declines ~1 Tg at a late-succession step
  (FL 858->857, ME 506->505, MN). This is gompit's density-dependent mortality
  capping over-accumulation -- intended, not a bug. The test now exempts gompit
  from the monotone-reserve check.
* **Real fix (60, all Nevada):** a degenerate `[0,0]` band around a positive mean.
  Nevada is almost entirely sparse pinyon-juniper woodland, so >90% of plots are
  near-zero biomass and the 10/90 plot-percentile band collapses below the mean.
  Fix: the merge now clamps the band to contain its value (`lo=min(lo,v)`,
  `hi=max(hi,v)`). Live.

After the fix the scan reports **0 violations**: scenario ordering, engine
ordering, anchor reconciliation, flux signs, and finiteness all hold across the
whole dashboard.

## Side finding: the calibrated engine is informative everywhere

The 0%-posterior-width states (GA/IN/MN) are not redundant with default -- their
calibrated reserve carbon sits well below default (GA -24%, IN -32%, MN -27%, ME
-13%, OR -26%, WA -24%, ID -14% at 2125). Calibration meaningfully moderates the
over-accumulation; the tight posterior just means low *parameter* uncertainty, a
well-constrained fit.

## Reproducibility

`run_fvs_perseus_pipeline.sh` chains the whole flow with the treeinit TPA+height
fix folded in at STAGE 1, so a future campaign starts from correct tree lists,
through aggregate -> managed -> merge -> ribbon -> this stress test (which must
report 0 violations before a push).
