# Acadian annual calibration table: how it was produced

`acd_annual_calibration.csv` holds annual diameter increment multipliers
(`dDBH.mult`) for the Acadian model, fit against real ACD FIA remeasurement.
This is the table the bridge patch reads in place of the FVS-NE `baimult`.

## Method

`fit_acd_annual_calibration.R` reads the ACD remeasurement records
(`fvs-conus/data/processed/acd/diameter_growth.csv`, 565,776 usable trees,
5 year remeasurement). For each tree it iterates the real Acadian `dDBH_fun`
annually over the tree's interval and finds the per species annual multiplier
that minimizes squared error between predicted and observed periodic diameter
growth. Species random effects come from `ddbh.fun.spp` in `AcadianGY.R`.

## Result

34 Acadian species fitted (the remaining valid codes are filled with 1.0 so the
bridge always finds a row). Mean diameter growth RMSE drops from 1.108 cm to
1.015 cm, an 8.4 percent reduction. Most multipliers are above 1, because the
uncalibrated Acadian equations under predict annual diameter growth on these
plots (for example balsam fir observed 1.30 vs predicted 0.97 cm per year, red
oak 1.98 vs 1.00). This is the quantitative form of the mismatch you observed.

## Assumptions and caveats (read before shipping)

1. Competition split: the records carry total `BAL` only, so the larger softwood
   term `BAL.SW` was set to total BAL and `BAL.HW` to zero. The `BAL.HW`
   coefficient is small (b5 = -0.0177), but for hardwoods this slightly inflates
   modeled competition and pushes their multipliers upward. With a proper
   softwood and hardwood split the hardwood multipliers would come down somewhat.
2. Climate site index was mapped from FVS site index in feet using the same
   piecewise map the bridge uses, defaulting to 12 when missing.
3. Competition is held at the start of period value across the interval, which
   matches how the model applies within period competition; it does not evolve
   BAL year to year inside this fit.
4. These are bridge calibration factors under the stated assumptions, not new
   biological parameters. Treat them as a working table to remove the level and
   species bias, and refine once a softwood and hardwood BAL split is available.
5. Height (`dHt.mult`) and mortality (`mort.mult`) are left at 1. Fit them with
   the same machinery on `HT_t1/HT_t2` and on `STATUSCD` survival when needed.

## How the bridge consumes it

Per `bridge_patch_baimult_scale.md`, `make_fvs_calib()` reads this CSV and sets
`dDBH.mult`, `dHt.mult`, `mort.mult` from it (keyed on the Acadian alpha code
`SP`), keeping only the size caps from FVS. This decouples Acadian calibration
from the FVS-NE multipliers.
