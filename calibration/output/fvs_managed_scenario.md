# FVS managed (harvest) scenario for PERSEUS

The reserve (no-harvest) FVS engines are the upper bound of forest carbon. The
managed (harvest) scenario is the realistic counterpart: it applies the
data-driven harvest regime from the `conus_hcs` rasters plus background
disturbance to each state's reserve trajectory.

## Rates (data-driven, per state)

`build_state_harvest_rates.R` samples the canonical harvest maps at the actual
FIA plot locations (so each state's rate reflects where its own plots sit), from
`conus_hcs/data/analytic/maps_conus_canonical`:

* `conus_expected_ba_removed_annual` -> `harvest_frac_yr` (expected fraction of
  standing basal area removed per year; occurrence x intensity in one layer).
* `conus_p_partial/clearcut/stand_replacement_annual` for context.
* `p_disturbance_2022.tif` -> a multi-year disturbance susceptibility probability
  (NOT annual), annualized as `a = 1-(1-p)^(1/W)`, W = 20 yr.

Sampled rates are realistic: harvest removes ~0.8 to 2.1 %/yr of basal area
(ME 1.7%, OR 1.8%, GA 2.0%, MN 0.85%); annualized disturbance 0.5 to 1.6 %/yr.

## Model (transparent, repeatable)

`fvs_managed_scenario.py` walks each reserve density path in 5-yr steps:

```
growth5  = d_reserve[t] - d_reserve[t-5]      # the engine's own increment
removed5 = 5 * (h + a) * d_managed[t-5]
d_managed[t] = max(0, d_managed[t-5] + growth5 - removed5)
harvest_flux[t] = h * d_managed[t]            # harvested carbon density
```

It uses the FVS engine's own growth increment and removes the harvested +
disturbed fraction each step, so the managed path tracks below reserve at a
working-forest steady state. Applied to `agc_live_total` and `agb_dry`; the
harvest carbon flux becomes `harvest_c_yr`.

## Result (2125, Mg C/ha; reserve -> managed)

| state | default | calibrated | gompit |
|-------|--------:|-----------:|-------:|
| ME | 269 -> 68 | 234 -> 57 | 179 -> 22 |
| GA | 170 -> 31 | 129 -> 22 | 110 -> 13 |
| MN | 201 -> 86 | 147 -> 61 | 106 -> 32 |
| OR | 235 -> 46 | 173 -> 31 | 171 -> 26 |

Managed standing carbon settles at ~15 to 45% of the no-harvest upper bound,
state-specific by harvest intensity, with gompit lowest (density-dependent
mortality plus harvest). Maine harvest flux ~4.7 Tg C/yr.

## On the dashboard

`fvs_perseus_merge.py` with `FVS_MANAGED_ROOT` set injects a `managed (harvest)`
bucket + `harvest_c_yr` for all three FVS national models, alongside the reserve
bucket. Reserve curves are unchanged.

## Caveats / refinements

* The decrement model is first-order (no explicit stand-age / rotation
  structure); the YC engine carries owner-rotation harvest for comparison.
* The disturbance annualization window (20 yr) is an assumption; disturbance is a
  secondary term (harvest dominates).
* Rates are climate-static (2024 harvest economics), per the agreed scope.
