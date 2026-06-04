# Plantation-aware management scenarios

Refinement to the FVS managed scenarios: management intensity is not uniform.
**Intensive** management (short-rotation clearcut) applies only to **plantations**
-- plots that FIADB flags as artificially regenerated (`STDORGCD = 1`) -- while
**extensive** (partial harvest) applies to natural stands. This matches how the
forest is actually worked and avoids over-stating intensive management.

## Plantation identification

`build_plantation_flag.py` reads the FIADB COND table: a plot is a plantation if
the majority of its forested area (`CONDPROP_UNADJ`) is in `STDORGCD = 1`
conditions. CONUS: **10.2%** of forest plots are plantations, concentrated where
expected -- GA 28%, MS 27%, OR 20%, vs ME 2.4%, MN low.

## Scenarios (per-plot, per-state)

From the empirically-sampled blended harvest rate `h` (conus_hcs
expected-BA-removed, which already mixes partial + clearcut over the landscape)
and the plantation area fraction `p`, split into extensive `h_ext` and intensive
`h_int = k*h_ext` (k = 1.9, the YC engine's intensive scaling), calibrated so the
plantation-weighted blend reproduces `h`:

    h = (1-p)*h_ext + p*h_int   ->   h_ext = h/((1-p)+p*k),  h_int = k*h_ext

This conserves the observed statewide removal while confining the heavier regime
to plantations. Buckets:

| bucket | regime |
|--------|--------|
| reserve (no harvest) | no harvest (upper bound) |
| managed (extensive) | every managed plot gets `h_ext` (all-partial counterfactual) |
| managed (harvest) | **realistic**: plantation -> `h_int`, natural -> `h_ext` |
| managed (intensive) | every managed plot gets `h_int` (all-intensive bound) |

## Result (2125 live AGC, Mg C/ha, default engine)

| state | plant % | reserve | extensive | **harvest** | intensive |
|-------|--------:|--------:|----------:|------------:|----------:|
| GA | 28% | 1647* | 38 | **33** | 18 |
| OR | 20% | -- | 51 | **47** | 32 |
| ME | 2.4% | -- | 69 | **68** | 35 |
| MN | low | -- | 88 | **87** | 57 |

(*GA reserve shown as state total Tg; density columns are Mg C/ha.) Where
plantations are common (GA, OR) the realistic managed path drops below extensive,
because the plantation fraction carries the heavy regime. Where forests are
natural (ME, MN) realistic and extensive coincide -- intensive only bites where
plantations exist, which is exactly the point.

## On the dashboard

`fvs_managed_v2.py` -> `managed_<ST>.csv` (metric, mgmt, year, value);
`fvs_perseus_merge.py` (with `FVS_MANAGED_ROOT`) injects all four buckets +
`harvest_c_yr` for every engine. Live. Reserve curves unchanged.

## Notes / next

* k=1.9 mirrors the YC intensive scaling; could be sharpened from the conus_hcs
  partial-vs-clearcut intensity layers directly (per-regime rates rather than a
  single ratio).
* Owner x plantation interaction (industrial plantations vs NIPF natural) is a
  natural extension -- the membership table already carries owner.
