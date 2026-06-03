# Scaling FVS projections into PERSEUS: scenarios, statewide estimates, TreeMap

Scoping memo for: (1) the three-way CONUS FVS comparison, (2) disturbance +
harvest scenario runs feeding statewide estimates to the PERSEUS dashboard, and
(3) a TreeMap2022 spatially-explicit track with a FIADB-vs-TreeMap comparison.
Grounded in what is already on Cardinal (June 2026).

## 0. What is already in place (so we build, not rebuild)

* **FIADB campaign running now.** SLURM array 11228741, 381 tasks = the full
  FIADB (~1.8M stands, 5000/batch), each stand projected 100 yr in **default
  and calibrated** FVS (AGB via NSBE). Idempotent reruns. This is the per-plot
  trajectory engine everything else reuses.
* **gompit "national" FVS** validated in-engine across 23 buildable variants
  (env or `GOMPMORT` keyword). This is the third arm.
* **PERSEUS dashboard data model** (`perseus_wire/api`, schema `perseus_api_v1`):
  per-**state** series JSON, 70 metrics (agb_dry, agc_live_total, vol_stem,
  harvest_c_yr, nbp_yr, mean_stand_age, old_forest_share_120, disturbance rates,
  diversity...), each metric -> scenario class (currently only
  `managed (harvest)`) -> engine (currently only **CEM**, with RCP45/85 + econ
  variants), time series `[year, mean, lo, hi]` 2004-2099. **FVS is not yet an
  engine.** So our work plugs in as new engines.
* **Disturbance rasters**: `p_disturbance_{2016,2020,2022}.tif`, 30 m CONUS
  (154179 x 97279), per-pixel **probability** 0.014-0.719, mean 0.178.
* **TreeMap2022**: `reference_rasters/TREEMAP/TM2022/TreeMap2022_CONUS.tif` (30 m,
  official USFS), `TreeMap2022_CONUS_Tree_Table.csv` (donor-plot tree lists),
  data dictionary; plus `TREEMAP_outputs_v5/attrs/2022/tm_id_plain.tif` = the
  **donor plot ID per pixel**, and `qmd/topht/fortypcd/ecoregion` attribute
  rasters.
* **Maine prototype**: `ME_AGB_Map_Comparison/` (fia, hexagons, lidar, processed
  maps) -- a seed for the multi-scale FIADB-vs-map comparison.

## 1. Three FVS engines into the dashboard (near-term, low cost)

Default (status quo), calibrated (Bayesian), gompit (national mortality) become
three FVS engines. Per arm: take each stand's AGB trajectory, convert to the
dashboard metrics (AGC = AGB x ~0.47; add stem volume and mean age from the FVS
summary tables), expand by FIA plot weight (TPA/EXPNS) to per-acre -> total,
aggregate to **state x year x metric**, emit as `series/<STATE>.json` engine
class `FVS`, scenario `potential (no management)` with a MC/uncertainty band.

* Default + calibrated: data is being produced now; aggregation is the only new
  piece (reuses NSBE + FIA expansion already in the harness).
* gompit arm: rerun the same campaign with `FVS_LIB_DIR=fvs_gompit/lib` +
  `FVS_GOMPIT=1`, config = default growth (gompit owns mortality). ~Same compute
  as one campaign arm.

Cost: LOW. Deliverable: FVS as a first-class engine alongside CEM, three
parameterizations, no-management baseline, all states.

## 2. Disturbance + harvest scenarios -> statewide estimates

The current runs are the **potential** (no disturbance, no harvest) trajectory.
To populate the dashboard's scenario axis:

**Disturbance (data in hand).** For each stand, sample its pixel value from
`p_disturbance_2022.tif` (annual/periodic disturbance probability). Per 5-yr
cycle draw a Bernoulli event; on an event apply an FVS disturbance pulse
(mortality fraction by disturbance agent via `FIXMORT`/fire-fuels extension, or
a partial stand-replacing event). Monte Carlo N draws -> the `[mean, lo, hi]`
band the dashboard expects. Maps directly to metrics `any_disturbance_rate_pct`,
`insect/disease/animal/human_rate_state_pct`.

**Harvest.** No standalone harvest-probability raster exists yet. Three options,
in increasing realism:
  (a) **BAU rotation** by state/forest-type (simple, matches CEM's
      `managed (harvest)` class for apples-to-apples) -- recommended first cut;
  (b) **FIA-derived** harvest rates from the human-disturbance component
      (`human_rate_state_pct`) and FIA TRTCD/DSTRBCD removal records;
  (c) **econ/stumpage-driven** using `perseus_wire/api/stumpage.json`.
Output `harvest_c_yr` flux + the post-harvest stock trajectory.

Cost: MEDIUM. Reuses the campaign harness; new work is the FVS scenario-keyword
injection, the raster sampling per stand, and the MC wrapper.

## 3. TreeMap2022 spatially-explicit -- the key feasibility result

Naive per-pixel FVS over ~5 billion forested 30 m pixels is infeasible. **But
TreeMap is an imputation: every pixel carries a donor FIA plot ID
(`tm_id_plain.tif`), and the donor set is essentially the FIA plots the campaign
is already projecting.** So:

> spatially-explicit TreeMap projection = run FVS **once per unique donor plot**
> (already happening in the FIADB campaign) + a raster **join** painting each
> 30 m pixel with its donor plot's trajectory.

Compute for the FVS side is therefore ~**zero marginal** over the running
campaign -- the cost is the tm_id -> plot-CN join and writing per-cycle AGB/AGC
rasters (~1-2 GB each x 20 cycles, the cspi `aligned_30m`/`conus_tiles` infra
already exists). Full independent per-pixel FVS is the infeasible path and is
**not** needed.

Cost: MEDIUM (raster join + storage), not the "billions of runs" people fear.
Recommend a 1-2 state pilot first (ME, where the comparison prototype exists)
before CONUS-wide rasters.

## 4. FIADB vs TreeMap comparison across spatial scales

* **FIADB**: design-based plot-expansion estimates -- statistically rigorous at
  state/region; no sub-county spatial detail.
* **TreeMap**: wall-to-wall 30 m, aggregable to any polygon; carries imputation
  error.

Comparison design: aggregate both to **state** (should agree if TreeMap is
unbiased -> validates the painting), then to **county / HUC / AOI** (TreeMap adds
detail FIADB cannot resolve; quantify where imputation bias emerges by comparing
TreeMap areal sums to FIA design estimates at progressively finer scales). The
Maine prototype (`ME_AGB_Map_Comparison`, with lidar) is the calibration anchor.
Deliverable: a multi-scale agreement curve -- where TreeMap adds value vs FIADB,
and the scale below which imputation error dominates.

## 5. Recommended sequence (cost-ordered, autopilot-able where marked)

1. [running] finish FIADB default + calibrated. **AUTO**
2. add gompit arm (rerun w/ gompit lib + env). **AUTO**, ~1 campaign-arm cost
3. aggregate 3 arms -> `perseus_api_v1` state series JSON (FVS engine). **AUTO**
4. disturbance scenario (raster-driven MC) on FIADB. MED; design then auto
5. harvest scenario -- start BAU to match CEM. MED; needs decision (3b)
6. TreeMap pilot (ME): tm_id -> CN join, paint, per-cycle rasters. MED
7. FIADB-vs-TreeMap multi-scale comparison (ME -> CONUS). MED

Steps 1-3 I can run on autopilot now. 4-7 each have one design decision below.

## 6. Decisions for Aaron (each unblocks a step)

* **Harvest definition** (step 5): BAU rotation by forest type (fast,
  CEM-comparable) vs FIA-derived removals vs econ/stumpage? Recommend BAU first.
* **Disturbance application** (step 4): mortality-pulse magnitude per agent --
  use FIA DSTRBCD-implied severities, or a single generic severity to start?
* **Climate** : CEM has RCP45/85. Do the FVS engines stay climate-static
  (cleaner attribution of the mortality/calibration difference), or do we add
  the FVS Climate-FVS extension to match CEM's scenarios?
* **TreeMap scope**: CONUS-wide now, or ME/Lake-States pilot first (recommended)?

## 7. Bottom line

All four threads are feasible and most of the data + harness already exist. The
two non-obvious wins: (a) the three FVS arms plug into PERSEUS as engines with
low marginal effort, and (b) **TreeMap spatially-explicit reuses the running
FIADB per-plot compute almost for free via the donor-plot join** -- so the
FIADB-vs-TreeMap comparison is within reach without a separate billion-run
campaign. The realism step that needs your input is the harvest/disturbance
scenario definition.
