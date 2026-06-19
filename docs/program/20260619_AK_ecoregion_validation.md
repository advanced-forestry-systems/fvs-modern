# FVS-Alaska validation against Canadian NFI, stratified by NA Level III ecoregion
2026-06-19. OODA autopilot run. Objective: validate (and quantify the calibration need for) the FVS Alaska
variant against MAGPlot BC remeasurement growth, stratified by North American CEC ecoregion rather than by a
crude coastal longitude cut, so the comparison uses BC areas ecologically analogous to SE Alaska.

`[MODE: Autonomous | STACK: R/data.table + pure-Python spatial join + Python FVS driver | HARDWARE: OSC Cardinal SLURM]`

## OBSERVE

- MAGPlot BC: 16,517 sites, 2,451 protocol-consistent remeasurement pairs (matched subplots and DBH tag limit
  across visits, plausible BA change), mean interval ~24 yr.
- Ecoregion layer: NA CEC Level III (NA_Eco_L3_WGS84), 2,548 polygons, joined to BC site lat/long.
- FVS-AK engine: standalone binary, in-database tree list, one stand per subprocess.

## ORIENT

Two compute vulnerabilities resolved this run:
1. The PROJ/GEOS/GDAL stack for R sf and Python geopandas is broken on the cluster (missing libproj.so.25).
   Bypassed with a pure-Python point-in-polygon join (pyshp reads the WGS84 polygons, matplotlib.path tests
   membership), no geospatial system libraries required.
2. FVS-AK under-expands certain BC tree lists. Diagnosed earlier: the engine scales the whole stand linearly by
   INV_PLOT_SIZE, and at the default value it under-counts sparse lists. Corrected per stand with a two-pass
   auto-scale (run once, measure the realized TPA, set INV_PLOT_SIZE so the initial stand matches the compiled
   TPA), then kept only stands whose initial BA reproduces the observed initial BA within tolerance
   (clean-ingestion filter). This removes mis-ingested stands rather than letting them bias the result.

## DECIDE

Run FVS-AK default on the clean-ingestion BC stands, compute bias against observed t2 by NA Level I ecoregion,
and report increment bias (projected dBA minus observed dBA, which cancels any residual tag-limit offset) and
standing-level bias. Marine West Coast Forest is the SE-Alaska-analog ecoregion (SE Alaska is itself Marine
West Coast Forest); the interior ecoregions are contrast groups.

## RESULT

![AK bias by ecoregion](fig_ak_ecoregion_bias.png)

| NA Level I ecoregion | n | interval | standing BA bias | BA increment bias | proj dBA | obs dBA (m2/ha) |
|---|---|---|---|---|---|---|
| Marine West Coast Forest (AK analog) | 11 | 22 yr | -35.6% | -74.1% | 3.43 | 13.25 |
| Northwestern Forested Mountains | 104 | 22 yr | -29.2% | -77.3% | 2.08 | 9.19 |
| Northern Forests | 17 | 20 yr | -36.3% | -90.1% | 0.96 | 9.75 |
| Taiga | 16 | 13 yr | -29.2% | -108.5% | -0.51 | 5.96 |
| North American Deserts | 21 | 35 yr | -12.9% | -9.3% | 8.36 | 9.22 |

Findings:

- The default FVS Alaska variant substantially under-predicts growth in its own analog ecoregion. On Marine
  West Coast Forest the projected basal-area increment is 3.43 m2/ha over 22 yr against an observed 13.25, an
  increment bias of -74%. So the under-prediction reported earlier on the un-stratified coastal set is not an
  artifact of fast-growing southern coastal sites; it holds in the ecoregion that actually matches SE Alaska.
- The under-prediction is broad across productive BC ecoregions: Northwestern Forested Mountains -77%
  (n=104, the best-sampled group), Northern Forests -90%, Taiga -109% (the variant projects net basal-area
  loss where the stands gained 6 m2/ha).
- The one ecoregion where FVS-AK tracks observed growth is North American Deserts (dry interior BC, e.g. the
  Thompson-Okanagan), increment bias only -9%. These low-productivity stands resemble the slow SE-Alaska
  forests the variant was parameterized on, which is why the defaults fit them.
- Implication for calibration: the correction is ecoregion (productivity) dependent, not a single multiplier.
  The implied basal-area-increment multiplier is about 3.9x for Marine West Coast Forest (13.25/3.43) and
  about 1.0x for dry interior. A constant BAIMULT would over-correct the dry interior.

Within Marine West Coast Forest the Level III split is the exact SE-Alaska analog (Coastal Western
Hemlock-Sitka Spruce Forests, n=3; Pacific and Nass Ranges, n=7), but those cells are small; the Level I
estimate (n=11) is the reportable analog number.

## Caveats

- The clean-ingestion filter passes only stands FVS-AK ingests correctly, so n per ecoregion is modest
  (Marine West Coast Forest n=11). The interior groups are better sampled. The direction and rough magnitude
  are robust; precise per-ecoregion multipliers need more clean stands or the FVS-AK expansion fix so the
  filter can be relaxed.
- Applying the correction in-engine is still open: the BAIMULT keyword did not change growth as placed
  (multipliers 1, 3, 5 gave identical output), so the calibration must go through the project's keyword
  pipeline or the config_loader calibrated path, not a raw BAIMULT line. That is the next step.

## Scripts and data

- eco_join_py.py (pure-Python ecoregion join), bc_site_ecoregion.csv (site to ecoregion lookup)
- ak_eco_par.py + sub_eco.sh (parallel SLURM validation), ak_eco_validate_results.csv
- fig_ak_ecoregion_bias.png
