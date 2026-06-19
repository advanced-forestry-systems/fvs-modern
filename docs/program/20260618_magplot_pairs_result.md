# MAGPlot remeasurement pairs for ACD and AK calibration
2026-06-18. The Canadian NFI MAGPlot package (open.canada.ca dataset 1a73441d) is staged and fully extracted
on Cardinal at ~/magplot/. This step builds stand-level remeasurement pairs from source for the two target
variants and confirms both are calibration-ready.

## Real MAGPlot schema (confirmed from source)

- magp_sites: magp_site_id, province (BC 16,517; QC 12,892; NB 7,609; NL 1,350; YT 327; NT 262; SK 171; MB 2),
  lat/long, elevation. NS and PE are NOT in the package, so the Acadian target is NB only.
- magp_trees (1.5 GB): magp_site_id, plot_type_id, tree_num, meas_num, tree_id, species_g (genus, e.g. ABIE),
  species_s (epithet, e.g. BAL), dbh (cm), height (m), stem_ha (expansion factor), tree_status.
- magp_tree_header: per (site, plot, meas_num) meas_year. meas_year is consistent within a site-measurement,
  which is the reliable join key for the interval (the plot_id string format differs between trees and header).
- tree_status differs by region: NB uses L / DS; BC uses LS / LF / DS. Live = status begins with "L".

## Variant mapping

- ACD (Acadian): NB (New Brunswick). 7,609 sites. Species are the Acadian set (ABIE BAL balsam fir, PICE RUB
  red spruce, PICE MAR black spruce, ACER RUB red maple, ACER SAC sugar maple, BETU ALL yellow birch,
  BETU PAP paper birch, THUJ OCC northern white cedar, etc.).
- AK (Alaska / coastal): BC (British Columbia). 16,517 sites. Top species PINU CON (lodgepole pine),
  PSEU MEN (Douglas-fir), TSUG HET (western hemlock), THUJ PLI (western redcedar), PICE GLA, ABIE LAS,
  PICE SIT (Sitka spruce), ABIE AMA (Pacific silver fir). The AK variant is the coastal-BC analog; interior
  BC species (lodgepole, Douglas-fir, larch, ponderosa) are out of the AK variant's native range and would map
  to ie/East-Cascades-type variants, so AK calibration restricts to the coastal subset.

## Pairs built

- ACD / NB: 263 stand pairs (pre-existing from the May sprint, verify/magplot_pairs.csv), 5-year intervals.
  Prior FVS Acadian validation (magplot_insource_v17_results.csv, n=262): BA bias -0.04% (essentially
  unbiased), BA R2 0.88, QMD +9% (12.95 vs 11.89 cm obs), TPH -7% (1675 vs 1807). The default Acadian engine
  already tracks NB MAGPlot well; calibration targets the residual QMD over-prediction.
- AK / BC: 8,236 stand pairs built this session (magplot_bc_pairs.csv). Mean interval 21.3 years (NFI long
  remeasurement), BA 40.5 -> 47.8 m2/ha, QMD_t2 23.1 cm, TPH_t1 2,201. Observed growth is calibration-ready;
  the AK engine projection/validation is the next compute step (adapt magplot_fvs_runner.py to the AK variant
  and the coastal-BC species crosswalk).

## Scripts

- build_bc_pairs.py: streams the 1.5 GB tree table, filters to BC sites, writes magp_trees_bc.csv (3.83M rows).
- bc_pairs_v2.py: builds the stand pairs (site x plot_type), joins meas_year on (site, meas_num), live = status
  starts with L and stem_ha > 0, BA/TPH/QMD from the stem_ha expansion.
