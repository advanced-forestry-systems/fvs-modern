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

## AK (BC) engine validation: pilot findings (status)

The FVS-AK machinery and species crosswalk are validated; the calibration number is not yet trustworthy
because of a treeinit ingestion-scaling detail specific to the AK/database path. Precisely:

- Species crosswalk: 100% of coastal-BC trees map to FIA species codes the AK variant accepts (TSUG.HET,
  THUJ.PLI, PICE.SIT, PSEU.MEN, ABIE.AMA, PINU.CON, etc.). No coverage loss.
- Clean pairs: 2,451 protocol-consistent BC remeasurement pairs (magplot_ak_bc_pairs_clean.csv), filtered to
  sites with the same subplot set and DBH tag limit across visits and plausible BA change (0.7 to 1.8x). Mean
  interval 24.7 y, BA 51.4 -> 60.6 m2/ha, median annual BA increment 0.34 m2/ha/yr, 16% mortality-driven
  decline. This is clean, realistic growth data, ready to calibrate against.
- The blocker: a single-stand trace shows FVS-AK reads all records (e.g. 86 of 86) but reports Tpa 60 from a
  430-TPA treeinit (about 7x low) and runs 1-year cycles despite TIMEINT 0 10. The same fvs2py/standalone
  treeinit path is level-unbiased for ACD/NE (NB validation: BA 17.3 vs 17.3 obs), so this is an AK-variant
  database-ingestion detail (tree_count expansion and cycle-length handling under region 10), not a data
  problem. It is one focused engine-debug fix away from a real AK bias number.

ACD remains the validated maritime reference: the prior sprint's NB MAGPlot validation is near-unbiased on BA
(magplot_insource_v17_results.csv: -0.04% BA bias, R2 0.88; QMD +9%, TPH -7%). The residual QMD over-prediction
is the ACD calibration target, consistent with the CONUS recalibration pattern.
