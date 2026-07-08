# FVS x PERSEUS: refinements and next steps

State of play: three national FVS engines (default, calibrated, gompit) are live
on the dashboard with reserve and managed (harvest+disturbance) scenarios; the
treeinit expansion bug is fixed; a height-imputed rerun is in flight to make the
2025 inventory year whole; CONUS TreeMap-spatial and the FIADB-vs-TreeMap
multi-scale comparison are done; landowner / ecoregion / state trends with
bootstrap uncertainty are built. From here, in rough priority order.

## Tier 1 -- finish the in-flight chain

1. **Re-anchor at 2025 after the height rerun.** With `treeinit_h` heights, the
   2025 point is whole, so drop the 2030 anchor and anchor the dashboard at the
   true inventory year. Re-aggregate `out_fvs_v3`/`out_gompit_v3`, re-merge,
   re-run strata-trends and the TreeMap comparison on v3.
2. **Publish the landowner / ecoregion / state trend layer** to the dashboard as
   a selectable breakdown (the data exists; needs a UI cut), with the bootstrap
   CI ribbons.

## Tier 2 -- deepen the uncertainty

3. **Bayesian posterior-draw uncertainty.** The calibration carries 500 posterior
   draws per component (config/uncertainty.py UncertaintyEngine). Propagating a
   sample of draws through the projection gives a parametric carbon CI that is
   more defensible than the percentile band, and separates parameter uncertainty
   from sampling and structural (engine-spread) uncertainty. This is the single
   biggest scientific upgrade.

   *Concrete design (scoped, ready to build):* the draws JSONs already exist
   (`config/calibrated/<variant>_draws.json`); `config_loader` supports
   `version="custom", custom_config=<path>`, and `run_fvs_projection` takes a
   config_version. Propagation: for each of N draws (start N=30), materialize the
   draw's parameter vector as a custom config JSON, run a per-state plot
   *subsample* (~80 plots/state is enough for the density CI) through the
   projection with that config, collect carbon. The ensemble across draws gives
   the parametric density CI per state/year; scale to population totals by the
   existing area model. Cost ~ 80 plots x 50 states x 30 draws x 2 scenarios,
   one SLURM array, comparable to one campaign arm. Only the calibrated variants
   have posteriors, so the CI attaches to the calibrated engine. Validate on one
   state (e.g. ME, 30 draws) before the full array.
4. **Forest-type stratum trends.** The membership table carries FORTYPCD; add a
   fourth aggregation scale (forest type) to the trend layer for type-specific
   trajectories.

## Tier 3 -- harvest and disturbance realism

5. **Stand-age / rotation-aware harvest.** The current managed model is a
   first-order density decrement. Couple the per-owner rotation logic the YC
   engine already uses (Industrial clearcut ~45 yr, NIPF partial 20 yr, etc.) so
   harvest responds to stand maturity, and reconcile FVS-managed against
   YC-managed per state.
6. **Resolve the disturbance temporal basis.** `p_disturbance` is a multi-year
   susceptibility probability; confirm its window (currently annualized over 20
   yr) and, if available, use the dated disturbance layers (2016/2020/2022) to
   fit an annual hazard rather than assume one.
7. **Harvested-wood-products carry-over.** Harvested carbon currently leaves the
   live pool; route `harvest_c_yr` into an HWP pool (the ecoregion-economics
   layer already models HWP/NPV) for a full managed-landscape carbon balance.

## Tier 4 -- spatial and scope

8. **TreeMap-spatial dashboard engine.** Add `fvs_*_treemap` as an explicit
   engine (like `yc_treemap_spatial_v1`) so users can toggle the area basis; the
   per-plot areas are already computed.
9. **Sub-state / county reporting.** The TreeMap pixel join supports county and
   ecoregion polygons directly -- the area basis where TreeMap most outperforms
   uniform expansion (state ratios already span 0.36 to 1.39).
10. **Climate-sensitive variant.** Current runs are climate-static (agreed
    scope). The Climate-FVS hooks or a CMIP-driven site-index shift would give a
    climate scenario axis.

## Tier 5 -- engineering hygiene

11. **Fold the height + TPA fixes upstream** into the treeinit extraction so
    future campaigns start correct (rather than patching post-hoc).
12. **NSBE height imputation in `compute_plot_agb`** as a belt-and-suspenders
    fallback, so any future missing-height trees never zero out.
13. **One-command pipeline** (treeinit fix -> campaign -> aggregate -> merge ->
    trends -> figures) as a Make/Snakemake target for reproducibility.
