# An effective strategy for integrating the CONUS calibration into fvs-modern

**Date:** 2026-06-11
**Scope:** the end to end plan to land the unified CONUS work in the fvs-modern engine: the
two-layer species-free tree model (both legs plus blend), the CSPI v7 site surface, ingrowth
and species composition, the new Garcia tradition stand level constraint layer, and
uncertainty, with a sequenced, gated rollout and a parity-tested engine path.

## 1. What fvs-modern becomes

A single calibrated variant, `conus_sf`, replaces the twenty regional variants. It has four
coordinated parts, all driven by one productivity index (BGI for growth, the CSPI surface
for the height and crown terms) and all carrying posterior draws for uncertainty:

1. **Tree level base equations**, six components (DG, HG, HT-DBH, HCB, survival, CR), each
   carrying both a species-dependent and a species-independent readout and a per-species
   shrinkage blend (`w = n/(n+kappa)`, stress-test defaults for kappa).
2. **Tree level modifiers**, the Layer 2 disturbance and treatment adjustments (common form,
   or trait-mediated for HCB and CR).
3. **Stand level constraint layer**, a compact age-independent annualized state model for top
   height, density, and basal area that disciplines the summed tree predictions.
4. **Ingrowth and composition**, a plot level recruitment count plus a trait-driven species
   composition allocation that seeds new trees.

## 2. The engine execution flow, one annual cycle (the heart of it)

For each plot, advancing from year t to t+1, for each posterior draw j:

1. **Summarize state** from the current tree list: N, BA, BAL, SDI, RD, QMD, CCH, and top
   height H, all metric.
2. **Stand level projection (constraint targets).** Advance the stand state one year with the
   transition equations to get target H(t+1), N(t+1), G(t+1). Each is a function only of the
   current state, the interval, and the productivity index, so it is annualized and path
   invariant.
3. **Tree level base prediction**, per tree: blend the two species readouts for the
   intercept, add the shared slopes, the site term, and the Layer 2 modifier, to get the
   annual DG and HG increments, the static HCB (and CR = 1 minus HCB/HT), and the annual
   survival probability.
4. **Constrain by disaggregation.** Scale the tree level diameter growth and survival so the
   summed basal area growth and the resulting density match the stand level G(t+1) and N(t+1)
   targets, using proportional allocation that preserves the tree level distribution shape
   (Ritchie and Hann). Top height H(t+1) anchors the height increments. The individual tree
   competition still decides who grows and who dies; the stand layer sets the totals.
5. **Ingrowth.** The recruitment count model gives the number of new trees; the composition
   model allocates them among species from overstory composition, site, disturbance, and
   shade-tolerance traits; seed them into the tree list at the ingrowth DBH threshold with
   sampled heights. This is the +N term that the stand level density equation expects.
6. **Advance** the tree list to t+1 and repeat.

Running steps 1 through 6 across draws j gives trajectory envelopes; the existing Bakuzis
aggregate turns them into median plus credible bands.

The new and load-bearing idea is steps 2 and 4: a small, robust stand model sets the
aggregate envelope and the tree model fills it in, so long horizon projections stay
biologically disciplined instead of drifting on compounded tree level error.

## 3. Config mapping (extends the conus_sf schema)

Each fvs-conus product maps to a block the engine reads:

| fvs-conus product | engine config block |
|---|---|
| six base bundles (`*_sf_*`) | `components[c].shared / species_independent / species_dependent / blend` |
| Layer 2 modifiers (62m) | `components[c].modifier` (common or trait_mediated) |
| CSPI v7 surface + QRF intervals | `site.surface = "cspi_v7"`, `site.uncertainty = "qrf"` |
| top height / density / basal area transitions | `stand_level.states[H,N,G]` with transition params + `disaggregation.rule` |
| ingrowth count + composition | `ingrowth.count` (negbinom or hurdle) + `ingrowth.composition` (multinomial/Dirichlet) |
| thinned posterior draws (all of the above) | `*.posterior` draw references |

The schema already carries the tree blocks; this strategy adds `stand_level`, `ingrowth`, and
a top-level `site` block, plus the disaggregation rule.

## 4. Stand level constraint mechanism, phased

- **Phase A, one way diagnostic (low risk, do first).** Run the stand transitions alongside
  the projection and report where the summed tree predictions diverge from the stand
  trajectory, using the existing `36_conus_benchmark.R` engine. No change to predictions; it
  only measures drift and validates the stand equations against FIA.
- **Phase B, two way constraint.** Turn on the disaggregation so the stand targets scale the
  tree predictions. Gate B on A: only constrain once the stand equations validate and the
  diagnostic shows the constraint reduces drift on held out plots.

## 5. Uncertainty, end to end

Two sources combined by Monte Carlo over the draw index: parameter draws for every component
(tree, modifier, stand, ingrowth) carried as thinned parquet referenced from the config, and
site draws from the CSPI v7 QRF predictive interval. The stand states H, N, G are propagated
through the same draw loop so the aggregate trajectory has credible bands, not just the tree
sums. This reuses `21_uncertainty_propagation.R` and the Bakuzis pipeline; only the source
loaders are new (`41_conus_sf_uncertainty.R`).

## 6. Parity and validation

- **Engine vs R reference parity.** The Fortran path must match the R reference predictor
  (`conus_sf_predict.R`) per tree to tolerance, for the MAP draw and several posterior draws,
  for each leg and the blend, and now also for the stand transitions and the disaggregation
  step (`test_conus_sf_parity.R`, extended).
- **Stand level validation.** Path invariance unit test (ten one-year steps equal one ten-year
  step) for each transition; held out FIA plot benchmark of H, N, G trajectories; self-thinning
  consistency against the SDImax surface.
- **System validation.** The Bakuzis biological-law checks and the FIA benchmark engine on
  full stand trajectories, calibrated vs default, with credible bands.

## 7. Sequenced milestones (gated)

1. **Lock inputs.** Finalize v7 as the site surface; finalize HG/HCB site term on the ΔLOO
   result (running). Confirm kappa defaults (done, stress-test values).
2. **Exporter.** Fill `40_export_conus_sf_config.R` with the confirmed bundle paths so the
   six tree components plus modifiers emit to `conus_sf.json` with draws. Parity test the tree
   path first, in isolation.
3. **Stand layer, fit and validate.** Top height GADA first (prototype running), then density
   tied to SDImax, then basal area. Wire as Phase A diagnostic. Gate the rest on top height
   validating.
4. **Ingrowth and composition.** Land the count (hurdle vs negbinom by held out ELPD) and the
   trait-driven composition; reconcile recruitment with the stand density term.
5. **Constrain.** Turn on Phase B disaggregation once the stand layer validates and reduces
   drift.
6. **Uncertainty and release.** Wire the draws and the QRF site intervals, run the Bakuzis
   credible-band validation, then production config write, GitHub merge of the feature branch,
   Zenodo new version, DOI backfill.

## 8. Decisions still open

- HG and HCB site term: pending the ΔLOO (running). The schema flips `site.cspi` per component
  when it lands.
- Disaggregation target: constrain on basal area growth and density (recommended) versus
  volume; basal area is the cleaner observable and ties to SDImax.
- Draw count vs engine runtime: 200 explicit draws to start, Cholesky tier if too slow.
- Ingrowth base: hurdle versus negative binomial, settled by held out ELPD.

## 9. Why this is the effective path

It reuses what exists (the calibration framework, the posterior-to-json pipeline, the
uncertainty propagation, the Bakuzis and FIA benchmark engines) rather than rebuilding, it
keeps one consistent productivity driver and one species architecture across all levels, it
adds the stand layer as a low-risk diagnostic before any constraint, and it gates every step
on a validation so the integration never advances on an unconfirmed piece. The result is a
single maintainable variant with disciplined long-horizon behavior and credible bands, which
is exactly what the carbon and yield applications downstream consume.
