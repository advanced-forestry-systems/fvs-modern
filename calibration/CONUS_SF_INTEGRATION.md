# Integrating the CONUS species-free equations into fvs-modern, with uncertainty and a species blend

**Date:** 2026-06-11
**Author:** drafted on autopilot for Aaron Weiskittel
**Inputs:** fvs-conus bundles + modifier writer (handoffs 06-05, 06-09), the fvs-modern
calibration framework (`CONUS_MODEL_SPECIFICATION.md`, `06_posterior_to_json.R`,
`21_uncertainty_propagation.R`, `BAKUZIS_UNCERTAINTY_README.md`), and Aaron's directive
to carry **both** species-dependent and species-independent equations with a per-species
blend.

## 1. What already exists on each side

**fvs-conus produces (the new unified CONUS model):**
six base components, each as a Stan posterior with a standard bundle
(`_sf_fixed`, `_sf_gamma`, `_sf_species`, `_sf_re_L1/L2/L3/FT`, `_sf_manifest`):
DG (v8 BGI), HG (v8rd), HT-DBH, HCB, Survival (crz), CR (CR2-direct). Plus a Layer 2
modifier writer (`62m`) that emits disturbance/treatment modifiers per component
(common form, or trait-mediated for HCB and CR), and the CSPI site surfaces
(v4 production at 1 km; v6 best accuracy; v7/v8 distilled to 30 m, v7 QRF in progress
for prediction intervals).

**fvs-modern consumes (the engine side):**
calibrated parameters are turned into engine JSON by `06_posterior_to_json.R`, the
Fortran engine in `src-converted/` reads that config, a Python shared-library bridge
(`fvs2py`) drives projections, and uncertainty is already propagated by
`21_uncertainty_propagation.R` (sample N posterior draws, run a projection per draw,
take percentile bands) and the Bakuzis credible-band pipeline (default + calibrated MAP
+ N draws over a 36-scenario grid). The species threshold in the spec is **>= 5,000
tree-period observations** for species-level estimation, with species below threshold
**pooled into groups**. That threshold is precisely the species-dependent vs
species-independent boundary, and "pool below threshold" is exactly where the blend goes.

The uncertainty machinery exists but is wired to the older per-variant species-dependent
calibration (DG, mortality, H-D). The job is to feed the new unified species-free
posteriors and modifiers through the same path, and to make both legs first-class.

## 2. The unifying idea: the trait term is the prior mean, the species delta is the data

The Stan models already use non-centered hierarchical intercepts,
`b0[s] = mu_b0 + sigma_b0 * z_b0[s]`. The species-free (Leg B) work replaces the scalar
`mu_b0` with a trait predictor, `mu_b0(s) = W[s] %*% gamma`. So the two legs are not rival
models, they are two readouts of one fit:

- **Species-independent (Leg B, `pure_sf`):** use `W[s] %*% gamma` alone. Defined for
  every species, including the long tail with no data.
- **Species-dependent (Leg A):** use the fitted per-species intercept `b0[s]` where the
  species is well sampled.
- **Blend (the `hybrid_source_map`):** `b0_used(s) = w_s * b0[s] + (1 - w_s) * (W[s] %*% gamma)`.

The weight `w_s` is the natural shrinkage weight, `w_s = n_s / (n_s + kappa)`, so a
well-sampled species leans on its own intercept and a rare species falls back to traits.
Setting `kappa -> 0` recovers pure Leg A where data exist; `kappa -> inf` recovers pure
Leg B; the stress test says the honest default is near the Leg B end for HG/HCB/HT-DBH and
mortality, and meaningfully toward Leg A for CR. This single weight reproduces every arm of
the 06-05 stress test rather than hard-coding a 0/1 source map, and it degrades gracefully.

This is the design the engine should carry: store both readouts plus `w_s`, blend at
prediction time.

## 3. Config schema (engine-facing)

Extend the calibrated-variant JSON with one block per component holding three parts:

```jsonc
"diameter_growth": {
  "form": "organon_swo_v8_bgi",
  "shared": { "b1": ..., "b2": ..., "K1": ..., "K2": ..., "...": ... },  // fixed slopes
  "species_independent": {                 // Leg B: trait term
    "trait_cols": ["wood_density","shade_tol","..."],
    "gamma":      [ ... ],                  // W %*% gamma gives the intercept
    "scale_mean": [ ... ], "scale_sd": [ ... ]
  },
  "species_dependent": {                    // Leg A: per-species fitted intercepts
    "spcd":  [131, 202, ...],
    "b0":    [ ... ],
    "n_obs": [ ... ]                        // drives the blend weight
  },
  "blend": { "rule": "shrinkage", "kappa": 1500, "min_n_legA": 500 },
  "site_term": { "source": "BGI", "cspi_version": "v4" },   // per the v6 screen
  "modifier": { "form": "common", "alpha_0": ..., "alpha": { "insect": ..., "...": ... } },
  "posterior": { "n_draws": 200, "path": "draws/diameter_growth_draws.parquet" }
}
```

Notes that fall straight out of the existing results:

- `site_term` per component matches the 06-09 CSPI v6 screen: DG and HT-DBH keep their
  current term (BGI, CSPI v4); HG and HCB switch to CSPI v6 **iff** the resubmitted
  ΔLOO (jobs 11505362 / 11505363) confirms it. The schema carries a `cspi_version` field
  so the surface is a swap, not a refit, once confirmed.
- `modifier` is the Layer 2 block from `62m`: `common` for DG/HG/ingrowth/mortality,
  `trait_mediated` (global alphas + gamma file, subset of disturbance types) for HCB and CR.
- `blend.kappa` is a single tunable; CR gets a smaller kappa (leans species-dependent),
  the rest a larger one (leans species-free). Defaults from the stress test.

## 4. Uncertainty quantification

Two independent sources, combined by Monte Carlo over draws:

1. **Parameter uncertainty.** Carry a thinned set of posterior draws (suggest 200) for
   each component: shared slopes, trait gammas, per-species z, and modifier alphas/gammas.
   Store as a compact `draws/*.parquet` referenced from the config. The engine, for draw
   `j`, applies coefficient set `j` to the whole projection. This is exactly what
   `21_uncertainty_propagation.R` already does for DG/mortality/H-D; the extension is to
   point it at the unified species-free draws and the modifier draws.
2. **Site-surface uncertainty.** The CSPI prediction is itself uncertain. The v7 QRF model
   (job 11503014, running now) gives a predictive interval per pixel. Draw the site index
   from that interval per plot per draw, so site error propagates alongside coefficient
   error rather than being treated as known.

Combined, a projection run over `j = 1..N` draws (coefficients `j` plus a site draw `j`)
yields trajectory envelopes. The Bakuzis aggregate and figure scripts already turn that
long CSV into median plus 95 percent bands, so the reporting end is reuse, not new code.

A cheaper tier, for when full draw propagation is too expensive in the engine: store the
mean and a Cholesky factor of the coefficient covariance per component and sample inside
the engine. Recommend starting with the explicit thinned draws (simpler, exact) and adding
the Cholesky path only if runtime demands it.

## 5. Concrete pieces to build (the scaffold)

1. **Exporter** `40_export_conus_sf_config.R` (companion to `06_posterior_to_json.R`):
   reads the six bundles + the modifier JSONs + the CSPI metadata, assembles the schema in
   section 3 for all components, writes `config/calibrated/conus_sf.json` plus the thinned
   `draws/*.parquet`. Emits both legs and the blend block.
2. **Blend + predict reference** `conus_sf_predict.R` / engine hook: given a tree row and a
   draw index, compute `base = blend(legA, legB)`, apply the modifier, derive CR from HCB.
   This is the reference implementation the Fortran path must match (parity test against it).
3. **Uncertainty extension** `41_conus_sf_uncertainty.R`: thin wrapper that registers the
   unified species-free draws + site draws with `21_uncertainty_propagation.R` and the
   Bakuzis pipeline.
4. **Parity test** `test_conus_sf_parity.R`: assert the Fortran engine and the R reference
   agree per tree to tolerance for a sample of FIA plots, for the MAP draw and a few
   posterior draws, for each leg and the blend.

The exporter and the predict reference are the load-bearing two; 3 and 4 reuse existing
machinery and guard the integration.

## 6. Open decisions for Aaron

- **kappa per component.** The stress test sets the direction (CR low, others high). Do you
  want a fitted kappa (minimize held-out error of the blend) or the stress-test defaults?
- **CSPI surface for production.** v6 (1 km, R2 0.75) vs v7/v8 (30 m, R2 0.69). The schema
  supports either via `cspi_version`; the v6 ΔLOO result decides whether HG/HCB move to v6.
- **Draw count vs runtime.** 200 draws is a starting point for the bands; the engine cost
  scales linearly. Confirm acceptable, or set the Cholesky tier.
- **Repo home.** Exporter and config live in fvs-modern `calibration/`; the bundles stay in
  fvs-conus. Confirm the scaffold should land as a draft PR on `feat/conus-sf-integration`.

## 7. Status of the unblock

The 06-09 v6 confirmatory ΔLOO jobs failed only because the 40k subset could not clear the
hard-coded 5,000-obs species threshold. Fixed: the HG driver now takes a `min_obs_species`
argument (default unchanged at 5,000), the HCB driver already had `min_sp`, and a fresh
150k subset drawn from the full v6-joined file feeds both. Resubmitted at threshold 500
(HG 43 species, HCB 57 species): jobs **11505362** (HG) and **11505363** (HCB). Their
elpd_diff vs SE decides whether `site_term.cspi_version` flips to v6 for those two
components. No production config, coefficients, or git history were touched.
