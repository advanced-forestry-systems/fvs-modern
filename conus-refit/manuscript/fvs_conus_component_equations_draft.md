# A unified Bayesian hierarchical family of forest growth component equations for the conterminous United States

Draft scaffold, 2026-06-19. Companion to the fvs-modern engine/calibration paper. This paper is the
component-equation methods paper; the engine paper covers the modernized simulator, variant keyword
recalibration, and deployment. The two are designed to be non-overlapping: this one introduces and evaluates
the equation forms; the engine paper consumes their predictions.

## Scope and contribution (distinct from the engine paper)

The engine paper (fvs-modern) recalibrates the existing FVS component models at the keyword level (per-species
basal-area-increment multipliers, max-SDI, density-dependent recruitment) and reports stand-level benchmark and
international validation. This paper introduces a NEW family of component equations fit from scratch on national
FIA remeasurement data: diameter growth, height-diameter, crown ratio, height growth, survival, and ingrowth,
each as a Bayesian hierarchical model with a common random-effects structure (species, measured functional
traits, forest-type group, and EPA ecoregion at Levels 1 to 3). The central scientific question is whether a
species-free formulation, in which species enters only through traits and the hierarchical grouping, matches a
species-specific formulation, which would let the equations generalize to species and regions with sparse data.

## Proposed structure

1. Introduction. Limits of legacy variant-specific component equations; the case for a unified,
   trait-informed, hierarchically pooled formulation; the species-specific versus species-free question.
2. Data. National FIA remeasurement pairs (CONUS), the condition-level disturbance screen, functional-trait
   covariates, the CSPI site-productivity covariate (the index itself is documented in the companion CSPI
   deposit and paper), ecoregion and forest-type-group assignment.
3. Model family. The shared hierarchical template: non-centered random intercepts for species (z_sp),
   ecoregion (z_L1/L2/L3), and forest-type group (z_FT); trait effects (W gamma); component-specific likelihoods
   (diameter growth, height-diameter via algebraic difference, crown ratio, height growth, survival as
   Bernoulli, ingrowth as negative binomial). Both species-specific and species-free instantiations.
4. Fitting and evaluation. cmdstanr, convergence diagnostics, county-hash spatially blocked held-out folds,
   leave-one-out and held-out ELPD comparing species-specific versus species-free.
5. Results.
   - Component fits and trait effects.
   - Species-specific versus species-free: the per-species deviation (sigma_sp, response scale) is well
     identified for mortality (about +53 percent odds), diameter growth (about +23 percent), and crown
     recession (about +29 percent), and weak for height growth, height-diameter, and height-to-crown-base.
   - Stand-level consequence (from the engine four-arm comparison): the species-specific growth signal gives no
     median stand-level gain, so species-free is competitive for stand projection while species-specific is
     retained for mortality, diameter growth, and crown recession.
   - Forest-type and ecoregion dependence of all components (both species-dependent and species-free forms now
     carry z_FT and z_L1/L2/L3).
6. Discussion. When species-free suffices and when species-specific is needed; transferability to data-sparse
   species and regions; coupling to the FVS engine (the companion paper).

## Material to assemble (already drafted in the repo)

- 20260511_manuscript_methods_results_v1.md (methods and first results)
- 20260513_manuscript_v2_scope.md (scope refinement)
- 20260513_section_3.7_draft.md (a results subsection)
- CONUS_MODEL_SPECIFICATION.md (the equation specification of record)
- BAKUZIS_UNCERTAINTY_README.md (uncertainty framing)
- Stan specifications in stan/ (including the new dg_speciesfree_v8_forest_eco and
  ingrowth_negbinom_v3_forest_eco added this program, which give every component both forest-type and
  ecoregion structure)

## Next drafting steps (autopilot queue)

1. Merge the three existing manuscript markdown fragments into a single ordered draft against the structure
   above, reconciling notation to CONUS_MODEL_SPECIFICATION.md.
2. Produce the species-specific versus species-free ELPD table and the sigma_sp forest plot as the headline
   result figure.
3. Pull the component fit tables from output/comparisons and output/variants.
4. Run the citation-integrity and writing-voice passes (manuscript-preparer) before the docx build.
