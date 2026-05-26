# Abstract (v2 draft)

*Draft for the v2 manuscript. ~250 words. Anchored on the species-free
trait-driven pivot result (HG and DG Kuehne) and the data-driven covariate
selection result (SICOND replaces climate_si). Two `[CITE: ...]` markers
remain for the Paperpile pass.*

---

The Forest Vegetation Simulator (FVS) is the dominant operational growth
and yield modeling system in U.S. forestry, but its parameterization has
relied on species random intercepts and a small set of climate covariates
chosen on biological priors rather than empirical evidence. We recalibrate
six FVS component submodels (diameter growth, height growth, height to
diameter, height to crown base, mortality, and crown recession) at CONUS
scale using hierarchical Bayesian methods on 433,291 FIA condition pairs,
and we test two methodological refinements that depart from the standard
species random-intercept parameterization.

First, we replace species random intercepts with a linear combination of
eight standardized species traits (wood density, specific leaf area, leaf
nitrogen, seed mass, maximum height, shade tolerance, drought tolerance,
cold tolerance). This species-free formulation produces residual standard
deviations statistically indistinguishable from the species random
intercept reference fit (HG difference 0.26%; DG Kuehne difference 3.98%,
both within the prespecified 5% functional-equivalence threshold) while
sharpening identification of the trait coefficients themselves (DG Kuehne:
8 of 8 trait coefficients exclude zero under the species-free fit versus
4 of 8 under the reference). The trait-coefficient signs and orderings
agree between architectures where both are well identified.

Second, an empirical covariate exploration on a 100,000-observation
subsample identifies the FIA Site Index (SICOND) as a substantially
stronger site-productivity predictor than the climate-derived site index
(SICOND univariate R²=0.044 versus climate_si R²=0.007), with stand age,
plantation status, and elevation jointly contributing 30% additional
explained variance beyond core tree and stand covariates. The
data-driven covariate set replaces the prior climate-only specification
in the production CONUS diameter growth model.

Together these refinements deliver a more interpretable, more
identifiable, and more transferable CONUS-scale FVS recalibration. The
species-free formulation makes the model extrapolate cleanly to species
not in the training set provided trait values are available; the
data-driven covariate selection grounds the model in empirically
demonstrated rather than assumed predictive structure. The complete
recalibration pipeline, posterior summaries, and integration tooling are
released as open-source software (fvs-modern, MIT licensed) with
continuous integration and DOI minting via Zenodo
[CITE: 10.5281/zenodo.19802673].

---

## Notes for the Discussion / Abstract finalization round

- Word count is at the upper bound; tighten by 20-30 words for typical
  journal limits.
- The 5 percent functional-equivalence threshold needs a citation
  ([CITE: pre-specified threshold reference]) — same marker that
  appears in `v2_section5_discussion.md`.
- Add Premer et al. 2025 / Joo et al. 2025 / Woodall et al. 2025 cross-
  references if the journal target wants positioning against the
  comparative-modeling and forest-carbon-uncertainty literature.
- The "more interpretable, more identifiable, more transferable" triad
  in the final paragraph is the v2 manuscript's headline claim and
  should be echoed in the Introduction Section 1.
- The HG full-data B1 refit (1.98M observations) is referenced in
  Section 5.7 but is not yet in hand. Update the Abstract to mention
  the production-scale result when the fit completes.
