# Section 6 (v2). Conclusion

*Draft v2 Conclusion. Replaces v1 Section 6 in its entirety. Anchored
on the three concrete contributions outlined in the v2 Introduction
(species-free traits, data-driven covariates, dual-architecture
runtime) and on the operational implications laid out in Section 5.*

---

The fvs-modern v2 release advances the FVS calibration program along
three concrete dimensions while preserving the open-source release model
and continuous integration infrastructure that v1 introduced.

The species-free trait-driven reparameterization is the most consequential
of the three. By replacing per-species random intercepts with a shared
trait-coefficient vector, we preserve residual variance to within the
prespecified 5 percent functional-equivalence threshold (HG difference
0.26 percent; DG Kuehne difference 3.98 percent) while sharpening
identification of the trait coefficients themselves (DG Kuehne: 8 of 8
γ exclude zero under the species-free fit vs 4 of 8 under the species
random-intercept reference). The signs and orderings of the trait
coefficients are robust to the architecture choice, so the biological
story does not depend on which parameterization is used. The new
architecture extrapolates cleanly to species that fall outside the
training set provided trait values are available, which addresses one
of the persistent practical limitations of growth and yield models in
the context of climate-driven range shifts and afforestation projects.

The data-driven covariate selection is the second contribution. The
covariate exploration on a 100,000-observation subsample identified
the FIA Site Index (SICOND) as substantially stronger than the climate-
derived site index that v1 had selected on biological priors, and
identified stand age, plantation status, and elevation as additional
covariates with substantive predictive contribution. The exploration
also flagged that the apparent negative coefficient on the linear
climate_si term under v4 was a variable-substitution artifact rather
than a real biological signal. The lesson generalizes: targeted
covariate exploration on a small subsample is computationally cheap
relative to a single Stan fit and can shift methodological direction
in ways that the fits themselves cannot. We recommend that the practice
become standard in operational Bayesian calibration of growth and
yield models.

The dual-architecture runtime is the third contribution. The
`config_loader.py` interface added two new parameter sources
(`ParameterSource.CONUS` for strict CONUS-Phase-4 read; `ParameterSource.HYBRID`
for per-component fallback to the per-variant fits) and the
`get_conus_runtime_block()` decomposer translates the new categories_conus
JSON block into the runtime form expected by the FVS Fortran shim. The
two architectures coexist in the same variant JSON, so applications can
switch between them through a configuration flag without re-fitting. We
expect the species-free architecture to become the default for
applications that prioritize interpretability and species-out
extrapolation, with the per-variant architecture remaining preferred for
applications that need the marginal calibration accuracy of per-species
intercepts within established training-set species.

Three threads remain open as the v2 manuscript ships. The full-data B1
height growth refit on all 1.98 million CONUS observations is queued and
will sharpen the trait coefficients beyond the pilot-scale evidence
reported here. The data-driven v6 production candidate
(`dg_kuehne2022_speciesfree_v6_data_driven.stan`) is built but not yet
fit; selection between v6 and the v4 architecture variant suite will be
by leave-one-out cross-validation once both are in hand. And the
held-out species predictive evaluation, which will quantify the
species-out extrapolation accuracy of the trait-driven architecture, is
in the evaluation pipeline. The v3 manuscript will report on these three
items.

The broader implication for the FVS community is that the operational
parameter set need not be a fixed point. Each of the three contributions
in this release emerged from analytical decisions made after two years
of operational use of v1; future operational use will continue to
reveal the methodological choices that matter. fvs-modern is structured
to absorb those choices: the open-source codebase, the runtime
configuration switching, the continuous integration matrix across 25
variants, and the documented integration plan from CONUS Phase 4 fits
through to per-variant JSON blocks together provide the infrastructure
for an FVS that improves continuously as the underlying science and the
underlying data evolve.

---

## Notes for the Conclusion finalization round

- The "we recommend the practice become standard" sentence in the
  data-driven-covariate paragraph is a strong claim; consider softening
  to "we suggest..." or adding a citation to similar empirical
  covariate-selection practice in related fields.
- The "dual-architecture runtime" framing in paragraph 4 may need to be
  reconciled with what the journal's audience considers a "release"
  versus a "research result." The runtime work is concrete software
  engineering; the trait-driven and covariate-selection work are
  research results. Three contributions of mixed kinds may need
  rebalancing in the abstract and intro to match the journal target.
- "We expect the species-free architecture to become the default"
  is a forward-looking statement; some journals will want this softened
  or attached to a specific timeline.
- The final paragraph echoes the v1 manuscript's positioning vis-a-vis
  Premer et al. 2025 and Woodall et al. 2025; check whether to maintain
  parallel framing.
