# Section 1 (v2). Introduction

*Draft v2 Introduction. The first three paragraphs are preserved from
v1 with minimal edits. The fourth paragraph (Itter / Joo / Premer /
Woodall positioning) is extended to set up the v2 methodological
contributions. The final paragraph is rewritten to enumerate v2's
three concrete additions over v1.*

---

The Forest Vegetation Simulator (FVS) has been the operational standard
for individual tree growth and yield projection across federally managed
forests in the United States for more than four decades
[CITE: Dixon 2002]. FVS is organized into 25 geographic variants, each
with species-specific growth, mortality, height-diameter, and crown
ratio equations parameterized against regional forest data. Over time,
variant-specific parameters have accumulated from diverse historical
calibrations, often using inventory data collected prior to 2005
[CITE: Weiskittel et al. 2011], and the codebase has evolved from its
Fortran 77 origins with limited structured testing or dependency
management. Two modernization gaps follow from this history.

First, the source code is difficult for the research community to
extend. The architectural review by Crookston and Dixon (2005) documents
a codebase organized around a shared base engine and 25 geographic
variants but rooted in fixed-form Fortran 77 with extensive COMMON block
usage and computed GOTO branching. Two decades of subsequent maintenance
have preserved that organization without modernizing it: approximately
2,054 source files remain in fixed form, with limited unit testing, no
continuous integration, and few entry points for non-Fortran programmers.
Forty years of legacy code have therefore accumulated alongside the
variants, and no released variant carries a built-in mechanism for
propagating parameter uncertainty through projections. These
characteristics inhibit reproducibility, community contribution, and
platform portability.

Second, the parameters are stale. Most default FVS variant parameters
were fit to data collected before the modern Forest Inventory and
Analysis (FIA) program adopted its annual panel design
[CITE: Bechtold and Patterson 2005]. Advances in Bayesian hierarchical
modeling [CITE: Gelman et al. 2013] and the availability of roughly two
decades of additional FIA remeasurements create an opportunity to refit
all variants with modern data and methods, and to do so in a way that
preserves full posterior distributions for downstream uncertainty
propagation. This second gap extends beyond the parameter values
themselves to the architectural choices that govern how the parameters
enter the equations. The standard FVS parameterization carries one
random intercept per species per variant per component; with 25 variants
and up to 100+ species per variant across six component models, the
parameter count balloons quickly and many of the per-species intercepts
are weakly identified for species with few observations. A
species-trait-driven reparameterization that pools information across
species through a shared trait-coefficient vector offers a path to
better identifiability and to clean species-out extrapolation.

Recent individual tree growth models have increasingly adopted Bayesian
methods for calibration [CITE: Clark et al. 2016; Bohn and Huth 2018],
and Itter, Finley, and Weiskittel (2025) recently formalized how
connecting growth and yield models to continuous forest inventory data
with hierarchical Bayesian methods enables principled uncertainty
accounting at the stand and regional scales. The Pacific Northwest
perspective by Joo, Temesgen, Frank, Weiskittel, and Reimer (2025)
further documents the operational stakes: site index estimation,
maximum stand density derivation, and error propagation in long-horizon
projections are the practical limitations that growth and yield models
must address to remain credible for management decisions. The recent
comparative review by Premer, Simons-Legaard, Daigneault, Hayes, Solarik,
and Weiskittel (2025) shows that projected forest carbon outcomes can
differ by a factor of two across modeling systems applied to the same
Maine plots, and the Forest Carbon Modeling Group statement
[CITE: Woodall et al. 2025] identifies model uncertainty quantification
as a top community priority for the next decade. Three operational
implications follow. First, parameter uncertainty must be propagated to
projections rather than presented as point estimates. Second, calibration
choices must be auditable so that downstream comparisons can attribute
differences to specific methodological decisions. Third, the
parameterization must extrapolate cleanly to species that fall outside
the original training set, both for new species entering historic
ranges as climate shifts and for application of regional fits to species
that were sparsely represented in the source data. No prior effort has
combined a comprehensive Bayesian recalibration of FVS with an open
source software release, a runtime uncertainty engine, a factorial
benchmark spanning long horizons, and a trait-driven reparameterization
that targets all three of these operational requirements.

This paper describes the integrated release, fvs-modern, and the v2
methodological refinements that extend it. Section 2 documents the
software architecture, Section 3 the calibration methodology and
benchmarking framework, Section 4 the headline results from the
v1 per-variant calibration, and Section 5 the v2 methodological
refinements that emerged from two years of operational use. Specifically,
v2 contributes three additions over the v1 manuscript. First, the
species-random-intercept architecture used throughout v1 is replaced
with a species-trait linear predictor; we validate the substitution on
two CONUS Phase 4 fits with different forms and dataset sizes (height
growth and Kuehne 2022 diameter growth) and find that the residual
standard deviation is preserved while the trait coefficients are
estimated more sharply. Second, the site productivity covariate set
is revised on empirical rather than biological-prior grounds: a targeted
covariate exploration identifies the FIA Site Index (SICOND) as
substantially stronger than the climate-derived site index that v1 used,
identifies plantation status as deserving main-effect representation,
and adds elevation as a topographic covariate. Third, the production
runtime is updated to support both the per-variant and the
species-free architectures through a single configuration switch, so
that downstream applications can choose the architecture appropriate to
their species coverage and uncertainty requirements without re-fitting.

We frame these three additions as the operational consequence of two
years of FVS-modern in production. Each is a methodological choice that
revealed itself through analysis rather than through the original
calibration design: the species-RE limitation surfaced when calibrated
fits extrapolated awkwardly to species not in the training set; the
climate_si limitation surfaced when the v4 architecture variant suite
produced an implausibly negative coefficient on the linear site term.
Section 5.6 develops the species-free result; Section 5.8 develops the
covariate-exploration result; Section 5.7 lays out the remaining open
threads for the v3 manuscript.
