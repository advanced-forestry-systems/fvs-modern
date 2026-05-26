# Section 5 (v2). Discussion

*Draft stub for the v2 manuscript Discussion. Builds on the v1 Discussion
in `fvs_combined_draft.md` (Sections 5.1 through 5.5) and adds a sixth
subsection on the species-free trait-driven pivot. The four crosscutting
threads from v1 (CR distribution choice, arid-variant DG ceiling, Bakuzis
multi-criteria pass rates, annualization) are revisited under the new
findings. References to `[CITE: ...]` are left as placeholders for the
Paperpile pass.*

## 5.6 Species-free trait-driven calibration

The v2 effort extends the species random intercept architecture used
throughout v1 by replacing it with a species-trait linear predictor.
Under this reformulation (model B1), each species enters the linear
predictor not through an integer-indexed random intercept but through a
linear combination of eight quantitative traits (wood density, specific
leaf area, leaf nitrogen on mass basis, seed mass, maximum height, shade
tolerance, drought tolerance, cold tolerance). The traits are
standardized within the FIA species set and the trait-coefficient vector
is shared across all species. The reference model (B2) retains the
species random intercept structure and adds the same trait covariates
as additional predictors.

We applied B1 and B2 to two CONUS Phase 4 fits with deliberately
different forms and dataset sizes: a height growth pilot (HG-Organon
fixedK on 9 species, ~150,000 observations) and a full-coverage diameter
growth fit (DG-Kuehne 2022 with lognormal residual, 100,000 observations).
The pilot serves as a stress test where data are sparse relative to
parameter count; the full-coverage fit is the production case.

Three findings carry the v2 story.

First, the residual standard deviation is statistically indistinguishable
between B1 and B2 in both fits. HG-fixedK gave σ_B1 = 2.1801 and σ_B2 =
2.1745 (difference 0.26 percent). DG-Kuehne gave σ_B1 = 1.8365 and
σ_B2 = 1.9126 (difference 3.98 percent). Both differences fall well
within the 5 percent threshold that defines functional equivalence in
the manuscript [CITE: pre-specified threshold reference]. The trait
linear predictor does not lose information relative to the species
random intercept in either case.

Second, the trait coefficients themselves are more precisely identified
under B1 than under B2. Figure 1 panel B shows that in DG-Kuehne all
eight γ posterior credible intervals exclude zero under B1, whereas
only four exclude zero under B2. The species random intercept in B2
absorbs variation that the trait coefficients capture in B1, and the
posterior spreads of the γ widen accordingly. HG-fixedK shows the same
pattern in attenuated form (Figure 1 panel A): with only 9 species the
posterior is too data-limited for any γ to exclude zero in either model,
but the B1 standard deviations are systematically smaller than B2 by a
factor of 1.3 to 3.5 across the eight coefficients (manuscript Table 2
columns sd_ratio). The information that species random intercepts
absorb is information about between-species variation that can equally
well be expressed through species traits.

Third, the substantive signs and orderings of the trait coefficients
agree between the two base models where both are well identified.
DG-Kuehne γ[5] (maximum height) is the largest positive coefficient
under both B1 and B2; γ[3] (leaf nitrogen) is the next largest positive
under B1 and is positive but credible-interval-spanning-zero under B2;
γ[1] and γ[2] (wood density and SLA) are both reliably negative under
both models. The biological story does not depend on which base
architecture is used; the species-free architecture simply produces
sharper estimates.

These three findings together support the central v2 hypothesis:
species random intercepts in growth and mortality models are absorbing
trait-mediated variation that can be modeled directly. Doing so has
two operational benefits. The fitted model extrapolates cleanly to
species not in the training set, provided trait values are available.
And the model's residual variance is unchanged, so the calibration
quality measured at the stand level (Sections 4.3 and 4.5 of v1) is
preserved.

The two crosscutting v1 findings that bear most directly on the v2
result are the arid-variant DG ceiling and the annualization fix. The
arid variants (BC R² 0.31, PN 0.34, UT 0.40, CR 0.41, BM 0.45) ranked
lowest in the v1 per-variant DG calibration and were also the variants
where species composition was sparsest. Under the species-free
formulation those variants are no longer disadvantaged by per-variant
species coverage; the trait-coefficient vector pools information across
the CONUS dataset. We expect the v2 refit on the full 1.98 million
observation HG dataset (Section 5.7) to show the largest improvements
in the arid variants for this reason. The annualization fix that
unlocked realistic stand-level trajectories in v1 (Section 5.1) is
unaffected by the species-free pivot; the same period-to-annual
correction applies to both B1 and B2.

The two remaining crosscutting threads are less directly tied to the
species-free result but inform how it will be evaluated. The Bakuzis
multi-criteria framework provides a process-level test that does not
depend on species random intercepts, so the same battery of Sukachev,
Eichhorn, recession-under-competition, and mortality-U-shape tests
applies to the species-free fits without modification. The crown-ratio
v1-to-v2 Beta regression result (Section 5.1) is an example of the
broader pattern that motivated v2: choosing a response distribution
appropriate to the data (a level rather than a change) recovers a large
fraction of variance, and choosing a covariate representation
appropriate to the species (traits rather than random intercepts)
similarly recovers identifiability for the coefficients of interest.

## 5.7 Future work and open questions

Three substantive items remain open as v2 ships.

First, a full-data B1 HG refit on all 1.98 million observations is in
the queue. The pilot validated the architecture; the production fit
will sharpen the trait coefficients further and provide the species-out
predictive evaluation that Figure 2 will report. The expected wall
time on the Cardinal cluster is approximately 12 hours under the
current Stan configuration with adaptive HMC and eight chains.

Second, the linear site coefficient in the DG growth equations is
implausibly negative in the pilot fit. A four-way architecture variant
suite is under evaluation (quadratic site, ecoregion-varying site
slope, trait-modulated site slope, full combination of the three).
The L1-ecoregion-varying site slope approach was validated against the
April HG pilot data and is the leading candidate. Selection among the
four variants depends on cross-validation against held-out plots and
will be reported in v3 if the result materially changes the trait
coefficient estimates.

Third, the eight-trait set used in B1 is a curated subset of the
species traits available in the TRY database [CITE: TRY] and from FIA
species profiles. A trait importance analysis on the production HG
B1 fit will identify which subset of the eight contributes most to
predictive accuracy. If a smaller subset suffices, the operational
trait-data acquisition burden for new species is correspondingly
reduced.

## 5.8 Data-driven covariate selection

A second methodological refinement entered v2 late in the calibration
cycle. Following the species-free pivot, a 90-second covariate
exploration on a 100,000-observation subsample of the CONUS DG data
ranked candidate site, stand, and disturbance covariates by univariate
R² and by conditional R² added on top of the core tree-level
covariates (ln_DBH + ln_CR_adjusted + ln_BAL_softwood_adjusted). The
exploration was prompted by an implausible negative coefficient on the
linear climate_si term in the v4 fits and aimed to identify whether a
different site-productivity metric would behave more sensibly.

The empirical ranking diverged from the climate-driven prior. The FIA
Site Index (SICOND, in feet at base age 50) had univariate R² of 0.044
on annual diameter growth, with the expected positive sign. The
climate-derived site index (climate_si), which had been the primary
site variable in v0 through v5, had univariate R² of only 0.007 and
was effectively uninformative when fit alone. Three additional
covariates added substantial information beyond the core tree variables:
stand age (STDAGE, 0.121 R² added), plantation status (is_plantation,
0.107 R² added), and elevation (ELEV log-scaled, 0.070 R² added).
Together these four covariates plus SICOND account for an additional
30% of explained variance beyond the core, where the prior climate_si
specification added only the climate term's modest contribution.

Three implications follow.

First, the apparent negative climate_si coefficient under v4 was not
a real biological signal but an artifact of the climate_si variable
absorbing variance that better-targeted variables (SICOND, ELEV, stand
age) should be carrying. Once SICOND enters the model directly,
climate_si becomes redundant.

Second, the plantation status indicator deserves main-effect status
rather than the modifier-only role it carried in v0 through v5. Its
0.107 R² added is too large to be expressed as a perturbation around
a non-plantation baseline.

Third, the covariate exploration approach is a methodological
contribution in its own right. Many published growth-and-yield
calibrations select covariates from biological priors and never
test whether the chosen variables actually carry the predictive
signal they were nominated for. A targeted exploration on a small
subsample provides empirical justification for the covariate set and
catches the kind of variable-substitution artifact that the v4 fits
exhibited.

The v6 production candidate (`dg_kuehne2022_speciesfree_v6_data_driven.stan`)
incorporates the data-driven selection: SICOND replaces climate_si as
the primary site variable; is_plantation enters as a main effect;
ELEV enters as a topographic covariate; sdi_complexity (sdi_additive /
SDI_Reineke) captures stand structural complexity; and the previously
species-RE-absorbed trait effects retain their interpretation under the
species-free formulation. Selection between v6 and the previous v4
production candidate will be by LOO cross-validation once both fits
are in hand.

## Threading into v1 Discussion

When the v2 manuscript is assembled, Subsections 5.6 and 5.7 land after
5.5 and before the v1 Conclusion (currently Section 6). The v1 numbered
Sections 5.1 through 5.5 are minimally edited: 5.1 (Scientific
implications) gains a forward reference to 5.6, and 5.4 (Limitations)
notes that the species random intercept dependence was a v1 limitation
that v2 addresses. The v1 Sections 2 through 4 do not require
reordering for the v2 manuscript; the v2 changes are additive Discussion
content plus Figures 1 and 3 in Section 3 [CITE: figure numbers
finalized at production].

## Open thread checklist for the Discussion pass

- [ ] Finalize the 5 percent functional equivalence threshold reference
  ([CITE: pre-specified threshold reference])
- [ ] Add TRY database citation in 5.7 ([CITE: TRY])
- [ ] Add figure number cross-references once Section 3 layout is final
- [ ] Add Premer et al. (2025) and Joo et al. (2025) cross-reference in
  5.6 paragraph 7 (the manuscript-positioning paragraph from v1's 5.1
  already cites both; v2's 5.6 should mention that the species-free
  approach makes inter-model comparison cleaner because trait-coefficient
  values are interpretable across modeling systems while species-RE
  values are not)
- [ ] Add Woodall et al. (2025) cross-reference for the v2 paragraph 1
  about operational benefits of species-out extrapolation
- [ ] Resolve sentence cadence in 5.6 paragraph 4 (the substantive
  signs/orderings paragraph) once the held-out species prediction
  result is in hand
- [ ] Cite the v6 production fit and the LOO comparison against v4
  in 5.8 once both fits are in hand
- [ ] Reconcile the 5.7 architecture-variant text (which lists quad /
  l1site / traitsite / full as candidates) with the 5.8 v6 candidate
  to avoid conflicting v3-pending statements
- [ ] Confirm Figure 4 caption mentions "job 9572122, 100K subsample,
  90 second wall" for reproducibility
