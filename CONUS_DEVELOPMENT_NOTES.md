# FVS-CONUS Development Notes

Saved from the April 9, 2026 planning session. These items need to be addressed
when active development begins on the `conus-variant` branch.

## Computational Strategy for National Fit

Fitting an ORGANON-family DG model with Cao/Weiskittel annualization loop to
~3.5M trees across CONUS with nested species/ecoregion RE. Plan:

1. **Greg's MLE fits as baseline**: Use nlsLM fixed-effects estimates as
   informative priors for the brms fit. This anchors the Bayesian model.
2. **ADVI initialization** with simple annualization (Option B) to get rough
   posteriors quickly (hours, not days)
3. **Targeted HMC on 10% stratified subsample** with full Stan annualization
   loop (Option A) to refine posteriors
4. **Evaluate** whether full HMC on the complete dataset is necessary or if
   the subsample posteriors are sufficiently stable
5. Per-variant Wykoff fits took 12 to 36 hours each on Cardinal; the CONUS
   ORGANON fit with the projection loop will be substantially more expensive
   per iteration. Budget 48 to 72 hours for the full fit on Cardinal.

## Random Effects Structure (DECIDED: nested, species > ecoregion)

Ecoregion nested within species: `(1 | species) + (1 | species:ecoregion)`
in brms/lme4 syntax. Rationale:

- A species' shade tolerance and growth strategy constrain how much the
  ecoregion effect can shift its growth rate. Nesting captures this.
- Reduces parameter count vs. fully crossed (no species x ecoregion cross
  for every combination)
- Species-level random slopes on competition terms (BAL effect) where
  sample size supports it
- Partial pooling handles data-sparse species x ecoregion cells naturally

## Annualization Strategy (DECIDED: Cao/Weiskittel method, all intervals)

### Background

Cao (2000) proposed fitting annual growth equations to periodic data by
iteratively projecting within the fitting objective. Weiskittel et al. (2007)
extended this to mixed-effects models via nlme. The key idea: instead of
naively dividing periodic increment by period length, you fit an annual growth
function and project it forward step by step, matching the TOTAL periodic
growth as the fitting target.

### The gr.hat() approach (from AnnualDiaGrowth.r)

The ORGANON-family DG equation predicts annual increment:

```
dg_annual = exp(b0 + b1*log(dbh+6) + b2*dbh + b3*log((cr+0.2)/1.2)
              + b4*log(si-1.37) + b5*bal/log(dbh+1) + b6*sqrt(ba))
```

For a tree observed over a period of `yip` years, gr.hat() projects forward
by iterating the annual equation year by year, linearly interpolating the
stand-level covariates (BAL, CR, BA) between start and end values:

```r
for (i in 1:yip) {
  gr <- annual_dg(dbh_current, bal_current, cr_current, si, ba_current, ...)
  dbh_current <- dbh_current + gr
  bal_current <- bal_current + (bal2 - bal1) / yip
  cr_current  <- cr_current  + (cr2  - cr1)  / yip
  ba_current  <- ba_current  + (ba2  - ba1)  / yip
}
predicted_periodic_growth <- dbh_current - dbh1
```

The loss function compares predicted_periodic_growth to observed diam_gr.

### Why this matters

Naive annualization (dg / yip) biases short-interval data low and long-interval
data high because growth compounds nonlinearly. Trees that grow faster early
in the period gain competitive advantage (larger dbh reduces BAL effect),
accelerating later growth. Naive annualization misses this feedback.

### Adaptation for brms

The challenge is that brms expects a formula-based model, not an iterative
projection loop. Three approaches for Bayesian fitting:

**Option A: Custom Stan function (recommended)**

Write the gr.hat() loop as a Stan user-defined function, then call it in
the brms nonlinear formula. This preserves the exact Cao/Weiskittel approach
within the Bayesian framework.

```r
# brms nonlinear formula (conceptual)
bf(diam_gr ~ gr_hat_stan(d1, bal1, bal2, cr1, cr2, ba1, ba2, si, yip,
                          b0, b1, b2, b3, b4, b5, b6),
   b0 ~ 1 + (1 | species) + (1 | species:ecoregion),
   b1 + b2 + b3 + b4 + b5 + b6 ~ 1,
   nl = TRUE)
```

The Stan function would implement the year-by-year projection identically
to the R gr.hat() function. This is the most faithful adaptation.

**Option B: Approximate with simple annualization for initial fit**

Use naive annualization (dg / yip) as the response variable for an initial
ADVI fit to get starting values, then refit with the full Stan loop for
final posteriors. Practical for getting the pipeline running quickly.

**Option C: Pre-annualize with Cao method in R, fit annual model**

Use the nlme or nlsLM fits (as Greg is doing) to get MLE annualized
predictions, then use those as the target for a brms annual model without
the projection loop. Less principled but computationally simpler.

Recommendation: Start with Option B to get infrastructure running, then
implement Option A for the final published model.

### Data decisions

- Use ALL available remeasurement combinations (5-year, 10-year, etc.)
- Do NOT filter to a single interval length
- The Cao/Weiskittel method handles variable intervals naturally via the
  projection loop (yip parameter)
- Include a plot-level random effect to account for correlation when the
  same plot has overlapping intervals (e.g., 0-5, 0-10 from same baseline)

### References

- Cao, Q.V. 2000. Prediction of annual diameter growth and survival for
  individual trees from periodic measurements. Forest Science 46:127-131.
- Weiskittel, A.R.; Garber, S.M.; Johnson, G.P.; Maguire, D.A.;
  Monserud, R.A. 2007. Annualized diameter and height growth equations for
  Pacific Northwest plantation-grown Douglas-fir, western hemlock, and red
  alder. Forest Ecology and Management 250:266-278.
- Hann, D.W.; Marshall, D.M.; Hanus, M.L. 2006. Reanalysis of the
  SMC-ORGANON equations for diameter-growth rate, height-growth rate, and
  mortality rate of Douglas-fir. Research Contribution 49. Oregon State.

## Greg Johnson's fvs_remodeling repo (github.com/gregjohnsonbiometrics/fvs_remodeling)

As of April 2026, Greg's repo uses:

- **Equation form**: ORGANON-family (Hann et al. 2006), fitted via
  minpack.lm::nlsLM (nonlinear least squares, fixed effects only)
- **Data**: 3.48M FIA observations, 84 species with >= 5,000 observations
- **Covariates**: DBH, CR, height, BAL, elevation, EMT
- **Ecoregion**: EPA L3 ecoregions joined spatially to plot coordinates;
  used for residual diagnostics but NOT yet as model covariates or RE
- **Known gap**: Spatial autocorrelation in DG residuals persists even
  with elevation + EMT. Height growth parameters sometimes have
  biologically implausible signs. CR coefficients often nonsignificant.
- **No Bayesian fitting yet**: All models are fixed-effects MLE

The CONUS brms pipeline needs to extend this by:

1. Adding species and ecoregion nested RE
2. Implementing the Cao/Weiskittel annualization within Stan
3. Adding climate covariates beyond elevation + EMT
4. Evaluating ClimateSI as a productivity proxy
5. Proper uncertainty quantification via posterior distributions

## FIA Data Pipeline for CONUS

- Script 01 currently fetches data per variant. CONUS needs a national fetch.
- **Deduplication step** keyed on CN/PLT_CN is essential: many plots appear in
  multiple variants' geographic boundaries.
- `plot_raster_lookup.csv` (1.9M records) needs to cover ALL CONUS FIA plots,
  not just those in current variant boundaries. Verify ClimateNA extraction
  covers the full plot network.

## Crown Ratio Redesign

Crown ratio modeling is the weakest component (R-squared 0.13 to 0.38). For
CONUS, explore:

- Nonlinear alternatives to the current linear model
- Allometric crown models derived from airborne lidar (GEDI, 3DEP)
- Empirical lookup tables by species group and competition class
- Crown competition factor (CCF) as an alternative state variable

## Ingrowth Model

Three options for CONUS:

a) **Empirical FIA-based**: Observed ingrowth rates by species group and stand
   condition. Simple, defensible, not climate-sensitive. Currently implemented.
b) **Establishment model**: Climate niche models (SDMs) combined with gap/
   disturbance detection.
c) **Hybrid**: Empirical rates modulated by climate suitability scores.

Start with (a) for initial release, develop (c) as a research priority.

## Species Code Unification

Traditional FVS uses sequential integer species codes per variant (species 1,
2, 3...) that differ between variants. CONUS uses FIA SPCD directly. The NSBE
equations already work this way. Fortran integration needs to handle:

- Dynamic species array allocation (no fixed MAXSP COMMON blocks)
- SPCD-based coefficient lookup instead of sequential index
- Graceful handling of species not in the training data (Jenkins fallback)

## Soil Data Integration

gSSURGO data is not yet integrated. For CONUS, soil properties (available water
capacity, depth to restrictive layer, clay content) could improve site
productivity and drought response predictions. Consider:

- Pre-extracting soil properties for all FIA plot locations
- Adding soil covariates to the DG and mortality models
- Evaluating marginal improvement over climate-only models

## Validation Framework

1. **Within-sample**: Cross-validated RMSE, R-squared, bias by ecoregion
2. **Out-of-sample**: Hold-out 20% of FIA plots for independent validation
3. **Variant boundary test**: Compare projections for plots near variant
   borders processed by (a) assigned variant, (b) adjacent variant,
   (c) FVS-CONUS
4. **Long-term realism**: 50-year projections checked against ecological
   bounds (Reineke SDI limits, expected mortality patterns)
5. **Species range edges**: Test performance at range margins where
   variant-specific fits have limited data

## Long-Term Model Comparisons (Bakuzis Matrix Update)

The v3 ablation benchmarks only cover short-term projections (~6.5 year
average FIA remeasurement). ClimateSI and SDImax contribute <0.1% at this
horizon. However, their effects compound nonlinearly over longer periods:

- SDImax onset at RDI=0.70 only triggers in dense stands and only after
  enough growth pushes RDI above the threshold. Over 50 to 100 years,
  cumulative mortality differences will be substantial.
- ClimateSI drives site-dependent growth rate differences that compound:
  a 5% annual DG difference becomes a 40%+ BA difference at 50 years.

Plan: Run 100-year projections for a representative set of initial
conditions (Bakuzis matrix: species x site x density combinations) under
both default FVS and calibrated (v3) to produce updated long-term
comparison figures.

Implementation:
- Use script 17 (stand_projection_engine.R) for long-term projections
- 4 species groups x 3 site classes x 3 initial densities = 36 scenarios
- Project at 5-year steps for 100 years under both default and calibrated
- Compare: BA trajectories, mortality onset timing, volume at rotation age
- Produce Bakuzis-style matrix figure (patchwork of trajectory plots)
