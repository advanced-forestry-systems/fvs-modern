# CONUS-Wide FVS Component Model Specification

**Authors**: Aaron Weiskittel, Greg Johnson
**Date**: 2026-04-11 (updated with K1/K2 estimation, RD, metric, dynamic CR)
**Status**: Draft specification for review

## 1. Overview

This document specifies four tree-level component equations for a unified CONUS-wide
FVS variant. Model forms are adopted from ORGANON (Hann et al. 2002, 2006; Hanus
et al. 2000) and extended with Bayesian hierarchical structure (species x ecodivision
random effects), climate covariates, relative density interactions, and annualized
likelihoods for periodic FIA remeasurement data.

**ALL UNITS METRIC**: DBH in cm, HT in m, BA in m2/ha, BAL in m2/ha.
FIA imperial data is converted at the data preparation stage.

**Site productivity**: ClimateSI (m), BGI (biomass growth index from asymptotic AGB
analysis), or climate PCA scores. Traditional FIA SICOND is NOT used.

The four components are:

1. **Diameter growth** (DG): annual inside-bark diameter increment (cm/yr)
2. **Height growth** (HG): annual height increment via potential x modifier (m/yr)
3. **Height to crown base** (HCB): static prediction of crown base position (m)
4. **Survival** (PS): annual survival probability (gompit with exposure)

Crown ratio is derived dynamically from HCB: **CR = 1 - HCB/HT**

All models are fit simultaneously to national FIA remeasurement data using CmdStan
via cmdstanr in R.

## 2. Data Sources

**Primary fitting data**: FIA periodic and annual remeasurement pairs from all 48
CONUS states. Trees >= 1.0 inch DBH with valid remeasurement intervals.

**Required columns per tree-period observation**:
- Tree identifiers: TRE_CN, PLT_CN, STATECD, COUNTYCD
- Species: SPCD (FIA species code, ~300 species CONUS-wide)
- Ecoregion: ECOSUBCD or ecodivision (Bailey's system, ~30 divisions)
- Size at t1: DBH1, HT1, CR1, HCB1 (derived as HT1*(1-CR1))
- Size at t2: DBH2, HT2, CR2, HCB2
- Alive/dead at t2: STATUS2 (1=alive, 2=dead)
- Competition at t1: BAL1, BA1, CCFL1, TPA1
- Competition at t2: BAL2, BA2, CCFL2, TPA2 (for annualization interpolation)
- Site: SI (SICOND from FIA, variant-specific), ClimateSI, climate PCA scores
- Topography: SLOPE, ASPECT, ELEV
- Climate: EMT (extreme minimum temperature), TD (temperature difference/continentality)
- Measurement interval: YEARS (t2 - t1)

**Species threshold**: >= 5,000 tree-period observations for species-level estimation.
Species below threshold are pooled into genus or hardwood/softwood functional groups
with informative priors from the group mean.

## 3. Hierarchical Structure

All four components share the same hierarchical framework:

```
ecodivision (d = 1, ..., D)   ~30 ecodivisions
  species (s = 1, ..., S)     ~150-300 species (depending on threshold)
    plot (p = 1, ..., P)       ~100,000+ FIA plots
      tree (i = 1, ..., N)     ~16M tree observations
```

**Random effects specification**:
- Species random intercepts on all components (non-centered parameterization)
- Ecodivision random intercepts on diameter growth and mortality
- Species x ecodivision interaction random intercept on diameter growth (where data support it)
- Plot-level random effects NOT included in Stan (too many groups for HMC); handled via
  calibration at prediction time using BLUP-style adjustment

The non-centered parameterization for species intercepts:
```
b0[s] = mu_b0 + sigma_b0 * z_b0[s]
z_b0[s] ~ Normal(0, 1)
```

And the crossed ecodivision effect:
```
b0_total[s,d] = mu_b0 + sigma_b0_sp * z_sp[s] + sigma_b0_eco * z_eco[d]
```

## 4. Component 1: Diameter Growth

### 4.1 Model Form (ORGANON SWO; Hann et al. 2002, 2006)

Annual diameter increment (all metric: DBH in cm, BA in m2/ha, BAL in m2/ha):
```
ln(DG) = b0[s,d] + b1*ln(DBH + K1) + b2*DBH^K2 + b3*ln((CR + 0.2) / 1.2)
         + b4*ln(SITE_PROD) + b5*(BAL / ln(DBH + K4)) + b6*sqrt(BA)
         + b7*CLIM1 + b8*CLIM2 + b9*RD + b10*RD*ln(BAL + 5.0)
```

where:
- b0[s,d] = species x ecodivision random intercept
- b1 through b10 = shared (fixed) slopes across all species
- **K1 = ESTIMATED** (prior: Normal(6.0, 2.0), lower bound 0.1; additive DBH offset)
- **K2 = ESTIMATED** (prior: Normal(1.0, 0.5), lower bound 0.01; power on DBH)
- K3 = 1.0 (fixed; power on BAL; all SWO species use 1.0)
- K4 = 2.7 (fixed; consistent across ORGANON variants)
- SITE_PROD = ClimateSI (m), BGI, or climate PCA score (NO traditional SI)
- CLIM1, CLIM2 = climate PCA scores (top 2 PCs from ClimatePCA.r)
- RD = Curtis relative density (BA / sqrt(QMD), metric)
- ln(BAL + 5.0) = log transformed BAL for interaction stability
- DG = exp(ln(DG)) with low CR adjustment: CRADJ = 1 if CR > 0.17, else
  1 - exp(-(25*CR)^2)

**Rationale for estimating K1 and K2**: The ORGANON SWO variants hardcode K1
at 5.0 or 6.0 and K2 at 1.0 depending on species group and reanalysis vintage.
For a CONUS wide model spanning hundreds of species, estimating these shape
constants from the data provides flexibility for the ln(DBH + K1) and DBH^K2
terms to adapt to the aggregate growth response. Priors centered on ORGANON
defaults keep estimates anchored to prior biological knowledge.

**Rationale for RD covariates**: Relative density captures stocking relative
to carrying capacity, complementing the absolute density measures (BA, BAL).
The RD x ln(BAL) interaction allows the competitive effect of BAL to vary
with overall stand density. This follows the mortality modeling example that
uses sqrt(BA*RD) as a density metric.

**Comparison with Greg's form**: Greg uses a compound term
`log((DBH+1)^2 / (CR*HT+1)^b3)` that merges size and crown into one predictor.
We retain the separated ORGANON form because: (a) it preserves identifiability of
the CR and DBH effects for diagnostic purposes, (b) the log((CR+0.2)/1.2)
transformation is bounded and well-behaved at CR extremes, and (c) it allows site
index and competition to enter independently. Greg's climate terms (EMT, ELEV) are
incorporated as additional covariates rather than replacing SI.

### 4.2 Annualization (Cao 2000; Weiskittel et al. 2007)

For a tree with measurement interval of T years and time-varying covariates
(BAL, CR, BA), annual predictions are accumulated iteratively:

```
Initialize: d_curr = DBH1, bal_curr = BAL1, cr_curr = CR1, ba_curr = BA1
Linear interpolation rates:
  d_bal = (BAL2 - BAL1) / T
  d_cr  = (CR2 - CR1) / T
  d_ba  = (BA2 - BA1) / T

For t = 1 to T:
  dg_t = exp(b0[s,d] + b1*ln(d_curr + K1) + ... + b6*sqrt(ba_curr) + ...)
  dg_t = dg_t * CRADJ(cr_curr)
  d_curr = d_curr + dg_t
  bal_curr = bal_curr + d_bal
  cr_curr  = cr_curr + d_cr
  ba_curr  = ba_curr + d_ba

Predicted periodic growth: DG_hat = d_curr - DBH1
Observed periodic growth:  DG_obs = DBH2 - DBH1
```

The likelihood is:
```
DG_obs ~ Normal(DG_hat, sigma)
```

or equivalently with weighting by predicted value:
```
(DG_obs - DG_hat) / sqrt(DG_hat) ~ Normal(0, sigma)
```

### 4.3 Priors

| Parameter | Prior | Justification |
|-----------|-------|---------------|
| mu_b0 | Normal(-5.0, 2.0) | Centered on ORGANON SWO DF intercept |
| sigma_b0_sp | Exponential(1) | Species variation in base growth rate |
| sigma_b0_eco | Exponential(2) | Ecodivision variation (tighter than species) |
| K1 | Normal(6.0, 2.0), lb=0.1 | Additive DBH offset; anchored to Hann et al. 2006 |
| K2 | Normal(1.0, 0.5), lb=0.01 | Power on DBH; most SWO species use 1.0 |
| b1 | Normal(0.8, 0.5) | Positive effect of size, ~0.84 in SWO |
| b2 | Normal(-0.04, 0.03) | Negative DBH^K2 effect |
| b3 | Normal(1.0, 0.5) | Positive CR effect, ~1.0 to 1.16 in SWO |
| b4 | Normal(0.8, 0.5) | Positive site productivity effect |
| b5 | Normal(-0.008, 0.005) | Negative BAL effect |
| b6 | Normal(-0.04, 0.03) | Negative stand density effect |
| b7, b8 | Normal(0, 1) | Climate effects (weakly informative) |
| b9 | Normal(-0.5, 0.5) | Negative RD effect (higher RD = less growth) |
| b10 | Normal(0, 0.3) | RD x ln(BAL) interaction (weakly informative) |
| sigma | Exponential(1) | Observation error |

## 5. Component 2: Height Growth

### 5.1 Model Form (ORGANON; Hann et al. 2002, 2006)

Height growth uses a potential x modifier approach:

**Potential height growth**: Derived from a Chapman-Richards site index curve
```
HG_pot = f(SI, HT, age)
```
where f() is species-specific from published site index equations. For a
CONUS-wide model without reliable age, we use an alternative formulation
following Greg's approach:

```
HG_pot = HT_max * b1_hg * b2_hg * (1 - exp(-b1_hg * HT))^(b2_hg - 1) * exp(-b1_hg * HT)
```

This is the derivative of the Chapman-Richards curve with respect to time,
expressed in terms of current height rather than age. HT_max is the asymptotic
height, treated as a function of site productivity.

**Modifier** (ORGANON HG_SWO form):
```
FCR = -P5 * (1 - CR)^P6 * exp(P7 * CCFL^0.5)
B0_hg = P1 * exp(P2 * CCFL)
B1_hg = exp(P3 * CCFL^P4)
MODIFIER = P8 * (B0_hg + (B1_hg - B0_hg) * exp(FCR))
```

**Annual height growth**:
```
HG = HG_pot * MODIFIER * CRADJ + b_clim1*CLIM1 + b_clim2*CLIM2
```

where CRADJ = 1 if CR > 0.17, else 1 - exp(-(25*CR)^2)

**Simplified alternative** (for initial implementation): Following Greg's
direct approach:
```
ln(HG) = b0_hg[s,d] + b1_hg * ln(HT_max) + b2_hg * b3_hg * CR^b4_hg
         * exp(-b2_hg*HT - b5_hg*CCFL - b6_hg*CCH^0.5
         + b7_hg*ELEV + b8_hg*TD^0.5 + b9_hg*EMT)
         + (b3_hg - 1) * ln(1 - exp(-b2_hg * HT))
```

### 5.2 Annualization

Identical iterative approach as diameter growth. For each year t in the
measurement period, predict annual HG, update HT_curr, and interpolate
time-varying covariates (CCFL, CCH, CR).

### 5.3 Hierarchical Structure

- b0_hg[s,d]: species x ecodivision random intercept
- Shape parameters (b2_hg through b4_hg): shared across species initially;
  species random slopes can be added if data support it
- HT_max: can be treated as a function of ClimateSI with species random effect

## 6. Component 3: Height to Crown Base

### 6.1 Model Form (ORGANON HCB_SWO; Hanus et al. 2000)

Static logistic prediction of HCB as a proportion of total height.
All metric: HT in m, DBH in cm, BA in m2/ha:
```
HCB = HT / (1 + exp(b0_hcb[s,d] + b1_hcb*HT + b2_hcb*CCFL + b3_hcb*ln(BA)
             + b4_hcb*(DBH/HT) + b5_hcb*SITE_PROD + b6_hcb*RD))
```

where:
- b0_hcb[s,d] = species x ecodivision random intercept
- CCFL = crown competition factor in larger trees
- BA = stand basal area (m2/ha)
- DBH/HT = form ratio (cm/m; captures taper differences)
- SITE_PROD = ClimateSI (m), BGI, or climate PCA (NO traditional SI)
- RD = Curtis relative density (BA/sqrt(QMD), metric)

This predicts HCB directly rather than modeling crown ratio change over time.
Crown ratio is then derived dynamically as **CR = 1 - HCB/HT**. This avoids
the R2 = 0.237 problem with the current direct CR change model.

### 6.2 Dynamic Crown Ratio Model

Crown ratio is NOT modeled directly. Instead, CR is derived from the static HCB
model evaluated at current tree and stand attributes:

```
CR(t) = 1 - HCB_pred(t) / HT(t)
```

For projection, crown recession (upward movement of the crown base) is computed
as the difference between predicted HCB at t2 and t1:

```
HCBG = max(HCB_pred(HT2, DBH2, CCFL2, BA2, RD2, ...) -
           HCB_pred(HT1, DBH1, CCFL1, BA1, RD1, ...), 0)
HCB2 = HCB1 + HCBG
CR2  = 1 - HCB2 / HT2
CR2  = clamp(CR2, 0.01, 1.0)
```

With the constraint that HCBG >= 0 (crown base does not move downward). This is
the approach used in ORGANON (crngrow.f90, lines 116-119). The dynamic CR feeds
back into the DG and HG equations at each annual time step during projection,
creating the key feedback loop: competition -> HCB -> CR -> growth -> competition.

The `predict_dynamic_cr()` function in 33_fit_hcb_organon.R implements this logic.

### 6.3 Fitting

The HCB model is fit as a static cross-sectional equation to all FIA trees with
measured crown ratio at time t1 (converting to HCB = HT * (1 - CR)). This avoids
requiring paired measurements and uses the full dataset.

Likelihood:
```
HCB_obs ~ Normal(HCB_pred, sigma_hcb)
```

## 7. Component 4: Survival (Mortality)

### 7.1 Model Form (ORGANON PM_SWO; Hann and Hanus 2001)

Complementary log-log (gompit) link with exposure offset for annualization.
All metric: DBH in cm, BA in m2/ha, BAL in m2/ha, site productivity in m:

```
cloglog(P_mort_annual) = b0_m[s,d] + b1_m*DBH + b2_m*DBH^2 + b3_m*CR
                         + b4_m*SITE_PROD + b5_m*ln(BAL + 5.0) + b6_m*CLIM1
                         + b7_m*sqrt(BA*RD)
```

Five-year (or T-year) survival probability:
```
P_surv_T = (1 - P_mort_annual)^T = exp(-exp(eta) * T)
```

where eta is the linear predictor above. Using the gompit link with
`offset = log(T)`:
```
cloglog(P_mort_T) = eta + log(T)
```

This follows from the uploaded Mortality Modeling Example.R (Gompit.e.m1),
which demonstrated the gompit with exposure outperforming alternatives.

### 7.2 Covariates

Following ORGANON SWO, adapted to metric units:
- DBH (cm): linear and quadratic (captures U shaped mortality; priors rescaled for cm)
- CR: crown ratio as vigor indicator (strong negative effect on mortality)
- SITE_PROD: ClimateSI (m), BGI, or climate PCA (NO traditional SI; no +4.5 offset)
- ln(BAL + 5.0): log transformed BAL in m2/ha (ORGANON PM_SWO transform;
  stabilizes competition effect at high BAL values)
- CLIM1: climate variable (EMT or drought index) for regional mortality patterns
- sqrt(BA*RD): relative density interaction from the Weiskittel mortality example;
  captures the joint effect of absolute density and stocking relative to maximum

The ORGANON SWO form also includes BAL*exp(b7*OG) for old growth interaction,
but OG is rarely available in FIA. We omit this term initially.

### 7.3 Hierarchical Structure

- b0_m[s,d]: species x ecodivision random intercept (non-centered)
- Ecodivision effects are especially important for mortality because background
  mortality rates vary substantially by region (fire, insects, drought)
- Species tolerance class can serve as an additional grouping level:
  `b0_m[s,d] = mu_tol[tolerance_class] + sigma_sp * z_sp[s] + sigma_eco * z_eco[d]`

### 7.4 Priors

| Parameter | Prior | Justification |
|-----------|-------|---------------|
| mu_b0_m | Normal(-3.0, 2.0) | Low base annual mortality (~5%) |
| sigma_sp_m | Exponential(1) | Species variation |
| sigma_eco_m | Exponential(2) | Ecodivision variation |
| b1_m | Normal(-0.06, 0.04) | Negative DBH in cm (rescaled from inches) |
| b2_m | Normal(0.0005, 0.0003) | Positive DBH^2 in cm^2 (U shape for large trees) |
| b3_m | Normal(-3.0, 1.0) | Strong negative CR effect |
| b4_m | Normal(0.01, 0.01) | Weak site productivity effect |
| b5_m | Normal(0.005, 0.003) | Positive ln(BAL+5) effect (competition) |
| b6_m | Normal(0, 0.5) | Climate effect (weakly informative) |
| b7_m | Normal(0.1, 0.1) | Positive sqrt(BA*RD) effect (density interaction) |

## 8. Site Productivity Strategy

Site productivity enters all four components. Traditional FIA SICOND is NOT used
because variant-specific definitions are inconsistent across CONUS and many plots
lack reliable site index. We fit parallel models with three alternatives:

1. **ClimateSI** (m): Random forest model predicting site index from ClimateNA
   normals (30+ derived variables). Trained on FIA plots with reliable measured SI.
   Provides a continuous, climate-based productivity surface for all CONUS plots.
   Model: `rf_NASI_ClimateNA_ECO_model.rds` on Cardinal.

2. **BGI** (Biomass Growth Index): Derived from the asymptotic AGB analysis on
   Cardinal (asym_agb_analysis pipeline). Uses Chapman-Richards asymptote with
   plot-level BLUPs. Captures realized productivity from observed biomass
   accumulation rather than height growth potential.

3. **Climate PCA**: First 4 principal components of ~30+ derived ClimateNA variables
   (MAT, MAP, EMT, EXT, MWMT, MCMT, MSP, FFP, DD5, DD_0, plus ratios and
   interactions). Rasters: `CONUS_Climate_Embedding_Top4.tif` on Cardinal.
   PCA model: `CONUS_Climate_PCA_Model.rds`.

All three versions use identical model structure; only the SITE_PROD variable changes.
Final selection based on cross-validated prediction accuracy at the stand level,
biological realism of parameter estimates, and spatial consistency of residuals.

**Relative density (RD)** serves as a complementary measure of site utilization.
Curtis RD (BA/sqrt(QMD)) is used across all components. Where SDImax estimates
are available from the brms quantile regression (brms.SDImax.1-27-24.csv), plot
specific RD = SDI/SDImax can be used for a more ecologically meaningful metric.

## 9. Implementation Plan

### Phase 1: Data Assembly (R script: 30_build_conus_dataset.R)
- Pull all FIA remeasurement pairs CONUS-wide
- Compute derived variables (BAL, CCFL, CR, HCB)
- Attach ecodivision codes via spatial join
- Attach climate normals (PRISM or ClimateNA)
- Export species observation counts; apply 5,000-obs threshold

### Phase 2: Model Fitting (Stan + R)
- 31_fit_dg_organon.R + organon_dg_conus.stan
- 32_fit_hg_organon.R + organon_hg_conus.stan
- 33_fit_hcb_organon.R + organon_hcb_conus.stan
- 34_fit_mortality_gompit.R + gompit_mortality_conus.stan
- Each script handles data preparation, Stan compilation, HMC sampling,
  convergence diagnostics, and posterior export

### Phase 3: Comparison and Validation
- 35_compare_site_productivity.R: head-to-head comparison of SI vs ClimateSI vs PCA
- 36_conus_benchmark.R: stand-level FIA benchmark engine adapted for single CONUS model
- 37_bakuzis_matrix_conus.R: biological realism checks across ecodivisions

### Phase 4: Integration
- 38_posterior_to_json_conus.R: export calibrated parameters as JSON config
- Update config_loader.py for CONUS variant
- Integration with projection engine (script 17) and benchmark engine (script 19)

## 10. References

- Cao, Q.V. 2000. Prediction of annual diameter growth and survival for individual
  trees from periodic measurements. For. Sci. 46:127-131.
- Hann, D.W.; Hanus, M.L. 2001. Enhanced mortality equations for trees in the mixed
  conifer zone of southwest Oregon. FRL Res. Contrib. 34.
- Hann, D.W.; Hanus, M.L. 2002. Enhanced diameter-growth-rate equations for
  undamaged and damaged trees in southwest Oregon. FRL Res. Contrib. 39.
- Hann, D.W.; Marshall, D.M.; Hanus, M.L. 2006. Reanalysis of the SMC-ORGANON
  equations for diameter-growth rate, height-growth rate, and mortality rate of
  Douglas-fir. FRL Res. Contrib. 49.
- Hanus, M.L.; Hann, D.W.; Marshall, D.D. 2000. Predicting height to crown base
  and crown recession in unthinned young Douglas-fir stands. FRL Res. Contrib. 29.
- Weiskittel, A.R.; Garber, S.M.; Johnson, G.P.; Maguire, D.A.; Monserud, R.A. 2007.
  Annualized diameter and height growth equations for Pacific Northwest plantation-grown
  Douglas-fir, western hemlock, and red alder. For. Ecol. Manage. 250:266-278.
- Weiskittel, A.R.; Hann, D.W.; Kershaw, J.A.; Vanclay, J.K. 2011. Forest Growth
  and Yield Modeling. Wiley. 415 p.
