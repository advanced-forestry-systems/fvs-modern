# Do García's principles apply to the tree-level equations? Yes, and annualization is the lens

**Date:** 2026-06-11
**Question (Aaron):** do García's state-space principles apply to our tree-level components, and
note that everything needs to be annualized.

Short answer: yes, several of his principles already operate at the tree level, one of them
(the global/local GADA decomposition) is essentially what our hierarchical Bayesian
architecture already does, and annualization, the principle you flagged, is mostly in place but
has two specific gaps and one deeper refinement. I audited the production Stan models directly.

## 1. Annualization audit of the six tree-level components

Annualization means the likelihood is written so a variable measurement interval is handled
correctly and projection is interval-free. Reading the actual Stan code:

| component | how the interval enters | annualized? |
|---|---|---|
| Diameter growth (Kuehne v8 BGI) | annual increment response, error scaled by sqrt(years) | yes |
| Height growth (unified v8rd) | annual increment `hg_obs_a`, likelihood sigma/sqrt(years) | yes |
| Survival (gompit) | cloglog on annual survival with exposure: P(survive T) = exp(-hazard*T) | yes, exemplary |
| Ingrowth (negbinom) | log(years) offset on the count | yes |
| Height to crown base | "Cross-sectional, no annualization" (static prediction) | no, static |
| Crown ratio (CR2-direct) | logit(CR2) from logit(CR1) and covariates, no interval term | no, gap |

So four of six are properly annualized, and survival is the most García-faithful of all: it is
an instantaneous hazard integrated over the interval, which is exactly a differential equation
solution. The two exceptions matter differently:

- **HCB is static, not a rate.** It predicts crown base from current attributes, and crown-base
  change is the difference of two static predictions. This is actually path invariant, because
  a pure function of the current state gives the same change regardless of the path taken, so
  it does not break annualization. But it is not a rate equation, and the constraint that crown
  base cannot move down is ad hoc. If you want full state-space consistency, HCB becomes a rate
  (or stays static but the constraint is replaced by a monotone formulation).
- **CR2-direct is the one genuine annualization gap.** It maps logit(CR2) from logit(CR1) plus
  covariates with no interval term, so a five-year and a fifteen-year crown-ratio change are
  treated the same. This biases the transition and should carry an exposure or interval term,
  the same fix used everywhere else. (If CR is instead derived from HCB and height as
  CR = 1 - HCB/HT, this is moot, but the stress test used the direct CR2, so the gap is real.)

## 2. The other García principles at the tree level

- **Global and local parameters (GADA): already done.** García splits parameters into global
  (common to all stands) and one or a few local (stand or site specific). Our tree-level
  architecture is exactly this: the trait-driven species-free term is the global structure, and
  the species and ecoregion random effects are the local parameters. The species random
  intercept is the GADA free parameter. So the hierarchical Bayesian calibration is already
  García-consistent on this axis, which is a satisfying convergence rather than a refinement.
- **Reducible SDE with measurement error: the strongest tree-level refinement.** Our height
  components (HG and HT-DBH) are built on tree height, which in FIA is frequently modeled or
  imputed and carries real measurement error. We fit annual height growth with a normal
  likelihood that treats the observed increment as exact, the same mistake that gave the stand
  top-height model negative increment skill. García's reducible SDE separates process noise from
  observation error. Applying it to HG and HT-DBH should sharpen them; DG is less affected
  because DBH is measured precisely. A demonstration fit of the top-height SDE with and without
  a measurement-error term is running now (job 11553482) to quantify how large that observation
  noise is.
- **Rate of change rather than difference.** DG, HG, survival, and ingrowth are rate-like;
  HCB and CR2 are differences. García prefers rate formulations because they extend cleanly to
  thinning and disturbance, which is also where our modifier layer operates, so reformulating
  HCB and CR as rates would compose more naturally with the modifiers.
- **Minimal, dimensionally consistent state.** The competition predictors (BA, BAL, SDI, RD,
  CCF, CCFL) are collinear. García's dimensional reduction argues for trimming to a minimal
  sufficient set, which would reduce parameter redundancy without losing signal.

## 3. Concrete tree-level refinements, prioritized

1. **Add an interval term to CR2-direct** so it is annualized like the rest. The one clear gap.
2. **Refit the height components (HG, HT-DBH) with a measurement-error model** in García's
   reducible-SDE spirit, since height carries observation error that the current exact-response
   likelihood ignores. The running top-height SDE quantifies the stakes.
3. **Test path invariance explicitly for DG and HG** (ten one-year steps equal one ten-year
   step), since the within-interval covariate change (the Cao annualization at prediction time)
   can break exact invariance even when the fit is annualized.
4. **Reformulate HCB (and optionally CR) as rate equations** if you want full state-space
   consistency, replacing the ad hoc crown-base monotonicity constraint with a monotone rate.
5. **Trim the collinear competition predictors** toward a minimal sufficient density state.
6. **Affirm the global/local structure**; the hierarchical species and ecoregion effects are
   already García's GADA decomposition, so no change there, just recognition.

## 3a. CR2 annualization fix, validated (job 11553509)

Implemented and tested the principled fix for the one annualization gap. The current CR2 model
maps logit(CR2) from logit(CR1) with no interval term, so it predicts a fixed change even over a
zero-length interval. The García-consistent replacement is an exponential approach to an
equilibrium crown ratio:

```
logit(CR2) = E + (logit(CR1) - E) * exp(-k * T_years)
E = b0 + trait + random effects + b1 dbh + b2 dbh^2 + b3 ba + b4 bal + b6 ln_csi   (equilibrium)
```

A quick validation fit on 200,000 trees (mean interval 6.5 years) confirms it works: the annual
approach rate k estimates to 0.064 per year (a crown-ratio adjustment half-life near 11 years,
biologically sensible), and the form is path invariant to machine precision (ten one-year steps
equal one ten-year step). At a zero interval it correctly returns CR2 = CR1, which the current
model cannot. Its in-sample R2 on logit(CR2) is 0.40 against 0.43 for the current form; the small
gap is the trivial autocorrelation the current b_cr1 term picks up, not real skill, and the
production test is a held-out ΔLOO on the Stan models. The annualized Stan model
(`crown_ratio_t2_annualized.stan`, with a `T_years` data input and `k` replacing `b_cr1`) is
written and ready; the remaining step is to add `T_years` to the CR2 fit driver and run the
current-versus-annualized ΔLOO.

## 4. The unifying point

Annualization is not just bookkeeping for FIA's variable intervals; it is the same property
(an autonomous rate integrated over any interval) that makes the stand-level state-space model
path invariant. The tree level already satisfies it for four of six components and embodies the
global/local GADA structure throughout. The two real gaps (CR2 annualization, height
measurement error) are the same two lessons the stand-level work surfaced, which is the sign
that one consistent García-style treatment, rates plus measurement error plus global/local
parameters, applies cleanly from the tree to the stand.
