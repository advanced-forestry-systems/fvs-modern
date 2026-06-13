# Oscar García's state-space modeling: review and refinements for the CONUS stand-level work

**Date:** 2026-06-11
**Purpose:** review García's age-independent state-space and GADA methodology and translate it
into concrete refinements for the CONUS stand-level models (top height, density, basal area)
and how they couple to the tree-level system.

## 1. García's framework, in his own terms

García builds whole-stand growth models as **dynamical systems**. The core ideas, drawn from
his state-space paper, his Stand Growth Models tutorial, the 1983 Biometrics height-growth
SDE, his mortality work, and the `resde` estimator:

- **Model the rate of change of a minimal state, not the input-output trajectory.** Describe
  the stand by a small sufficient state vector and specify how it evolves, dy/dt = f(y), rather
  than fitting yield as a function of age. His canonical whole-stand state is **top (dominant)
  height H, stems per hectare N, and a size variable (basal area or volume)**. Everything else
  (merchantable volume, the diameter distribution) is recovered from the state. This is exactly
  the three states we prototyped.
- **Bertalanffy-Richards as the deterministic core, written as an autonomous ODE.** Because the
  equation is autonomous (no explicit age term), projection is **age independent, annualized,
  and path invariant** by construction. The transition functions follow the generalized
  algebraic difference approach (GADA). We already reproduce the path-invariance property to
  machine precision, so the foundation is right.
- **Local and global parameters (the GADA free parameter).** Most parameters are common to all
  stands (global); one or a few are stand specific (local) and carry site or initial-condition
  effects. Site quality is therefore not an external covariate but an emergent local parameter.
  Top height acts as the biological clock in place of age.
- **Reducible stochastic differential equations with measurement error.** García does not fit
  the deterministic curve by least squares on the observed values. He writes growth as an SDE
  reducible by a change of variables to a linear SDE, carrying **both process noise
  (environmental variability) and observation error**, and estimates it by maximum likelihood
  on possibly noisy, unevenly spaced observations. His `resde` R package does this for
  univariate reducible SDEs with global and local parameters, by fixed or mixed effects.
- **Density and self-thinning modeled relative to size, not age.** García argues mortality and
  density should be modeled as change relative to change in size (the trajectory in N versus
  size space bending onto the Reineke self-thinning line), formulated as a rate so it handles
  thinning and disturbance, rather than as a fixed age-based trajectory.
- **Dimensional analysis and dimensionality reduction.** He uses dimensional consistency to
  cut free parameters and keep the state minimal.

## 2. What this says about our current stand-level models

Our three prototypes are structurally correct (the right three states, autonomous forms,
path invariant), but three of García's principles point at specific upgrades, two of which we
have already tested today.

### 2.1 Couple the states into one system (the biggest structural refinement)

We fit H, N, and BA as three independent equations. García's model is a **coupled** system:
the states drive each other. Today's results prove the point. Basal area carrying capacity is
not a function of a height site index at all; when we drove the BA maximum with CSPI it came
out the wrong sign and pinned at its bound. Refit with the maximum driven by **SDImax** (a
size-density, stand-state quantity), the sign is correct and sensible:

```
Gmax = 0.036 * SDImax ,  BA2 = BA1 + (Gmax - BA1)(1 - exp(-k*dt)) ,  k = 0.053/yr
```

So BA dynamics are governed by the stand's own carrying capacity (SDImax) and should take the
height increment and density as rate drivers, exactly García's coupling. The refactor is to
write the three as one state-space system with shared parameters rather than three regressions.

### 2.2 Estimate with the SDE likelihood that separates measurement error (the fix for top height)

Our top-height increment skill was negative (R2 -0.11), and the reason is now clear in García's
terms: we treated a noisy plot-level top-height estimate as exact and fit it by least squares,
so measurement error swamped the small real increment. García's SDE explicitly models
observation error alongside process noise. This is the principled fix, not more covariates.
`resde` is now installed on Cardinal; the next top-height fit should use it, with the site or
asymptote as the **local parameter** (mixed effects) rather than fixed ecoregion offsets.

### 2.3 Model density relative to size, and decompose into ingrowth plus self-thinning

Our density self-thinning works well and improved with a García-style power form on relative
density:

```
rate = 0.45 * RD^2.73 per year  (RD = SDI/SDImax) ,  R2 on ln(N) = 0.60
```

This is the **self-thinning (mortality) limb**: gentle below the boundary (about 1.7%/yr at
RD 0.3), accelerating sharply onto the Reineke line (45%/yr at RD 1). Per García, net density
is then ingrowth minus self-thinning, so this couples directly to the negbinom ingrowth count:
ingrowth adds stems from below, this term removes them near the boundary. That decomposition
is the right way to reconcile the recruitment model with the stand-level density equation.

## 3. Concrete refinements, prioritized

1. **Refactor the stand layer into one coupled state-space system** for (H, N, BA), states
   driving each other, with the local or site parameter carried by H and the carrying capacity
   by SDImax. Stop fitting three independent regressions.
2. **Adopt the reducible-SDE estimator (`resde`) with measurement error** for the height (and
   ultimately all) state equations. This is the direct fix for the top-height increment problem
   and is García's own method and software.
3. **Make site a local (mixed-effects) parameter** predicted from CSPI and traits, rather than
   a fixed covariate. This ties the stand layer to the site surface and the trait architecture
   we already built, and is faithful to the GADA local-parameter idea.
4. **Density = ingrowth (negbinom) minus self-thinning (RD power form)**, fit and reconciled
   together rather than as separate stories.
5. **Drive basal area by SDImax carrying capacity, with H increment and density in the rate**,
   confirmed by today's refit (g1 positive).
6. **Keep the state minimal and dimensionally consistent**; recover volume and the diameter
   distribution from the state, as García does, which also gives a second route (his diameter
   distribution recovery) alongside our Ritchie-Hann disaggregation to push stand state down to
   the tree list.
7. **Preserve base-age and path invariance** (already satisfied) and keep the autonomous
   Bertalanffy-Richards forms.

## 4. What we changed today, and what is next

Done today: confirmed the path-invariance property, fixed the basal-area driver (SDImax, not
CSPI), and moved density to the accelerating power-form self-thinning rate. Installed `resde`.

Next, on the same cadence: refit top height as a reducible SDE with measurement error and a
local site parameter (the principled increment fix), then assemble the three into a single
coupled state-space module that takes CSPI and SDImax as inputs and emits the (H, N, BA)
trajectory with credible bands, wired as the Phase A diagnostic against `36_conus_benchmark.R`
before any constrained disaggregation.

## Key references

García 1983, A stochastic differential equation model for the height growth of forest stands,
Biometrics 39:1059-1072. García 1994, The state-space approach in growth modelling, Canadian
Journal of Forest Research 24:1894-1903. García, Stand Growth Models: Theory and Practice
(tutorial). García, Forest Stands as Dynamical Systems: An Introduction. García 2009, A simple
and effective forest stand mortality model. García 2019, Estimating reducible stochastic
differential equations by conversion to a least-squares problem, Computational Statistics 34;
and the `resde` R package (CRAN).
