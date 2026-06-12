# Ingrowth, species composition, and a stand level constraint layer: suggestions

**Date:** 2026-06-11
**Prompted by Aaron:** v8 status, the ingrowth and species composition work, and whether to
build simple age independent yet annualized stand level equations (Garcia tradition) for
basal area, stem density, and top height to constrain the tree level growth and mortality.

## 1. CSPI v8: yes it exists, and v7 should still win

v8 is built (v7 stack plus five Sentinel 2 bands plus four MOD17 GPP/NPP layers). It buys
almost nothing over v7: OOB R2 0.702 vs 0.700, plot blocked CV 0.688 vs 0.686, and v8 is
still carried by v6_distill, elevation, and soil while the spectral terms sit low in
importance. The v7 distillation screen (job 11505475) then showed the v7 surface already
explains more component residual variance than observed site index. Recommendation stands:
lock v7 as the deployable 30 m production surface, keep v6 (1 km) as the accuracy reference,
and retire the v8 remote sensing branch unless a specific application needs sub stand
spectral detail. The cost of the v8 data pipeline is not justified by 0.001 R2.

## 2. Ingrowth and species composition

What exists: a plot level negative binomial recruitment count model (`35_fit_ingrowth_negbinom`
through v4) with EPA L1 to L3 hierarchical random effects, where v4 added top height (HT40)
and rd_sdimax and found top height far more predictive than stand age; a hurdle variant
(`35b`); and a first species composition model (`36_fit_ingrowth_species_composition`). The
manuscript treats ingrowth as the stand level seventh component and the modifier table keeps
it on the common (not trait mediated) form.

Refinements worth making:

- **Make composition a compositional model, trait driven.** The count needs no traits, but
  *which* species recruit is silvically driven, so this is where the trait block earns its
  place for ingrowth. Predict species or group shares with a hierarchical multinomial logit
  or Dirichlet multinomial conditioned on the standing overstory composition (the seed
  source), site and climate (BGI or CSPI), disturbance state, and shade tolerance and seed
  dispersal traits. The count model says how many; the composition model says who; their
  product gives recruits by species.
- **Keep the hurdle / zero handling.** Many FIA plots have zero recruits over a short
  interval. The plain negbinom can under handle excess zeros at annual resolution; the
  hurdle (35b) or a zero inflated negbinom is the safer base. Compare by held out ELPD.
- **Top height as the recruitment clock is correct and points straight at section 3.**
  Recruitment falling as the stand grows in height and closes is exactly the age independent
  behavior Garcia models capture, so the ingrowth driver and the stand level density
  equation below should share the same top height and self thinning logic rather than be
  fit in isolation.

## 3. A stand level constraint layer, Garcia tradition (the main suggestion)

**Worth doing. It is a real gap and a high value one.** The system today is six individual
tree components plus a plot level ingrowth count. Summed over a long projection, small
biases in DG, HG, mortality, and ingrowth compound, and the aggregate stand trajectory
(total basal area, density, top height) can drift away from what stands actually do. A small
stand level model gives an independent, robust prediction of that aggregate trajectory that
you can use to discipline the tree level engine. This matters most for the carbon and yield
applications (PERSEUS, the CBM intercomparisons) where the stand aggregate, not the
individual tree, is the quantity of interest.

### 3.1 A minimal state vector and three transition functions

Follow Garcia's state space view: describe the stand by a small state vector and evolve each
element with an autonomous transition function `y2 = G(y1, t2 - t1, site)` that is age
independent and, because it is a transition on the interval, annualized and path invariant by
construction. Three states:

- **Top height H** (use the HT40 you already compute). A Chapman Richards GADA in algebraic
  difference form, `H2 = f(H1, t2 - t1)`, which removes age by letting H reference itself
  (Garcia 1983 for the stochastic height SDE, the GADA for the difference form). This is also
  the canonical site reading and ties directly to the CSPI and HG work: top height growth is
  the height based site signal, so H and the site surface should be mutually consistent.
- **Density N** (stems per ha above the merch threshold). A mortality and self thinning
  transition driven by H increment rather than age, bounded by the maximum size density line
  you already estimate (SDImax_brms, Reineke). N is the net of two processes the tree model
  handles separately: ingrowth adds from below, self thinning removes near the boundary, so
  the stand level N equation is the natural place to reconcile the negbinom recruitment with
  the tree level mortality.
- **Basal area G**. A growth transition driven by H (the clock), N, and the productivity
  index (BGI or CSPI), age independent, `G2 = f(G1, H, N, t2 - t1, site)`.

Fit all three to the same FIA remeasurement pairs aggregated to the plot, with the same
ecoregion and forest type random effects and the same productivity driver as the tree level
system, so the two levels are consistent rather than two separate stories.

### 3.2 How it constrains the tree level model

Start one way, move to two way once validated:

1. **Diagnostic (one way, low risk, do first).** Run the stand level equations alongside the
   tree level projection and flag where the summed tree predictions diverge from the stand
   trajectory. You already have `36_conus_benchmark.R` as a stand level FIA benchmark engine,
   so this is mostly wiring, and it immediately tells you where the individual tree model
   drifts.
2. **Constrained disaggregation (two way, higher value).** Use the stand level BA growth and
   density trajectory as the envelope, then disaggregate to trees with a proportional
   allocation that preserves the tree level distribution shape while scaling the sum to the
   stand total (Ritchie and Hann style disaggregation). The individual tree competition
   dynamics still set who grows and who dies; the stand equation sets the total. This is the
   classic way whole stand and individual tree models are reconciled and it stabilizes long
   horizon behavior without discarding tree level detail.

### 3.3 Refinements and cautions to build in from the start

- **Path invariance is the annualization guarantee. Test it explicitly:** ten one year steps
  must equal one ten year step. Garcia transition functions satisfy this by design; verify it
  numerically as a unit test so the annual engine and any multi year call agree.
- **Compatibility (Clutter 1963).** If you fit both a BA growth and a BA yield surface, make
  the growth equation the derivative of the yield equation so the two are mathematically
  consistent. The state space form gives this for free because it is one dynamical system;
  keep it that way rather than fitting growth and yield independently.
- **Self thinning consistency.** The stand level N trajectory must not cross the SDImax
  boundary, and the tree level mortality should agree with the stand level self thinning at
  that boundary. Tie the density transition to the existing SDImax_brms surface so both
  levels share one carrying capacity.
- **Reconcile ingrowth at the stand level.** The recruitment term in dN/dt is exactly the
  negbinom count, and the species composition model allocates it. Fit them to be consistent
  rather than letting plot level ingrowth and stand level density imply different recruitment.
- **Keep it small.** Three states, a handful of parameters each, nonlinear forms that
  extrapolate sensibly. Resist adding covariates the tree level model already carries; the
  stand layer earns its keep by being simple, robust, and hard to break, not by matching the
  tree model's resolution.

### 3.4 Is it worth it

Yes. The data, the productivity driver, the SDImax surface, the stand benchmark engine, and
the Bakuzis stand level uncertainty harness already exist, so the marginal cost is three
compatible equations fit to data in hand. The payoff is disciplined long horizon projection,
a principled home for self thinning and site limits, a natural reconciliation point for
ingrowth and mortality, and exactly the aggregate quantities the carbon and yield
applications consume. It also strengthens the manuscript narrative: a unified tree level
calibration constrained by a compact, age independent, annualized stand level model is a
cleaner and more defensible system than tree equations alone.

## 4. Suggested sequencing

1. Finish what is running (the two ΔLOO jobs, the v7 QRF merge) and lock v7. No new launches
   until those land.
2. Ingrowth: switch composition to the trait driven compositional model and settle the
   hurdle vs negbinom base by held out ELPD. Small, self contained.
3. Stand level: fit top height (GADA difference form) first, since it is the cleanest and
   anchors the other two and the site work. Then density (tied to SDImax) and basal area.
   Wire as a one way diagnostic against `36_conus_benchmark.R` before any constrained
   disaggregation.
4. Carry all three stand states through the same posterior draw uncertainty path as the tree
   components so the aggregate trajectory has credible bands.

## Key references

Garcia 1983 (stochastic differential equation height growth), Garcia 1994 and the GADA
(age independent, self referencing site equations), Clutter 1963 (growth yield
compatibility), Reineke 1933 and the size density boundary, Ritchie and Hann (disaggregation
of stand to tree).
