# Briefing for the USFS FVS team: CONUS-wide component equations

**From:** A. Weiskittel, Center for Research on Sustainable Forests, University of Maine, with collaborators
**Date:** 14 June 2026
**Purpose:** brief the FVS staff on a program to build a single set of CONUS-wide component equations,
where it stands against the regional variants and against Greg Johnson's CONUS equations, and where
collaboration would be most useful. A companion technical report covers maximum SDI in depth.

## 1. What we are building and why

A single CONUS-unified set of FVS component equations, calibrated on FIA remeasurement, intended to
predict short- and long-term stand development for all US species at least as accurately as the 20
regional variants, with calibrated uncertainty. The motivation is the same one the FVS staff knows
well: 20 separately maintained variants with different functional forms and species coverage are hard
to keep current, leave rare species thinly supported, and make consistent national analyses (carbon,
fuels, climate) awkward. A unified, trait-aware set addresses all three without discarding what the
regional variants do well.

Two design axes organize the work, and they are the useful way to think about it:

- **Tree-level vs stand-level.** Tree-level equations predict each tree's growth and survival;
  stand-level equations constrain the whole-stand trajectory (density, basal area, top height) so
  long-term projections respect the size-density limit.
- **Species-dependent vs species-independent.** A per-species leg where data are sufficient, and a
  trait-driven "species-free" leg keyed on functional traits and climate that covers rare and
  unsampled species. The two are blended per species by a shrinkage weight, so well-sampled species
  keep their own fit and rare species fall back gracefully to the trait form.

## 2. Where it stands

**Tree-level components are essentially built.** Diameter growth (Kuehne form, BGI-driven), height
growth (ORGANON form), height-diameter, height-to-crown-base, crown ratio, survival (gompit with
exposure), and ingrowth (count plus a trait-based species-composition model). Three of these carry a
validated refinement awaiting a production refit (annualized HCB and crown ratio, a relative-size
senescence term in survival).

**Stand-level constraint equations are prototyped** (age-independent, annualized, in the spirit of
García's state-space approach) but not yet integrated as a constraint layer. The density equation is
built on relative density, which is why maximum SDI (below) is central rather than a side topic.

**Benchmark verdict, first clean result.** On a clean Northeast benchmark (FIA remeasurement, year-0
reproducing observed stand state exactly), the unified equations already beat the native NE and ACD
variants on basal area: about -0.6 percent bias versus +12 to +13 percent for the regional variants.
Density was the one weakness, and we traced it entirely to the maximum-SDI setting, not the growth or
mortality equations.

**Versus Greg Johnson's CONUS equations (diameter growth, Douglas-fir).** Our species-dependent DG is
competitive on RMSE (0.091 vs 0.097 cm/yr), with Greg better calibrated on bias and explained variance.
The notable result is the species-free leg: predicting Douglas-fir from traits alone, having never seen
a Douglas-fir, it reaches RMSE 0.118, within about 22 percent of a model fit on 156,000 Douglas-fir
observations. That is the coverage case where Greg's species-specific system, and the regional variants
for rare species, have no model at all.

## 3. The maximum-SDI finding (companion report has the detail)

We think this is the most immediately useful result for the FVS staff. The FVS species-weighted maximum
SDI is biased about 28 percent high against an FIA-derived maximum and has almost no plot-level skill,
and a localized, FIA-derived maximum (from our TreeMap-based 30 m CONUS surface, Zenodo
10.5281/zenodo.19509367) predicts observed self-thinning about 85 percent better, nationally and in the
West.

The operational nuance matters and we want to be candid about it. We ran paired FVS projections, default
versus localized maximum, on three variants (Pacific Northwest, Central Rockies, Southern). In the
Pacific Northwest the localized maximum cut density error by about a quarter in dense stands; in CR and
SN the raw drop-in was neutral-to-harmful. A scale diagnostic resolved why: the FIA-derived maximum
carries a useful spatial pattern and a level, and once the level is calibrated (about 0.9 times the raw
value in CR) the localized maximum beats the native ceiling with near-zero bias. So the spatial
information is good everywhere; only the level needs to be made consistent with the variant's mortality
calibration. The recommendation is therefore specific: adopt the localized spatial pattern and calibrate
its level jointly with the density-dependent mortality response, which is automatic in a unified CONUS
fit and a one-parameter regional adjustment when retrofitting a native variant.

## 4. What is candidly not done

The engine integration is the largest remaining task. Today the benchmark drives the real FVS engine
through multipliers and the SDIMAX keyword; the trait-driven species-free equations are not yet injected
as the engine's growth functions. Until they are, the strongest claims rest on the offline component
comparisons (competitive with Greg, species-free within 22 percent) and the maximum-SDI result, all of
which are real and validated. We are not claiming a finished unified variant; we are claiming a
validated component set and a clear, bounded path to the engine.

## 5. Where the FVS staff's input would be most valuable

1. **Maximum SDI.** Whether the FVS team would consider the FIA-derived maximum-SDI surface as a
   recommended source for the density-limit value, paired with a check that each variant's mortality is
   calibrated to it. The surface and plot table are available now.
2. **Benchmark design.** The most useful variants and FIA stratification for the head-to-head, and
   whether the FVS team has a preferred held-out evaluation protocol we should match so results are
   directly comparable to FVS's own validation.
3. **Engine integration path.** The cleanest mechanism to inject CONUS component equations into the FVS
   engine for evaluation (keyword, alternate equation hooks, or a parallel build), so a unified variant
   can be tested inside FVS rather than only offline.
4. **Definitions.** Confirming the SDI convention (metric, summation) so relative density is computed
   consistently between our surface and FVS internals.

## 6. One-paragraph summary

We have a CONUS-wide component-equation set that is built at the tree level, prototyped at the stand
level, already beats the regional NE and ACD variants on basal area, is competitive with Greg Johnson's
diameter growth and remarkably close even with zero species data, and carries a validated, model-agnostic
maximum-SDI improvement whose one important caveat (consistency with the mortality response) we have
quantified. The remaining work is engine integration and a broader benchmark. We would value the FVS
team's guidance on maximum SDI, benchmark design, and the integration path, and we are glad to share
data, code, and the companion technical report.
