# fvs-conus performance stress test: species-specific versus species-free component equations

Date 2026-06-18. Branch holoros/fvs-modern conus-sf-integration-2026-05-21. Autonomous OODA run.
Question: what does fvs-conus gain or lose by dropping the species-specific component equations for the
species-free, trait-driven formulation. Answered at two levels: the component fits (how much per-species
deviation each component carries) and the stand projection (whether that deviation changes stand outcomes).

## Method

Each fvs-conus component is a hierarchical Bayesian fit with a species random intercept on top of the
trait-driven, site-resolved (CSPI) mean. The species-specific contribution is the species random-effect SD,
sigma_sp, on the linear-predictor scale (log link for the size and growth components, logit for mortality).
The species-free formulation sets that term to zero and predicts from traits alone. The typical
per-species deviation on the response scale is exp(sigma_sp) - 1 (a multiplicative effect; on odds for
mortality). A component's species term is well identified when sigma_sp mean / sd is at least 2, and its
species effects are real when a meaningful fraction of species have an intercept credibly different from
zero (|mean| / sd > 2). Inputs are the banked production fits (cspi_traits1 species intercepts and fixed
summaries) in fvs-conus/output/conus. Reproduce with species_stress.R.

## Component-level result

| component | n species | sigma_sp | per-species deviation | identified (mean/sd) | species credibly nonzero | verdict |
|---|---|---|---|---|---|---|
| Mortality | 36 | 0.424 | +52.8% on odds | 7.1 | 78% | keep species-specific |
| Crown recession | 24 | 0.258 | +29.4% | 5.1 | 58% | keep species-specific |
| Diameter growth | 75 | 0.207 | +23.0% | 10.9 | 41% | keep species-specific |
| Height growth | 9 | 0.160 | +17.3% | 2.2 | 0% | species-free adequate |
| Height-diameter | 8 | 0.269 | +30.8% | 1.1 | 0% | species-free adequate (term not identified) |
| Height to crown base | 8 | 0.132 | +14.1% | 1.1 | 25% | species-free adequate |

![Figure 1. Typical per-species deviation captured by the species-specific term, by component. Green is a well identified species term (sigma_sp mean/sd at least 2); orange is weakly identified.](fig_species_stress.png)

Reading. The species term earns its place for three components. Mortality carries the largest and best
identified species signal: a 53 percent per-species spread on the odds, with 78 percent of species credibly
different from the trait mean. Diameter growth is the most sharply identified term (mean/sd = 10.9) at a
23 percent per-species spread, and crown recession is similar at 29 percent. For the other three the case
is weak: height growth and height-diameter have zero species credibly nonzero, and the height-diameter
sigma_sp is barely identified (mean/sd = 1.1, the SD nearly equals the mean), so its 31 percent figure is
noise rather than a supported effect. Height to crown base is small and weakly identified. For those three,
the species-free formulation loses little.

## Projection-level result (from the committed four-arm)

The component spread does not by itself say what happens to a stand projection. The single-framework
four-arm benchmark (fourarm_abcd, 21 variants, out-of-sample) ran the engine with and without the
per-species growth signal expressed as keyword multipliers. At the stand level the species-specific growth
signal did not improve the median bias and slightly worsened it:

| metric | without species growth (default) | with species growth (emulated) |
|---|---|---|
| Basal area | 6.1 | 6.7 |
| Trees per hectare | 16.3 | 17.2 |
| QMD | 12.5 | 15.1 |
| Merch volume | 12.0 | 13.0 |

Median absolute percent bias. The per-species growth deviations that are real at the tree level largely
cancel when aggregated to the stand, and on disturbance-clean undisturbed plots the default growth is
already near unbiased, so rescaling per species mostly adds noise to the stand total.

## Integrated verdict

The two levels are consistent and complementary. Species identity captures real, well identified
per-tree behavior for mortality, diameter growth, and crown recession, and little that is supported for
height growth, height-diameter, and height to crown base. But those per-species size and growth deviations
wash out at the stand level, so the species-free formulation reproduces stand basal area, density,
diameter, and volume about as well as the species-specific one, while gaining parsimony and
transferability to unmeasured species.

Recommendation. Use the species-free, trait-driven components as the default for CONUS-wide stand
projection, where they are competitive and transfer to the roughly 232 FIA species without a fitted
curve. Retain the species-specific term where the target is the per-species process itself, foremost
mortality (53 percent odds spread, 78 percent of species supported) and secondarily diameter growth and
crown recession, since those drive species composition and stand dynamics that the stand-total metrics
hide. Drop the species term from height growth, height-diameter, and height to crown base, where it is not
identified.

[DATA_STATE]: six component fits read; sigma_sp and identifiability extracted per component; projection-level four-arm medians joined.
[OUTCOME_VERIFICATION]: species term well identified (mean/sd>=2) for mortality, CR, DG, HG; not identified for HD, HCB. Stand-level four-arm shows no median gain from species-specific growth (QMD 12.5 vs 15.1).
[IMPACT_UTILITY]: manuscript-preparer (fvs-conus v2 model-selection section), data-curator (component model cards). Defer Zenodo until OOS-honest.
[NEXT_AUTONOMOUS_STEP]: re-fit height growth, height-diameter, and HCB as species-free and confirm held-out ELPD is within 2 SE of the species-specific fits, formalizing the drop.
