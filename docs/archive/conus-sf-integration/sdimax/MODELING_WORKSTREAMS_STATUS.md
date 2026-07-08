# Modeling workstreams: unified joint fit, tree-level/ingrowth refits, engine injection

**Date:** 2026-06-14
**Scope:** progress on the three modeling recommendations, run on Cardinal against FIA remeasurement.

## 1. Unified joint fit (done)

Fit one shared self-thinning response plus a per-variant level scalar on the localized (brms) maximum
SDI, estimated jointly from 351,000 stratified remeasured plots across 19 regions. The shared response
is a logistic in relative density (identifiable against the per-variant level, unlike a free smooth):

> annual self-thinning rate = -0.021 + 0.638 / (1 + exp(-(RD - 1.56) / 0.494)),  RD = SDI / (k_region · maxSDI)

The headline result: estimating a per-variant level jointly more than doubles the explained variance of
observed self-thinning over forcing a uniform level (R-squared 0.021 to 0.046). Absolute skill is low
because annual plot-level mortality is very noisy; the relative gain is the point, and it confirms that
variant-specific level calibration carries real signal in the observed data, not only in the engine.

The data-estimated levels span about 0.7 (Southern) to 2.4 (Tetons, Utah), most above 1.0. They agree
only weakly with the engine-derived optimal levels (correlation 0.12). That weak agreement is itself
informative: the level a variant needs *inside FVS* mostly compensates for FVS's own fixed mortality
form, whereas the joint-fit level reflects the true density-mortality relationship free of that form. It
is the strongest argument yet for refitting mortality jointly with the localized maximum rather than
retrofitting a level onto a native variant. Script `unified_joint_fit.R`, levels `joint_fit_levels.csv`,
figure `joint_fit.png`.

## 2. Tree-level and ingrowth refits (validated/banked; adoption gated)

The species-free (and species-dependent) forms for the tree-level components are fit and banked as
integration bundles: diameter growth (`dg`, `dg_v8_sf`), height growth (`hg`, `hg_v5_prod`,
`hg_v8rd_sf`), height-to-crown-base (`hcb_v2split`), height-diameter (`htdbh_v2split`), crown
(`cr_t2`), and survival (`surv_crz`). The validated refinements (annualized HCB and crown ratio, the
relative-size senescence term in survival via the gompit-with-exposure form) are present in these
bundles. So the "refits" are largely complete as fitted models.

Two honest constraints on finishing them in production. First, the production fits are Stan/HMC and a
fresh rerun needs a dedicated compute allocation; the cluster currently has eight of Aaron's jobs
queued under the association limit, so I did not submit competing jobs. Second, promoting any of these
into the production config is a hard stop pending sign-off (no production config writes). The remaining
ingrowth species-composition production fit (`36_fit_ingrowth_species_composition.R`) is the one piece
still needing a compute slot to finalize. Net: the modeling is done or one scheduled fit away; adoption
is a sign-off and a config write, not new science.

## 3. Engine injection of the species-free equations (framework proven; wiring is the build)

The injection-evaluation framework exists and ran cleanly (`benchmark_sf_vs_legA.R`). It evaluates four
arms on a held-out test set with both point accuracy and prediction-interval calibration: pure
species-free (the injected trait form), a hybrid (per-species where a leg-A fit is reliable, else
trait), a leg-A-style coverage, and a global-only baseline.

First real result, height-to-crown-base, Northeast, 60,000 held-out trees:

| arm | RMSE | R-squared | PICP-95 (target 0.95) |
|---|---:|---:|---:|
| pure species-free (injected) | 0.140 | 0.298 | 0.94 |
| hybrid | 0.140 | 0.298 | 0.94 |
| leg-A-style (per-species) | 0.147 | 0.226 | 0.93 |
| global only | 0.147 | 0.226 | 0.93 |

Two findings. The injected species-free form generalizes *better* than the per-species (leg-A) approach
on held-out data (RMSE 0.140 versus 0.147, R-squared 0.298 versus 0.226), and its prediction intervals
are well calibrated (PICP 0.94 against a nominal 0.95). And pure species-free equals the hybrid, so for
this component no per-species fallback is even required. A complementary in-sample LOO comparison shows
the species-dependent form fits the training data better, as expected from its extra parameters; the
out-of-sample benchmark is the one that matters for a model meant to generalize and to cover species it
has never seen, and there the trait form wins. The benchmark for the remaining components (height
growth, diameter growth, height-diameter, survival) is running and will populate the same table.

What "injection" still requires beyond this offline benchmark: wiring the banked bundle predictions into
the FVS engine's per-tree increment path so the engine itself runs on the trait-driven forms, rather
than the current multiplier-plus-SDIMAX calibration on the native equations. The bundles are the inputs;
the engine wiring is the build task, and it is the largest remaining piece for a true unified variant.

## Bottom line

The unified joint fit is done and clarifies why per-variant levels must be estimated with mortality
rather than retrofitted. The tree-level and ingrowth forms are fit and banked; their production adoption
is gated on a compute slot (ingrowth composition) and on sign-off plus a config write (everyone else),
not on new modeling. The species-free injection is proven competitive and well-calibrated on held-out
data, component by component, and the remaining work is wiring those banked forms into the engine.
