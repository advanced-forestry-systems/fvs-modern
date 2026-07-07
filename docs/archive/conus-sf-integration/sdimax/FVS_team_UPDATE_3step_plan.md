# A staged plan to refine FVS: update for the FVS team

**From:** A. Weiskittel and collaborators, Center for Research on Sustainable Forests, University of Maine
**Date:** 14 June 2026
**Purpose:** a single update on our FVS refinement work, framed as one plan with three steps that move
from near-term, deployable improvements to a long-term modernization. Each step stands alone and
improves FVS; together they modernize it end to end. A detailed technical report and the data are
available on request.

## The plan in one paragraph

We are pursuing three steps. **Step 1** layers FIA-derived calibration modifiers onto the existing FVS
engine, the basis of our fvs-modern work, and it already beats the regional Northeast and Acadian
variants on basal area. **Step 2** replaces the species-weighted maximum SDI with a localized,
FIA-derived one, jointly calibrated with mortality, to improve long-term density and self-thinning; we
have now validated this across all 20 CONUS variants. **Step 3** refits the component equations as a
single trait-aware CONUS variant covering all species, retiring the 20 separately maintained regional
variants; the components are built and engine integration is the next task.

## Step 1. Calibrate the engine with FIA modifiers (deployable now)

Keep the FVS engine and species tables, and layer FIA-derived calibration on top: mortality, growth,
and density-limit modifiers (MORTMULT, BAIMULT, SDIMAX) fit to FIA remeasurement. There are no new
equations to certify; it improves accuracy with the engine as-is. On a clean Northeast benchmark (FIA
remeasurement, year-0 reproducing the observed stand state exactly), the calibrated engine is nearly
unbiased on basal area (about -0.6 percent) where the native NE and ACD variants run +12 to +13
percent. This is the immediately useful, low-risk step.

## Step 2. Refine the maximum SDI (improves long-term predictions)

The density limit governs self-thinning over a rotation, so getting it right is central to long-term
accuracy. The FVS species-weighted maximum is biased about 28 percent high against an FIA-derived
maximum and has near-zero plot-level skill. Because maximum SDI is not observable, we tested it
predictively: relative density built from a localized, FIA-derived maximum predicts observed
self-thinning about 85 percent better than species-weighting, in every region (82,130 remeasured plots).

Sweeping the level applied to the localized maximum through the real FVS mortality response, for every
CONUS variant, gives two results. The level-calibrated localized maximum matches or beats the native
species-weighted maximum in 19 of the 20 variants, with the largest gains where the native level was
most off (Klamath Mountains 75 to 54 percent density RMSE, Acadian 53 to 37, Northeast 35 to 31). And
the optimal level is strongly variant-specific, spanning 0.6 to 2.0 with a median near 1.2. The lesson
is to set the density limit from a localized, FIA-derived maximum and calibrate its level jointly with
each variant's mortality, rather than dropping a common maximum in at a uniform level.

## Step 3. Refit the component equations as one CONUS variant (long-term)

The long-term goal is a single trait-aware CONUS variant that covers all species and retires the 20
regional variants. It is organized on two axes: tree-level versus stand-level, and species-dependent
(per-species fits) versus species-independent (a trait-driven, species-free form), blended per species
so well-sampled species keep their own fit and rare species fall back to traits. Early evidence is
encouraging: on Douglas-fir diameter growth our species-dependent fit is competitive with Greg
Johnson's species-specific model (RMSE 0.091 versus 0.097 cm/yr), and the species-free leg reaches RMSE
0.118 predicting Douglas-fir from traits alone, having never seen one, within about 22 percent, exactly
the rare-species case the regional variants do not cover. The components are built and the stand-level
constraint layer is prototyped; injecting the trait-driven equations into the engine is the next task.

## Status and where FVS-team input would help

| step | status |
|---|---|
| 1. FIA calibration modifiers | working and benchmarked; deployable on the current engine |
| 2. Maximum SDI | validated across all 20 variants; ready to adopt jointly with mortality |
| 3. CONUS variant | components built, stand-level prototyped; engine injection is next |

We would value the FVS team's guidance on three things: the recommended source for the maximum-SDI
value and the SDI convention to match so relative density is computed consistently; the cleanest path
to inject CONUS component equations into the engine for evaluation; and the benchmark protocol and
held-out evaluation you would want us to match so results are directly comparable to FVS's own
validation. We are glad to share the detailed technical report, the data, and the code.
