# Relevance of Batista et al. (2026) to the CONUS FVS program

**Paper:** Batista, Birch, Dickinson, Hoffman, Lutz, Miesel (2026). Benchmarking and calibrating FVS
diameter growth predictions with tree-rings and forest inventory data in Sierra Nevada mixed-conifer
forests. *Forest Ecology and Management* 618:123981. Online 11 June 2026.
**Date of this note:** 2026-06-15

## One line

An independent group reaches the same qualitative conclusions and recommends our Step 3, which is
useful directional corroboration; but its diameter-growth benchmark rests on tree rings, a biased
growth reference, and this program already calibrated diameter growth on FIA remeasurement, the better
reference, so it adds nothing quantitative to the growth calibration already in place.

## Two corrections to an over-enthusiastic first read

1. **Tree rings are a biased measure of tree growth.** Increment cores come from trees that survived to
   be measured, so the sample over-represents fast-growing large survivors and underrepresents the
   suppressed trees that died small. Back-calculating diameter from survivor cores reconstructs a
   growth trajectory that is systematically too high for the population (the survivorship, or
   large-tree selection, bias; e.g. Brienen et al. 2012). The benchmark in this paper is therefore
   biased upward, and the reported 41 percent underprediction is partly an artifact of that, not a
   clean property of FVS-WS. It should not be read as a quantitative magnitude.
2. **We already calibrated FVS diameter growth on FIA remeasurement.** FIA remeasurement directly
   measures diameter change on every tagged tree over the interval, including trees that die. It is
   population-representative and actually measured, so it is the less-biased reference for the growth
   channel, and the calibration is already done here. The external study uses a weaker reference for a
   step we have completed.

So the paper is directional support, not validation, and not new for the growth channel.

## What they did and found

They benchmarked the FVS Western Sierra Nevada (WS) variant diameter growth against tree-ring
reconstructed DBH plus FIA over 40 years, 1,016 trees in 128 plots, four conifers. The default WS
variant underpredicts diameter growth by about 41 percent (mean bias -40.6 percent, RMSE 17.2 cm, 66
percent relative), worst for large trees and in dense, high-QMD stands. A species-specific
multiplicative calibration (Random Forest multipliers, 0.6 to 2.3) cut RMSE to 12 percent, bias to -3.4
percent, and raised the share of predictions within 15 percent from 10 to 68 percent. The error is
driven by tree size and stand structure (DBH, QMD, density), not by climate or topography.

## Why it matters to us, step by step

- **Step 1 (FIA calibration modifiers).** This is independent validation, with ground-truth tree rings,
  of exactly our near-term step: native FVS variant growth is badly biased, and a multiplicative
  calibration (their RF multipliers, our GROWMULT/BAIMULT fit to FIA) corrects it dramatically. Their
  sevenfold gain in within-15-percent accuracy is the same kind of result our calibrated engine shows
  on basal area. It is the strongest external evidence yet that Step 1 is real and deployable.
- **Step 2 (density limit).** They find the error scales with QMD and stand density, the same
  structural axis our maximum-SDI and self-thinning work targets, and they flag dense, high-QMD,
  late-successional stands as where the default is most wrong. That is precisely where a correct
  density limit matters most.
- **Step 3 (refit the equations).** Their explicit conclusion is that calibration helps but that "more
  flexible size-dependent equations or refined competition indices may further enhance model
  performance." That is our Step 3 thesis stated by an independent group. Our trait-driven species-free
  diameter growth (competitive with Greg Johnson, structure-based: DBH, BAL, BA, CR) is exactly the
  flexible size-dependent form they call for.

## Specific connections

- The WS variant is one of the 20 we swept; it was among our higher-error variants, consistent with
  their finding that WS is poorly calibrated. Their result is diameter growth specific; ours was the
  density (self-thinning) channel, so the two are complementary views of the same miscalibrated variant.
- Their error structure (size and stand density dominate, climate and topography secondary) aligns with
  our structure-based species-free growth form, and is a useful caution against over-weighting climate
  for the diameter-growth channel specifically.
- They cite the same cluster of FVS-bias studies (Bagdon et al. 2021, Canavan and Ramm 2000, Dickinson
  et al. 2019, Ex and Smith 2014, Herbert et al. 2023, Leites et al. 2009, Pokharel and Froese 2008),
  plus Russell, Weiskittel and Kershaw 2013 and Giebink et al. 2022, so this program already sits in
  their reference frame.

## Methodological addition worth adopting

Their tree-ring reconstruction is a ground-truth validation for diameter growth that our FIA
remeasurement benchmark does not have. Pairing tree-ring increment data with the species-free diameter
growth equation, where cores are available, would give an independent, higher-resolution check on the
growth channel, especially for large trees where they show the default is most biased and where FIA
remeasurement intervals are coarse. This is a concrete future validation, not a change to the current
plan.

## Net

The paper is directional support: it shows another group sees FVS variants as systematically
structure-biased, endorses the multiplicative-calibration philosophy, and explicitly calls for flexible
size-dependent equations (our Step 3). It is not quantitative validation, because its diameter-growth
benchmark rests on a biased growth reference (survivor tree rings), and it does not add to the
diameter-growth calibration this program already completed on FIA remeasurement, which is the
less-biased and population-representative reference. In the FVS-team materials it is cited as related
work with these qualifications, not as independent validation.
