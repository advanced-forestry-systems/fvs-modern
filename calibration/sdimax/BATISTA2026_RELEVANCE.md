# Relevance of Batista et al. (2026) to the CONUS FVS program

**Paper:** Batista, Birch, Dickinson, Hoffman, Lutz, Miesel (2026). Benchmarking and calibrating FVS
diameter growth predictions with tree-rings and forest inventory data in Sierra Nevada mixed-conifer
forests. *Forest Ecology and Management* 618:123981. Online 11 June 2026.
**Date of this note:** 2026-06-15

## One line

An independent group, using ground-truth tree rings, reaches the same conclusions our program reaches
and recommends the same fix, which is strong external corroboration of all three steps.

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

This paper corroborates the three-step framing end to end and recommends our Step 3 explicitly. It
belongs in the FVS-team materials as external, independent, very recent support, and it strengthens the
case that systematic structure-dependent bias in FVS variants is general and that the fix is
recalibration now plus flexible size-dependent equations next.
