# FIA -> ORGANON crown-group crosswalk: refinement tested, coarse proxy retained

The gompit projection computes cch (crown closure at tree tip) with an ORGANON
crown-closure port, then affine-maps it onto the gompit cch scale
(cch = CCH_A + CCH_B*cch_hat). The crosswalk that assigns each FIA species to an
ORGANON crown group was a coarse softwood/hardwood proxy (FIA<300 -> group 1 DF,
else group 16 RA). The variant-specific gompit effect (EC/SO/NC ~0 vs SN -78)
raised the question of whether a finer, genus/crown-form crosswalk would improve
the cch fit. This was tested directly.

## Test

`calibration/python/refine_cch_crosswalk.py` builds a genus/crown-form crosswalk
over all 18 ORGANON SWO groups, recomputes cch_hat on the held cch validation
sample (113k trees, 4000 plots, the same set used to fit the original affine
map), and compares Spearman correlation with the panel's stored CCH1.

| crosswalk | n | Pearson | Spearman |
|-----------|--:|--------:|---------:|
| **coarse (softwood/hardwood)** | 113,465 | 0.844 | **0.925** |
| conifer-refined, hardwood coarse | 110k | 0.762 | 0.875 |
| full genus (18 groups) | 108,672 | 0.690 | 0.853 |

## Result: the coarse proxy is empirically best; no change adopted

Both refinements **degrade** the cch reproduction. The uniform softwood/hardwood
proxy reproduces the stored CCH1 better than genus-specific ORGANON crowns,
because the ORGANON crown-width equations are Pacific-Northwest-specific and,
applied to the full CONUS species set (especially ~140 eastern hardwoods forced
onto group-14 Oregon white oak geometry), introduce species-to-species crown
variation that does not match how CCH1 was generated. For a cch proxy whose role
is to rank-order crown closure before an affine map, uniformity wins.

Two consequences:

1. The validated coarse proxy (GGRP in `gompmort.f90`; CCH_A=0.062, CCH_B=0.0036)
   is retained unchanged. It is the empirically optimal choice, not a limitation
   -- which justifies the simple proxy in the manuscript.
2. The variant spread in gompit effect (EC ~0 ... SN -78) is therefore **not** a
   crosswalk artifact. It is a genuine property of the gompit model's cch
   sensitivity across stand structures (dense wet PNW and dense southern stands
   have high crown closure -> larger gompit mortality; dry open stands less).

The flagged "group-map refinement" item is thus resolved on the evidence: tested,
rejected, coarse proxy kept.
