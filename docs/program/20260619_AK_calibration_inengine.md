# FVS-Alaska in-engine growth calibration via BAIMULT
2026-06-19. Closes the AK calibration: applies the growth correction in the engine and confirms the bias drops.

## The keyword mechanism (resolved)

The BAIMULT keyword appeared inert in earlier tests. Reading the FVS source (base/mults.f90, base/keywds.f90)
showed why: MULTS reads the schedule cycle from a separate field, so BAIMULT needs three fields, not one:

```
BAIMULT  <cycle>  <species>  <multiplier>
```

with cycle 0 = all cycles, species 0 = all species. The earlier one-field line placed the multiplier in the
species slot, so FVS applied no multiplier. With the corrected three-field form (cycle 0, species 0, M) the
multiplier scales the large-tree diameter increment as intended; effect saturates at high M because dense
stands hit the self-thinning (max SDI) limit.

## Calibration result (BC MAGPlot, by NA Level I ecoregion)

Default FVS-AK vs the best whole-number BAIMULT, basal-area increment bias:

| NA Level I ecoregion | n | default increment bias | best BAIMULT | calibrated bias |
|---|---|---|---|---|
| Marine West Coast Forest (AK analog) | 11 | -49.4% | 2x | -6.4% |
| Taiga | 16 | -73.7% | 2x | -17.4% |
| Northern Forests | 17 | -41.3% | 2x | +18.7% (2x overshoots) |
| Northwestern Forested Mountains | 51 | -13.2% | 1x | -13.2% |
| North American Deserts | 21 | +39.7% | 1x | +39.7% |

Findings:

- The Alaska variant calibrates well in its own analog ecoregion: a 2x basal-area-increment multiplier reduces
  the Marine West Coast Forest bias from -49% to -6%, essentially removing the under-prediction.
- The correction is ecoregion (productivity) dependent, as the validation predicted. Productive coastal and
  boreal forest (Marine West Coast Forest, Taiga) want about 2x; dry interior (North American Deserts,
  Northwestern Forested Mountains) want about 1x, and a blanket 2x would over-correct them (Northern Forests
  flips from -41% to +19% at 2x). A single global multiplier is the wrong model; the AK calibration should be
  applied per ecoregion or per productivity class.
- Per-ecoregion point estimates carry sampling noise (n = 11 to 51, clean-ingestion subset, and the default
  Marine West Coast Forest bias is -49% here against -74% in the validation sample). The robust, repeatable
  conclusion is the direction and the rough 2x-coastal / 1x-dry split, not a precise multiplier per cell.

## Status

In-engine AK calibration is demonstrated end to end: keyword mechanism fixed, multiplier applied through the
standalone engine, bias measured before and after. The production step is to encode the ecoregion-dependent
multiplier (about 2x for Marine West Coast Forest and Taiga, 1x for dry interior) into config/calibrated/ak.json
and the keyword pipeline, then validate on a held-out MAGPlot fold. The current ak.json carries near-1.0
diameter-growth multipliers (calibrated against FIA Alaska), which do not reflect the coastal-BC
under-prediction; this MAGPlot result is the basis for updating it.

Scripts: diagnostics_2026-06-16/magplot/ak_calib_par.py, ak_calib_results.csv.
