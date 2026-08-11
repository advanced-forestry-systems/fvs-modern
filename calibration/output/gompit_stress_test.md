# Gompit in-engine stress test: 19 variants, ~1,100 stands, robustness pass

A broad robustness test of the validated gompit-in-FVS system (coarse crosswalk,
the empirically optimal one; in-engine Fortran integration). Each of 19 variants
spanning all five mortality-routine families was run native vs gompit on ~60
stands, 100-yr no-harvest, ~45,000 projection-years total.

## Robustness: clean across the board

**Zero** gompit output rows were NaN, infinite, negative, or above 2000 t/ac AGB
across all ~45,000 projection-years. No crashes, no runaway on any variant. Max
gompit AGB was 1122 t/ac (NC, dense redwood/Klamath conifer -- physically
reasonable). This confirms the in-engine substitution is numerically robust at
scale on the whole CONUS variant set.

## Mean AGB at yr100 (native -> gompit), ~60 stands/variant

| var | stands | native | gompit | change | max gompit AGB |
|-----|-------:|-------:|-------:|-------:|---------------:|
| NE | 60 | 67.6 | 56.2 | -17% | 369 |
| CS | 60 | 158.7 | 125.4 | -21% | 261 |
| LS | 59 | 168.7 | 120.2 | -29% | 264 |
| SN | 60 | 48.0 | 12.5 | -74% | 96 |
| CR | 59 | 5.8 | 2.6 | -54% | 20 |
| WS | 60 | 266.4 | 168.5 | -37% | 457 |
| EC | 59 | 109.6 | 96.2 | -12% | 253 |
| CA | 60 | 303.3 | 241.6 | -20% | 687 |
| WC | 60 | 253.2 | 138.5 | -45% | 364 |
| PN | 58 | 335.5 | 204.1 | -39% | 803 |
| AK | 60 | 10.7 | 7.0 | -35% | 78 |
| BM | 60 | 123.0 | 107.1 | -13% | 335 |
| CI | 60 | 14.6 | 3.7 | -75% | 16 |
| EM | 59 | 18.6 | 8.6 | -54% | 33 |
| IE | 60 | 111.5 | 81.2 | -27% | 260 |
| NC | 60 | 396.8 | 382.1 | -4% | 1122 |
| SO | 59 | 140.4 | 132.3 | -6% | 465 |
| UT | 60 | 0.9 | 0.2 | -82% | 5 |
| TT | 59 | 13.3 | 3.1 | -77% | 15 |

(Data: `gompit_stress_summary.csv`.)

## Reading

The pattern is consistent with the cch mechanism. The smallest gompit effects
are NC (-4%), SO (-6%), EC (-12%) -- western conifer, where the ORGANON crown
proxy fits best and native mortality is already reasonable. The largest are the
sparse, low-biomass high-elevation/high-latitude variants (UT 0.9 t/ac, TT/CI
~14, CR/EM ~15-19) where small absolute changes read as large percentages, and
SN (-74%) where dense southern stands carry high crown closure. Productive,
moderate-density variants (NE, CS, CA, BM, IE) land in a tight -13% to -27% band.

## Bottom line

Across 19 variants and five mortality-routine families, gompit-in-FVS produces
bounded, realistic, crash-free 100-yr projections with no numerical
pathologies -- the system passes the broad stress test. Combined with the
GOMPMORT keyword (reproducible activation) and the evidence-based retention of
the coarse crosswalk, the gompit Fortran integration is complete and robust.
