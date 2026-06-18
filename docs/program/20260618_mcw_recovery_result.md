# Crown width / MCW recovery from variant CCF coefficients (Marshall)
2026-06-18. Item 3 of the integration roadmap. Recovers maximum crown width (MCW) from each FVS variant's
CCF coefficients so a CONUS consistent MCW per species can be selected and fed into the crown / CCF
competition term of both the calibrated engine and the fvs-conus crown component.

## Method

FVS stores per tree CCF as `CCFT = RD1 + RD2*D + RD3*D^2` (D in inches, large tree branch). Under the
open grown identity `CCF = k*MCW^2` with `k = 0.001803026` and MCW in feet, the linear form gives

    B0 = sqrt(RD1/k),  B1 = sqrt(RD3/k),  MCW = B0 + B1*D.

Power form variants store `CCFT = RDA*D^RDB`, giving `A1 = sqrt(RDA/k)`, `A2 = RDB/2`, `MCW = A1*D^A2`.
Alaska stores a direct MCW power equation `MCW = B1*D^B2` mapped by EQMAP. All three are handled by
`mcw_recovery.py`, which parses each `src-converted/<variant>/ccfcal.f90` and applies the formulas.

A cross term consistency ratio `RD2 / (2*sqrt(RD1*RD3))` is reported per row. Values near 1.0 confirm the
CCF quadratic is close to a perfect square, so the linear MCW recovery is faithful; departures (e.g. cr
Douglas fir at 1.06) flag species where the published CCF curve is not exactly a squared linear MCW and
the recovered B0/B1 are an approximation of the endpoints.

## Coverage

15 of 17 variants with a `ccfcal.f90` parse directly (333 variant x species curves): ak, bm, ca?, ci, cr,
ec, em, ie, kt, nc, pn, so, tt, ut, ws, wc. Two delegate CCF to a separate crown width routine and are not
yet captured here:

- `ca` calls `R5CRWD` then `CCFT = CRWD5^2 * 0.001803`.
- `ls` calls `CWCALC` then `CCFT = 0.001803 * TEMCW^2`.

For both, MCW is computed directly inside the named routine (R5CRWD, CWCALC) rather than from coefficient
arrays in ccfcal; read those routines to extract the curve. The four eastern variants without a variant
`ccfcal.f90` (ne, acd, cs, sn) inherit the base crown width path; the Acadian region should use the
Russell and Weiskittel (2010) MCW equations per the roadmap.

## Headline: cross variant spread per species (MCW at 20 in DBH)

The same species carries materially different open grown crown widths across variants, which is the
inconsistency the unification targets:

| species | n variants | MCW20 min (ft) | MCW20 max (ft) | range (ft) | variants |
|---|---|---|---|---|---|
| MH (mountain hemlock) | 4 | 28.05 | 43.88 | 15.83 | cr;pn;so;wc |
| AF (subalpine fir) | 4 | 18.85 | 34.05 | 15.20 | cr;pn;so;ut |
| LP (lodgepole pine) | 3 | 31.72 | 43.50 | 11.78 | pn;so;ut |
| WF (white fir) | 4 | 26.33 | 37.05 | 10.72 | cr;pn;so;ut |
| DF (Douglas fir) | 5 | 31.78 | 36.79 | 5.01 | cr;pn;so;tt;ut |

Douglas fir spans five distinct curves across variants, confirming the Marshall note. Full detail in
`mcw_by_variant_species.csv` (per curve) and `mcw_cross_variant_spread.csv` (per species, sorted by range).

## Next steps

1. Capture ca and ls by reading R5CRWD and CWCALC; add eastern variants (ne, acd, cs, sn) from the base
   crown width path, with Russell and Weiskittel (2010) for Acadian.
2. Map the variant species abbreviations to FIA SPCD; for the ~232 FIA species with no direct FVS curve,
   assign by genus, then by softwood / hardwood.
3. Select one CONUS consistent MCW per species (candidate rule: the cross variant median curve, or a
   preferred published source where one exists), and feed it into the crown / CCF competition term of the
   calibrated engine and the fvs-conus crown component. This is the lever the defective CRNMULT keyword
   could not reach.

Script: `mcw_recovery.py`. Outputs: `mcw_by_variant_species.csv`, `mcw_cross_variant_spread.csv`.
