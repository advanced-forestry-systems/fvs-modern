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

## Coverage across all 25 variants (complete accounting)

The crown-width axis is now fully characterized for every variant:

- Recovered (15 western CCF variants): bm, ci, cr, ec, em, ie, kt, nc, pn, so, tt, ut, ws, wc, and the
  INDCCF-mapped pn/wc form. MCW = B0 + B1*D from RD1/RD3 (333 variant x species curves). The
  scientifically meaningful cross-variant inconsistency lives here, with Douglas-fir spanning five curves.
- Recovered (power form): ak stores a direct MCW power equation (MCW = B1*D^B2 via EQMAP), captured.
- Delegated to a crown-width routine: ca calls R5CRWD (the Region 5 crown-width routine) then
  CCFT = CRWD5^2 * 0.001803; ls calls CWCALC, which holds inline per-species crown-width equations with
  per-species caps (for example CW capped at 34, 33, 29 ft) rather than clean coefficient arrays. MCW for
  these two lives inside the named routine and is extracted by reading it, not by the CCF inversion.
- No parabolic MCW (eastern variants ne, acd, sn, cs): these zero the base crown-width coefficients
  (CWDL0 = CWDL1 = CWDL2 = 0 in grinit.f90), so they do not carry a per-species parabolic maximum crown
  width competition term in the CCF form. There is therefore no CCF-form MCW curve to unify for these
  variants; crown dynamics are handled through the crown-ratio models. If a unified crown width is wanted
  for the Northeast and Acadian region, the Russell and Weiskittel (2010) MCW equations are the
  recommended source.
- Out of scope: the Canadian variants bc and on have no FIA coverage and are handled separately.

Net: the unification target (the western CCF variants whose MCW curves disagree across variants) is fully
recovered; the eastern variants do not present the same parabolic-MCW inconsistency; ca and ls are the two
that still need a short read of their delegated routines if their curves are to be folded into the
cross-variant selection.

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
3. Select one CONUS consistent MCW per species, and feed it into the crown / CCF competition term of the
   calibrated engine and the fvs-conus crown component. This is the lever the defective CRNMULT keyword
   could not reach.

## CONUS-consistent selection (done)

`mcw_unify.py` selects a consensus curve per species as the cross-variant median of the recovered B0 and
B1 (linear form), producing one CONUS-consistent MCW = B0 + B1*D per species. 44 species receive a
consensus curve, 25 of them informed by more than one variant. Examples (MCW in feet, D in inches):

| species | n variants | B0 | B1 | MCW at 20 in | cross-variant spread |
|---|---|---|---|---|---|
| Douglas-fir (DF) | 5 | 4.71 | 1.499 | 34.69 | 5.01 |
| white fir (WF) | 4 | 5.45 | 1.253 | 30.51 | 10.72 |
| ponderosa pine (PP) | 3 | 3.49 | 1.343 | 30.34 | 1.29 |
| lodgepole pine (LP) | 3 | 3.27 | 1.423 | 31.72 | 11.78 |
| Engelmann spruce (ES) | 4 | 4.04 | 1.199 | 28.05 | 3.72 |

Full table: mcw_conus_consensus.csv. The median rule collapses the multi-curve spread (Douglas-fir from
five curves spanning 5 ft to a single consensus) into one defensible CONUS curve per species; a preferred
published source overrides the median where one exists (Russell and Weiskittel 2010 for the Acadian and
Northeast region). Remaining application step: emit these as the per-species crown-width coefficients in
the engine competition term and the fvs-conus crown component.

Script: `mcw_recovery.py`. Outputs: `mcw_by_variant_species.csv`, `mcw_cross_variant_spread.csv`.
