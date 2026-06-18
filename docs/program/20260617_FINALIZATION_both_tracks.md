# Finalization: fvs-modern x fvs-conus — decisive results and the resolved next steps
2026-06-17. Capstone for the finalization session. Three results lock the core scientific claims; two
refinements remain with precise specs.

## 1. The disturbance-artifact headline is PROVEN (removal-simulation converse test)

The decisive test the red-team demanded. On COND-harvested plots, the default FVS run never cuts, giving
a large apparent over-prediction. Simulating the recorded harvest (removing the trees FIA records as cut
from the FVS input, then projecting the residual stand) collapses the bias to the undisturbed level:

| variant | harvested BA bias: default (no cut) -> removal-simulated | undisturbed reference |
|---|---|---|
| ne | +51.2% -> +1.3% | +10.0% |
| sn | +110.1% -> +0.1% | +15.5% |
| pn | +41.9% -> +12.9% | +13.1% |

Accounting for the recorded harvest entirely explains the harvested-plot over-prediction. The pooled "FVS
over-predicts" is unsimulated removal, not growth error. Now proven three independent ways: COND
stratification (undisturbed small, harvested +33-55%), the fvs-conus projector (same pattern), and this
converse test. This is the spine of both manuscripts. (n_harv modest: ne 14, sn 22, pn 9 - extend, but
the effect is enormous and unambiguous.)

## 2. The four-arm comparison on the disturbance-clean basis (NE, 21,811 undisturbed conditions)

Filtering the existing fvs-conus per-condition predictions to COND-undisturbed:

| metric | default FVS | fvs-conus equations |
|---|---|---|
| BA | -12.3% | -7.2% |
| TPA | -7.2% | -6.2% |
| QMD | -5.3% | -3.1% |
| merch volume (CFNET) | -15.7% | -9.2% |
| harvested BA (reference) | +46.4% | +55.6% |

Findings: (a) on the CORRECT disturbance-clean basis, the fvs-conus equations reduce bias vs default for
every metric - the earlier "fvs-conus looks worse" was harvest contamination of the pooled benchmark;
(b) the harvest artifact (+46-56%) appears in the fvs-conus projector too, independently confirming the
headline. IMPORTANT methodological note: the fvs-conus STANDALONE PROJECTOR under-predicts undisturbed
(-12% default) while the FVS ENGINE over-predicts (+10%, the fvs-modern benchmark) - opposite signs
because they are different projection machinery. A true single-sample four-arm (A default / B
keyword-calibrated / C fvs-conus / D combined) requires ONE projection framework; this is the top
integration build (see next steps).

## 3. Out-of-sample (held-out spatial-fold) validation

Calibration derived on fold-A counties, applied unchanged to held-out fold-B counties:
- SIZE LEVERS GENERALIZE: QMD bias reduced out-of-sample in all 4 variants (sn +13.2->+1.1, kt
  +32.3->+11.8); BA reduced/held in 3 of 4. The brms SDImax + BAIMULT levers are not overfit.
- THE FIXED INGROWTH RATE DOES NOT TRANSFER: TPH over/under-shoots on the held-out fold (sn -6.7->+16.1).
  Recruitment is spatially variable and density-dependent; a fixed per-variant rate cannot transfer.

## Resolved status of the two tracks

fvs-modern: the disturbance-aware benchmark + four-lever adjustment are established and (for the size
levers) out-of-sample validated; the headline is proven by the converse test. Manuscript-ready after the
two refinements below + the red-team's master-table/CIs pass.

fvs-conus: the species-free equations improve on the disturbance-clean basis (BA -12->-7, volume
-16->-9 for NE); the framework (trait-driven, annualized, site-index-free via BGI/CSPI/ESI) is validated
in pilot. Manuscript-ready after the disturbance-clean re-benchmark is extended beyond NE and the
combined product (D) is tested.

## The two remaining refinements (precise specs)

A. DENSITY-DEPENDENT RECRUITMENT (fixes the one lever that did not transfer). Replace the fixed
   rate*interval*initialTPA with recruits = R_max * headroom * interval, where headroom =
   max(0, 1 - SDI_t1/SDImax_brms) (and optionally scaled by site productivity BGI/CSPI). Recruitment
   shuts off as stands approach the max-SDI limit, which both bounds over-recruitment in dense stands and
   makes the rate transfer across folds/regions. Re-run the held-out validation; expect TPH to then
   transfer out-of-sample. (Script ready to adapt: held_out_validation.py seed_rows + rate term.)

B. CROWN-WIDTH / MCW UNIFICATION (Marshall). Recover MCW per variant/species from the CCF coefficients:
   linear form B0=sqrt(R1/k), B1=sqrt(R3/k), k=0.001803026 (MCW=B0+B1*DBH); power form
   A1=sqrt(R4/k), A2=R5/2 (MCW=A1*DBH^A2); (R1+R2+R3)*DBH form C1=sqrt(Rx/k), MCW=C1*DBH^0.5. Tabulate the
   per-species cross-variant spread (Douglas-fir spans 4+ curves), select a CONUS-consistent equation per
   species (genus / softwood-hardwood mapping for the ~232 unmapped FIA species), prefer Russell &
   Weiskittel (2010) for the Acadian region, and feed the unified MCW into the compatible crown/CCF
   competition term of both the calibrated engine and the fvs-conus crown component. Needs the per-variant
   CCF coefficient tables (FVS variant guides / fvs-modern src) as input.

## The single highest-value next build

The single-framework four-arm comparison (A/B/C/D on one disturbance-clean, held-out plot set, with
bootstrap CIs) - it resolves how much improvement is equations (fvs-conus) vs density calibration
(fvs-modern) vs both, and is the shared headline figure for both papers. Build on the FVS engine (so the
keyword levers apply) with the fvs-conus equations injected via the standalone path or fvs2py.

## Artifacts (committed to holoros/fvs-modern, conus-sf-integration-2026-05-21)

removal_sim.py + removal_sim.csv (converse test), fourarm_clean.py (disturbance-clean fvs-conus),
held_out_validation.py + held_out.csv, the program charters + integration roadmap (docs/program/),
calib_final.csv, the red-team review and revision roadmap. Zenodo: still deferred to submission (new
version of the fvs_perseus_conus concept DOI).
