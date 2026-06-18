# FVS modernization program: state of play
2026-06-18. Single source of truth for the two integrated tracks (fvs-modern, fvs-conus). Supersedes the
scattered handoff files for current status. Branch: holoros/fvs-modern conus-sf-integration-2026-05-21.
Compute: OSC Cardinal, account PUOM0008, user crsfaaron.

## 1. The program in one paragraph

fvs-modern recalibrates the existing FVS engine and parameters (including a brms site-specific maximum
stand density index) with uncertainty; fvs-conus fits CONUS-wide species-free, trait-driven, annualized
equations. The settled scientific result is that the widely cited FVS over-prediction is mostly a
benchmark-design artifact of unsimulated harvest, not a growth-equation bias, proven three independent
ways. The two tracks are complementary: the fvs-modern keyword calibration fixes the size and density
metrics (quadratic mean diameter, trees per hectare), and the fvs-conus equations fix the level and
scatter metrics (basal area, volume). The work is framed as a disturbance-aware benchmark plus a prototype
adjustment layer, validated out-of-sample for the size levers.

## 2. What is established (do not redo)

- Disturbance artifact, three ways: COND stratification (undisturbed median basal area bias +1.8 percent
  across 22 variants, pooled +14, harvested +42); removal-simulation converse test (harvested bias
  collapses to undisturbed when recorded harvest is simulated: ne +51.2 to +1.3, sn +110.1 to +0.1, pn
  +41.9 to +12.9); and the fvs-conus projector showing the same pattern.
- Size and density levers transfer out-of-sample. On spatially held-out folds across eight variants, the
  keyword adjustment cut median absolute quadratic mean diameter bias from 15.7 to 2.2 percent, basal area
  11.3 to 7.6, net merch volume 13.4 to 8.3, all with bootstrap CIs.
- Density-dependent recruitment removes the fixed-rate failure: the Southern out-of-sample TPH
  over-correction (+16.1) is gone (now +1.4); Pacific Northwest remains the under-corrected exception.
- Crown-width inconsistency quantified: across the 15 western CCF variants, the same species carries
  materially different open-grown crown widths (Douglas-fir spans five curves).

## 3. Roadmap item status

| Item | Status | Commit |
|---|---|---|
| 1. Density-dependent recruitment + held-out re-run | Complete | 8ac8dcd, 514bf6a |
| 2. Four-arm comparison (engine A/B, projector C, figure) | Complete; arm D running | 3fe327b, 411cb27, 8a61a78 |
| 3. Crown-width / MCW unification | Recovery complete (all 25 variants characterized) | 514bf6a |
| 4. Red-team master table (CIs, brms match rate, LS note, volume spec) | Complete | 0ded299, 411cb27 |
| 5. Manuscript (fvs-modern) disturbance-aware + four-arm integration | Complete | fe26e85 |
| Arm D + all-variant four-arm (A/B/C/D in-engine) | Running (jobs 11745408 -> 11745409) | 8b05de5 |

## 4. Key artifacts (on the branch, diagnostics_2026-06-16/ and docs/program/)

- Code: held_out_validation.py (density-dependent recruitment), fourarm_engine.py (engine A/B with CIs),
  fourarm_projector.py (projector A-prime/C with CIs), fourarm_abcd.py (in-engine A/B/C/D), calib_ne.py
  (per-species BAIMULT), mcw_recovery.py (MCW from CCF), brms_match_rate.py, make_fourarm_figure.py.
- Results: held_out_density_dependent_20260618.csv, fourarm_engine_20260618.csv,
  fourarm_projector_NE_20260618.csv, brms_match_rate_20260618.csv, mcw_by_variant_species.csv,
  mcw_cross_variant_spread.csv, fourarm_abcd_20260618.csv (pending the running job).
- Notes (docs/program/): 20260618_held_out_density_dependent_result.md, 20260618_fourarm_result.md,
  20260618_master_results_table.md, 20260618_mcw_recovery_result.md, this file.
- Figure: fourarm_headline_20260618.png (within-framework complementarity).

## 5. Crown-width coverage across all 25 variants

Recovered for the 15 western CCF variants plus ak (power form). ca and ls delegate to R5CRWD and CWCALC
(a short read folds them in). The eastern variants ne, acd, sn, cs zero the base crown-width coefficients
in grinit and so carry no parabolic per-species MCW to unify; Russell and Weiskittel (2010) is the
recommended source if a unified Northeast and Acadian crown width is wanted. Canadian bc and on are out of
scope (no FIA). See 20260618_mcw_recovery_result.md.

## 6. What is proven vs prototype vs blocked

- Proven, out-of-sample: the disturbance artifact; the brms maxSDI and BAIMULT size levers.
- Prototype: density-dependent recruitment (transfers for most variants, not Pacific Northwest).
- Emulation pending in-engine injection: arm C and arm D, currently per-species BAIMULT standing in for
  the fvs-conus growth equations.
- Blocked on one dependency: a true in-engine arm C and D require fvs2py in-process tree loading. Route A
  (build fvs2py injection, maintainer-level) or Route B (multiplier emulation, running now). Recommended:
  submit on Route B with honest labeling, scope Route A as the next maintainer milestone.

## 7. Active jobs and automation

- Job 11745408 (gen_baimult): per-species BAIMULT for 21 variants. Running, near complete.
- Job 11745409 (fourarm_abcd): in-engine A/B/C/D across 22 FIA variants, NSAMP=600. Queued on success of
  the above.
- Scheduled task harvest-armd-fourarm fires 2026-06-18 13:00 EDT: pulls the four-arm result, computes the
  complementarity verdict (does arm D stack the growth and density gains), regenerates the figure as a
  four-arm version, writes 20260618_armD_result.md, commits, and surfaces it.

## 8. Completion strategy

- Phase 1 (now to one week): harvest arm D and refresh the master table and figure; raise thin-sample
  variants (Utah was skipped at n=8) over the fold threshold; extend removal-sim and held-out to larger n.
- Phase 2 (one to three weeks): finish MCW for ca and ls; develop the site-resolved recruitment form so
  Pacific Northwest TPH transfers; volume gross/net sensitivity and biomass via the Fire and Fuels
  Extension.
- Phase 3 (three to six weeks): decide Route A (fvs2py); finish fvs-conus manuscript v2 (Section 3.7,
  Discussion, Figure 1); prune the fvs-modern branch sprawl and merge the integration branch; put
  fvs-conus under version control; reserve the Zenodo deposit as a new version of the fvs_perseus_conus
  concept DOI once validation is out-of-sample-honest.

## 9. Open decisions for the PI

1. Approve or defer Route A (fvs2py in-engine injection). Determines whether the first submission uses the
   emulated or the true four-arm.
2. Submission target and timeline for the fvs-modern combined manuscript and the fvs-conus v2 manuscript.
3. Whether to fold ca and ls into the crown-width unification now or defer with the western set as the
   primary result.

## 10. Commit log this session (branch conus-sf-integration-2026-05-21)

8ac8dcd, 514bf6a, 48a56f2, 3fe327b, 0ded299, 411cb27, 8a61a78, fe26e85, 8b05de5, plus this state document.
