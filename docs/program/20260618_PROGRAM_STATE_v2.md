# FVS modernization program: state of play (v2, end of 2026-06-18)
Single source of truth, superseding the morning 20260618_PROGRAM_STATE.md. Branch
holoros/fvs-modern conus-sf-integration-2026-05-21. Compute on OSC Cardinal (PUOM0008).

## Deployment readiness

- Full all-variant build: clean. build_fvs_libraries.sh builds every variant to a shared library (27 .so).
- Regression suite: 42 passed, 0 failed, 2 skipped (rFVS env only), after fixing a relative-path bug in
  run_regression_tests.sh (it cd's per test, so a relative bin/test dir broke every test after the first).
  Fix: resolve FVS_BIN and TEST_DIR to absolute paths. Deployment-ready.
- Release + Zenodo runbook staged (20260618_release_deposit_runbook.md) with a drafted deposit metadata
  JSON; gated on a human review of the regression result, then merge to main, tag v2026.06.1, and a Zenodo
  new version of concept DOI 10.5281/zenodo.19802673.

## fvs-conus components: forest type + ecoregion (both modes)

Every component now depends on forest type and ecoregion in the species-dependent (unified) and
species-free forms. The unified model family (DG v8c/v9, HG, HCB, crown ratio, HT-DBH, survival) already
carried species + traits + forest type (z_FT) + ecoregion (z_L1/L2/L3). Filled the two gaps this session:
species-free DG (v8_forest_eco) and ingrowth (v3_forest_eco), both compile. Production fits and the paired
held-out ELPD (species-specific vs species-free with forest type + ecoregion) are the next compute step.

## Species-specific vs species-free stress test

Component sigma_sp (response-scale per-species deviation): well identified for mortality (+53% odds),
diameter growth (+23%, mean/sd 10.9), crown recession (+29%); weak for height growth, height-diameter,
HCB. At the stand level (four-arm) the species-specific growth signal gives no median gain, so species-free
is competitive for stand projection while species-specific is retained for mortality, DG, CR.

## Four-arm comparison and the four-way comparison

- Four-arm (default / keyword-calibrated / fvs-conus / combined): engine arms A/B done OOS with CIs;
  headline is QMD median |bias| 15.7 -> 2.2; arm D does not stack over B; true in-engine C/D await Route A.
- Four-way comparison (variants x species x ecoregion x landowner, disturbance-clean): harness running
  (job 11759257), aggregation staged, scheduled harvest 9 PM writes the report and figures.

## Crown-width unification (MCW)

Recovered MCW from CCF coefficients for 15 western variants; one CONUS-consistent consensus curve per
species selected (mcw_conus_consensus.csv). Feed into the crown/CCF competition term is the open step.

## New variants

- ADK (Adirondack): scaffolded from Acadian, FVSadk.so builds and loads. Fixed the add_variant.sh gap
  (it did not create the FVSadk_sourceList.txt build manifest). Calibration runs through the in-process
  .so path (shares Route A) or against NY Adirondack FIA once the variant can be run.
- ACD and AK: to be calibrated against Canadian NFI (MAGPlot) data: Maritimes (NB/NS/PE) for Acadian,
  coastal BC for the AK variant. GATED on staging the MAGPlot data package on Cardinal (the one external
  step; cannot be fetched from this environment).

## Route A (in-engine fvs-conus injection): breakthrough, one detail left

- fvs2py now reads/writes per-tree state (get_tree_attr/set_tree_attr) and loads a stand in memory via
  add_trees (fvsAddTrees). This eliminated the extree segfault.
- A surgical ERR= patch to prtexm.f90 lets the in-process example-tree scratch read fail gracefully, so the
  engine now projects an in-memory stand to a full summary.
- Remaining: the added trees register by count (TPA) but their size (DBH/HT) does not reach the BA
  computation; fvsAddTrees does not run the full per-tree inventory initialization the DB/treeinit path
  does. Focused FVS tree-init source work: populate the derived inventory arrays after the add, or read/set
  sizes at the cycle stage where the arrays are live (not stop point 7). Once done, the fvs-conus DG/HD
  predictions wire into the stop-point-5 loop for the true in-engine arms C and D.

## Open decisions for the PI

1. Stage the MAGPlot data on Cardinal to unblock ACD and AK calibration (one wget; see the calibration
   plan doc).
2. Review the four-way comparison report (posts tonight), then approve the release: merge to main, tag
   v2026.06.1, Zenodo new version.
3. Schedule the Route A tree-init finish as a focused FVS-source session.

## Commit log (2026-06-18, this program)
8ac8dcd density-dependent recruitment; 514bf6a OOS result + MCW recovery; 48a56f2 four-arm harness;
3fe327b/0ded299/411cb27 four-arm + master table; 8a61a78 figure; fe26e85 manuscript; 8b05de5 arm-D harness;
eae8b4f arm-D result; c32d40c/60960d6/9680389/4a2fa16 Route A diagnosis; 1b51a2b MCW consensus;
9e4eae7/bd10a12 forest-type+ecoregion; c7b324f build+regression deployment-ready; 5d1acc6 Canadian-NFI plan;
72406d0/1aba13c ADK; 4266a47/3415706/9f84242 Route A in-process projection works, tree-init gap documented.

## 2026-06-18 evening update

- Four-way comparison COMPLETE and committed (39dcf8a, 8c1764c): calib_4way_margins.csv plus
  fig_4way_landowner.png and fig_4way_ecoregion.png. Variant margin is full (951 conditions, n 40-112/variant);
  calibration cuts QMD bias across 11 of 13 variants (kt 26.3->3.9, ec 14.7->-2.2, sn 14.9->5.7, ne 8.7->2.9,
  acd 9.9->5.0). Ecoregion (L2) and landowner margins reported on the 123-condition disturbance-clean joined
  subset; species margin is thin there and is covered by the species-stress test. Join gate lowered 0.2 -> 0.1
  (the cspiv6 disturbance-clean universe legitimately carries ~13% of harness conditions).
- MAGPlot staged and fully ingested (c457b61). Confirmed real schema. ACD target = NB (NS/PE not in package);
  AK target = BC. Built remeasurement pairs from source: ACD/NB 263 (prior, validated), AK/BC 2,451 clean
  protocol-consistent pairs (from 6,845 raw; filtered for matching subplots/tag-limit across visits).
- ACD/NB is the validated maritime reference: near-unbiased BA (-0.04%, R2 0.88), QMD +9% is the calibration
  target. AK/BC machinery and species crosswalk validated (100% coverage); FVS-AK projection blocked on a
  treeinit TPA-scaling + cycle-length ingestion bug under region 10 (reads ~7x low TPA, 1-year cycles). One
  focused engine-debug fix from a real AK bias number. Route A tree-init work and this AK ingestion fix are the
  same class of focused FVS-source task.
