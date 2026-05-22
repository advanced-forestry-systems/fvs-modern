# FVS-ACD session final handoff

2026-05-22. Wrap-up of the FVS-ACD divergence, calibration, 12.3.6 adoption, and
Canadian MAGPlot work. All artifacts live in `calibration/acd_divergence_2026-05-21/`
on `holoros/fvs-modern` main; the same set is mirrored in the selected fvs-online
folder.

## What is done and merged

- PR #25: divergence diagnosis, annualized calibration framework, decomposition
  harness, real-model attribution, fitted calibration table, bridge patch + helper.
- PR #26: AcadianGY 12.3.6 validation (equivalence off == 12.3.5, direction on).
- PR #27: FIA vs Canadian MAGPlot cross-validation matrix (#128).
- PR #28: comprehensive 8/8 test pass; tag-fixed 12.3.6 deployed as the working
  `fvsOL/inst/extdata/AcadianGY.R` (12.3.5 backed up).
- PR #29: wiring the annual diameter calibration into the model worsens stand BA.
- This handoff + the ingrowth experiment (final commit).

## Key findings

1. AcadianGY.R and FVS-ACD diverge in the customRun bridge, not the equations
   (the fork's AcadianGY.R is byte-identical to the 12.3.5 reference). Drivers:
   FVS-NE multipliers injected on the wrong model/scale, HT/CR dubbing disabled,
   CSI source mismatch, option/ingrowth differences. Annual vs periodic is not a
   driver (the bridge already steps annually).

2. 12.3.6 = canonical 12.3.5 + an in-source opt-in size-dependent mortality
   correction (`ops$MORTCAL`), survivors-only so ingrowth is preserved. Off by
   default and byte-for-byte equivalent to 12.3.5 when off (verified on synthetic
   and real NB data). Version tag corrected to AcadianV12.3.6.

3. MORTCAL is FIA-specific (#128). Maine FIA: off +15.4%, on +8.6% (it halves a
   real over-projection). Canadian NB MAGPlot: off -0.04% (already unbiased), on
   -6.51% (it over-corrects). Keep MORTCAL off for Canada/MAGPlot.

4. The annualized diameter calibration is a TREE-LEVEL objective and must not be
   used to fix stand BA. Applied alone it worsens the FIA over-projection
   (+15.4% to +20.2%); with MORTCAL it is +13.3%, still worse than MORTCAL alone
   (+8.6%). The model relies on compensating errors; calibrating one component
   breaks the balance.

5. Ingrowth experiment (this round): toggling INGROWTH on vs off produced
   IDENTICAL FIA results (+10.9%, QMD 5.15 vs observed 4.97). The model recruits
   essentially zero new trees over the 5 to 10 year interval in this harness, so
   the residual over-projection is a mean-size (QMD) effect: the model's trees
   run about 3.6 percent larger than observed because observed stands gain small
   ingrowth the model does not produce. This is the concrete form of #127.

6. Canadian MAGPlot to Fortran FVS is unblocked. The route is the standalone
   FVSne/FVSacd binary + DATABASE keyfile + SQLite, one stand per subprocess
   (`magplot_fvs_runner.py`). The fvs2py route was the blocker (needs Python
   >= 3.11 and its run() never completed the DATABASE run). FVS-NE and FVS-ACD
   produce distinct projections on NB data.

## Remaining work, prioritized

1. Recruitment / ingrowth submodel (the #1 lever). The residual after MORTCAL is
   driven by the model not producing the observed small-tree ingrowth (finding
   5). Make AcadianGY recruit realistic small trees over a remeasurement interval
   and re-check QMD and BA. This is the genuine stand-level fix; diameter and
   mortality multipliers cannot close it.
2. If both tree-level accuracy and stand BA are wanted, calibrate jointly at the
   stand level so diameter, mortality, and ingrowth corrections are mutually
   consistent, not independent per-component multipliers.
3. Production polish: make the standalone FVS summary horizon deterministic;
   per-stand site index/elevation for the MAGPlot runner; stabilize the Cardinal
   upstream_fvs_check checkout (it keeps vanishing); propagate the version-tag fix
   to the canonical 12.3.6 copies (HRF, seven_islands).

## Reproduce (Cardinal)

    module load gcc/12.3.0 R/4.4.0
    cd ~/acadgy_fia_verify
    Rscript cardinal_acadgy_insource_v16.R     # FIA: canonical vs MORTCAL
    Rscript cardinal_acadgy_calib_v17.R        # FIA: + diameter calibration
    Rscript cardinal_acadgy_ingrowth_v18.R     # FIA: ingrowth on vs off
    cd ~/magplot_verify && Rscript cardinal_magplot_insource_v16.R  # NB cross-val
    # Fortran NB:  module load python/3.12; python3 ~/magplot_fvs_runner.py
