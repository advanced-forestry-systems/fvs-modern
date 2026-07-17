# NE board-foot volume status (feat/r9-intl-boardfoot-20260717)

Last updated 2026-07-17. Worktree: /fs/scratch/PUOM0008/crsfaaron/fvs-modern-volume

## Bottom line
The committed R9 International 1/4-in board-foot change (commit 29541ef) is
CORRECT-BUT-INERT for the cohort scorecard. It is latent infrastructure for the
legacy (LFIANVB=.FALSE.) board-foot path and can stay committed. It does NOT
touch the board foot the scorecard actually scores. The real large-tree
board-foot over-prediction lives in the NSVB estimator, not r9vol/r9clark/METHB.

## LFIANVB finding (definitive, with evidence)
LFIANVB = .TRUE. at runtime for every NE scorecard stand.
- Set unconditionally at initialization: src-converted/ne/grinit.f90:163
  `LFIANVB = .TRUE.  ! 2026-05-16 holoros fork: NSVB (Westfall et al. 2024) default`.
  It is NOT gated on the FIAVBC keyword.
- The scorecard harness KEYFILE (run_5model_scorecard_cohort_vol.py) emits only
  DATABASE/SUMMARY 2 + calibration keywords; it never emits FIAVBC. That is
  irrelevant because grinit already forces LFIANVB=.TRUE.
- The earlier suspicion that LFIANVB was FALSE (because FIAVBC is never invoked)
  is REFUTED. The earlier handoff claim that the binary uses LFIANVB=.TRUE. is
  CORRECT.

## Which path feeds the scored BdFt (definitive, with evidence)
The scored BdFt column = last-cycle `BdFt` from FVS_Summary2 (harness lines
458 / 853). FVS_Summary2 is written by dbssumry2 from the IOSUM/ISUMARY summary.
Under LFIANVB=.TRUE. the per-tree board foot is produced by the NSVB path:
  fvsvol.f90 (LFIANVB branch) -> VOLINITNVB -> nsvb.f90 NVB_CalcLOGVOL, which
  accumulates board foot log-by-log from an NSVB taper profile:
    VOL(2)  += LOGVOL(1,I)  via CALL SCRIB   (Scribner Decimal C)
    VOL(10) += BFINT        via CALL INTL14  (International 1/4-in)
  (nsvb.f90 ~lines 1062-1079).
The scorecard's BdFt is the International 1/4 slot (NSVB national reporting
convention; FIA VOLBFNET is also International, so the comparison is
like-for-like International vs International).

The METHB->TVOL(10)/TVOL(2) selection in fvsvol.f90 (lines 519-524) and the
r9clark/r9vol (900CLKE / 90xDVEE) equations are NOT on this path. The classic
SPCBV/BFV summary (vls/vols.f90) is bypassed for FVS_Summary2 under NSVB.

## Proof that METHB / r9vol are inert (A/B, genuinely distinct binaries)
Three NE binaries built from identical source except the METHB default:
- bin-r9vol/FVSne          METHB=9 (International), + r9vol INTL14 edit. md5 37f7457...
- bin-r9vol-base-true/FVSne METHB=6 (Scribner),     + r9vol INTL14 edit. md5 46e8c3a...
  (built by reverting grinit METHB 9->6, compiling, then git-checkout restore;
   grinit is back to METHB=9, git clean.)
The two binaries have DIFFERENT md5 (genuinely distinct in the METHB path), yet
the 500-pair cohort scorecards are BYTE-IDENTICAL (md5 0103bbe... for both
before_scribner and after_intl). BdFt aggregate bias is unchanged in every QMD
tercile:
    tercile     BdFt (Scribner METHB=6)   BdFt (Intl METHB=9)
    ALL T1               +32.1%                 +32.1%
    ALL T2               +39.2%                 +39.2%
    ALL T3 (large)       +54.1%                 +54.1%
    SW  T3 (large)       +75.2%                 +75.2%
    HW  T3 (large)       +50.6%                 +50.6%
Cubic (MCuFt) likewise unchanged. Conclusion: the committed change has zero
effect on the scored board foot. (The prior before_scribner run was also invalid
for a second reason: its "baseline" binary was rebuilt from the same METHB=9
source and was never actually METHB=6.)

## What the committed r9vol fix does
- ne/grinit.f90: METHB default 6 -> 9 (selects TVOL(10) in fvsvol.f90). Live only
  on the legacy non-NSVB path; inert here.
- volume/NVEL/r9vol.f90: replaces the dead VOL(10)=VOL(2) Scribner alias in the
  90xDVEE (Gevorkiantz, METHB==5) path with a genuine INTL14 log accumulation.
  NE never uses the DVEE path (uses 900CLKE Clark), so this is dead code for NE.
Both are harmless and correct for the paths they touch. Leave committed.

## The real over-prediction and the corrected next step
The large-tree board-foot over-prediction is board-foot-SPECIFIC and lives in the
NSVB estimator (nsvb.f90). Evidence: for large softwood (SW T3) cubic over-
predicts only +29% but board foot +75%; the board-foot equation amplifies large-
diameter volume far more than cubic. The driver is the NSVB taper profile
(NVB_CalcDiaAtHT, coefficients a,b from tables4.inc feeding TCUFT) and its 16-ft
log segmentation into INTL14, not the board rule.

Scoped plan (larger job, do NOT bolt onto this branch blindly):
1. Confirm the exact IOSUM/FVS_Summary2 board-foot slot (VOL(10) vs VOL(2)) and
   the tree->summary wiring under LFIANVB (add a one-shot trace in a throwaway
   build; ~5 min). Expected VOL(10) International.
2. Extract per-tree NSVB VOL(10) vs FIA VOLBFNET for matched large-diameter
   trees (DBH tercile T3, split SW/HW) and quantify the per-tree bias vs DBH/HT.
3. Localize: (a) merch spec mismatch (BFTOPD / BFMIND / stump vs FIA), (b) the
   NSVB taper DIB profile at large DBH, or (c) the cubic (TCUFT) that seeds the
   segmentation. Fix at the smallest responsible layer; refit/adjust only that.
4. Build a test binary, A/B rescore with rescore/run_rescore_base_true.py pattern
   (swap the binary path), confirm SW/HW T3 board-foot bias shrinks and cubic is
   unchanged.

Turnkey stub: rescore/run_rescore_base_true.py (rescore one binary into a tag) and
rescore/terciles.py (tercile A/B) are the reusable harness; point them at the new
test binary vs bin-r9vol.
