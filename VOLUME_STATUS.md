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

## 2026-07-17 (pm): tree-level decomposition + responsible layer (DEFINITIVE)

### Method
Isolated the estimator from projection with a cycle-0 identity test: sample large
live ME FIA trees (known DBH, HT, SPCD, VOLBFNET, VOLCSNET, VOLCFNET), build
single-tree NE stands, run FVS one cycle (no growth), read per-tree NSVB board
foot VOL(10) and merch cubic from FVS_Summary2, compare to FIA's own NSVB values
at IDENTICAL inputs. Harness: rescore/pertree_estimator.py, pertree_ctype.py,
trace_run.py; evidence in rescore/FINDINGS_pertree.txt and maxlen_sweep_results.txt.

### Tree-level ratio vs DBH (FVS/FIA, identical inputs, default R9 mrules MAXLEN=8ft)
SW  BF_ratio 1.05 (10-13") -> 1.35 (23-40");  CF_ratio 0.99 -> 1.14
HW  BF_ratio 0.94 (10-13") -> 1.42 (23-40");  CF_ratio 1.11 -> 1.28
DBH>=16 overall:  SW BdFt +28.4% / cubic +11.0% ;  HW BdFt +31.8% / cubic +20.9%
Board-foot-per-cubic (SW 23-40"): FVS 6.43 vs FIA 5.43 bf/cf (+18%).
=> divergence is board-foot-specific, size-increasing, present at the TREE level
   with zero projection. Small SW cubic matches FIA (0.99) -> same NSVB cubic
   estimator confirmed; the split is entirely in the board-foot derivation.

### Which layer (each candidate tested by A/B on identical-input binaries)
(a) CTYPE convention (FIA 'I' vs FVS 'F' branch in nsvb.f90): INERT. VOL(10) is
    byte-identical F vs I (SW26" 835 bf both). Ruled out.
(b) Merch top MTOPP (7.6 SW / 9.6 HW): conservative, NOT the driver. A lower FIA
    top (7.0/9.0) would lengthen the sawlog and make FVS bigger, wrong direction.
(c) Log segmentation length (Region-9 mrules MAXLEN). THIS IS THE LAYER.
    MAXLEN sweep, DBH>=16, identical inputs (cubic pinned at SW 1.110 / HW 1.209):
        MAXLEN(ft)   SW BdFt ratio   HW BdFt ratio
           8            1.284           1.318   (production default)
          16            1.247           1.284   (standard International log)
          32            1.190           1.196
          64            1.023           0.960   (whole-bole, no bucking)
    At whole-bole the same-tree board-foot excess collapses to ~0 (SW +2%, HW -4%)
    while cubic is untouched. So the ENTIRE same-tree board-foot-specific excess is
    the 8-ft bucking feeding NVEL INTL14: each short log gets a fresh 0.5"/4ft
    internal-taper credit off a real (flatter-at-base) small-end diameter, over-
    crediting board feet, and the effect grows with diameter/log count.

### How much each layer explains (large SW, board-foot terms)
- Scorecard aggregate large-tercile SW BdFt bias was +75%. Of that:
  ~ +28 pts = same-tree 8-ft-segmentation vs FIA whole-bole board-foot method;
  the shared cubic component (~+11 pts) is folded in; the remaining ~+40 pts is
  PROJECTION (FVS grows large trees larger than observed; the board rule's
  nonlinearity amplifies size error far more than cubic). Cubic same-tree is only
  +11% (SW), so most of the cubic scorecard excess is also projection.
- Standard 16-ft logs would trim only ~4 pts of the same-tree excess; even 16-ft
  leaves +25%. FIA's VOLBFNET matches ~whole-bole, NOT 8-ft or 16-ft bucking.

### Determination: NO board-foot CODE BUG (inherent method difference)
NVEL INTL14 and the nsvb.f90 accumulation are correct: no units error, no double
count, correct International 1/4 slot (VOL(10)), correct sawlog domain. The
large-tree board-foot divergence is inherent to two different board-foot METHODS:
FVS computes operational International 1/4 by bucking the sawlog into 8-ft logs (a
legitimate Region-9 merch standard, MDL=CLK/NVB), while FIA's NSVB VOLBFNET is a
whole-bole International integration of the same taper (no short-log bucking).
Per the task's own guidance ("if the divergence is inherent ... do NOT force a
change"), no production edit was made. Throwaway trace + env-gated MAXLEN/MTOP/
CTYPE binaries (bin-r9vol-trace, bin-r9vol-seg) were used only to localize; the
source edits were reverted and the worktree source is pristine (git clean).

### Actionable next step (scorecard fairness, a REPORTING choice not a growth fix)
To make the scored BdFt like-for-like against FIA VOLBFNET, compute the scored
board foot with FIA's whole-bole International method (demonstrated ML64 run:
large-tree same-tree BdFt bias -> ~0, cubic unchanged). Do NOT hard-set MAXLEN to
a non-physical value in the operational merch path; instead either (i) score BdFt
via the whole-bole integral in the summary path when comparing to FIA, or (ii)
document the 8-ft-log convention so the residual board-foot column is read as
operational, not FIA-matched. The remaining, larger board-foot excess after this
alignment is the growth-model projection signal (large trees grown too big),
which is where board-foot accuracy should next be pursued.

---

## FIA whole-bole International 1/4-in board foot: env-gated scored option (2026-07-17)

Branch feat/r9-intl-boardfoot-20260717. Flag FVS_BF_FIA_WHOLEBOLE (env-gated,
default operational 8-ft when unset).

### Which routine actually carries the SCORED board foot (empirically settled)
The scored FVS_Summary2 BdFt is NVEL VOL(10) from the Region-9 Clark path:
fvsvol.f90 BFVOL entry -> BBFV = TVOL(10) (METHB=9) -> VOLINITNVB/VOLINIT ->
R9CLARK -> subroutine r9bdft (r9clark.f90). r9bdft is the routine that builds
the board foot: it sums the closed-form International 1/4-in polynomial over
operational Region-9 8-ft logs. The earlier premise that scored BdFt came from
nsvb.f90 NVB_CalcLOGVOL VOL(10) was WRONG for the board foot: that nsvb VOL(10)
is the CUBIC NSVB call and is INERT to the scored board number. Settled by a
decisive same-binary A/B on the flag (below): flipping FVS_BF_FIA_WHOLEBOLE,
which gates ONLY the VOL(10) accumulation inside r9bdft, moves the scored BdFt.
If r9bdft were off the scored path the number could not move.

### Decisive identity test (cycle-0, 200 large ME FIA trees, same wholebole binary)
Same binary bin-r9vol-wholebole/FVSne; only the env flag differs. FVS/FIA ratio
of the scored BdFt vs FIA VOLBFNET, and MCuFt vs VOLCFNET (cubic control).

| group        |  n  | BdFt ratio op-8ft | BdFt ratio whole-bole | CF ratio (both) |
|--------------|-----|-------------------|-----------------------|-----------------|
| SW DBH>=16   | 60  | 1.284             | 0.967                 | 1.110           |
| HW DBH>=16   | 60  | 1.318             | 0.958                 | 1.209           |
| ALL DBH>=16  | 120 | 1.301             | 0.963                 | 1.159           |

The scored BdFt MOVES (op-8ft +30% -> whole-bole -3.7% same-tree) => r9bdft IS
on the scored path. With the flag on the scored board foot lands within ~3-4% of
FIA VOLBFNET (the tightest a single 0.5in/4ft International integration reaches
vs FIA's exact NSVB-taper integration). CF ratio byte-identical both legs: the
change is board-foot-only, cubic untouched.

### Implementation
r9clark.f90 subroutine r9bdft: env flag FVS_BF_FIA_WHOLEBOLE (GETENV, cached
SAVE). When =1 the per-8-ft-log VOL(10) accumulation is suppressed and VOL(10) is
set from ONE whole-bole International polynomial evaluation (length = summed
sawlog lengths + (NUMSEG-1)*TRIM inter-log trims; small-end DIB = merch-top log
DIB), single round to nearest 5 bf. TRIM passed in via a new r9bdft argument.
VOL(2) Scribner, cubic VOL(4)/VOL(7), and the per-log LOGVOL array UNCHANGED.
Default (flag unset) preserves the operational Region-9 8-ft-log board foot
bit-for-bit. Worktree scorecard harness copy
rescore/run_5model_scorecard_cohort_vol.py points ne at bin-r9vol-wholebole and
now os.environ.setdefault('FVS_BF_FIA_WHOLEBOLE','1') so future scored runs use
FIA parity by default; the A/B drivers pop the var for their op-8ft leg.

### Cohort scorecard, BdFt aggregate bias % by QMD tercile (500 pairs, variant ne)
BEFORE = operational 8-ft logs, AFTER = FIA whole-bole. MCuFt identical
before/after in every cell (cubic unchanged). QMD tercile edges (fvs_base
BEFORE): T1 <= 9.54 < T2 <= 11.39 < T3.

| model         | tercile   |  n | BdFt before% | BdFt after% | dBdFt | MCuFt (bef=aft)% |
|---------------|-----------|----|--------------|-------------|-------|------------------|
| fvs_base      | T1_small  | 77 | 32.1         | 37.6        | +5.5  | 29.5             |
| fvs_base      | T2_mid    | 96 | 39.2         | 40.0        | +0.8  | 29.1             |
| fvs_base      | T3_large  | 89 | 54.1         | 39.2        | -14.9 | 41.2             |
| fvs_base      | T3 SW     | 13 | 75.2         | 56.6        | -18.6 | 29.4             |
| fvs_base      | T3 HW     | 76 | 50.6         | 36.3        | -14.3 | 43.4             |
| fvs_regional  | T3_large  | 90 | 50.0         | 35.9        | -14.1 | 35.8             |
| fvs_regional  | T3 SW     | 13 | 69.4         | 50.9        | -18.5 | 23.5             |
| fvs_regional  | T3 HW     | 77 | 46.8         | 33.4        | -13.4 | 38.1             |
| organon       | T3_large  | 98 | 39.6         | 25.4        | -14.3 | 24.5             |
| organon       | T3 SW     | 14 | 63.8         | 47.8        | -16.0 | 18.8             |
| organon       | T3 HW     | 84 | 35.5         | 21.6        | -14.0 | 25.6             |
| conus_spdep   | T3_large  | 97 | 40.9         | 27.1        | -13.8 | 25.2             |
| conus_spdep   | T3 SW     | 14 | 62.6         | 45.9        | -16.7 | 17.6             |
| conus_spdep   | T3 HW     | 83 | 37.3         | 23.9        | -13.3 | 26.6             |
| conus_climate | T3_large  | 98 | 39.1         | 24.7        | -14.4 | 24.2             |
| conus_climate | T3 SW     | 14 | 59.8         | 42.0        | -17.8 | 16.3             |
| conus_climate | T3 HW     | 84 | 35.7         | 21.8        | -13.9 | 25.7             |

(a) Same-tree method artifact removed: large-tercile BdFt drops 14-15 pts per
    model and now tracks the cubic level (AFTER BdFt_T3 ~= MCuFt_T3 every model).
(b) Cubic unchanged: MCuFt identical before/after in all cells.
(c) Residual large-tree BdFt bias = TRUE PROJECTION SIGNAL (large trees grown too
    big; the board rule no longer amplifies it): fvs_base +39.2%, fvs_regional
    +35.9%, organon +25.4%, conus_spdep +27.1%, conus_climate +24.7% (T3_large,
    matching cubic +24-41%). Small trees (T1) rise ~4-6 pts under whole-bole
    (short boles, minor absolute board foot); the large tercile is where the fix
    acts.

### Caveat
Whole-bole International sits ~3-4% BELOW FIA VOLBFNET at large DBH (0.963 ALL),
the inherent residual between a single 0.5in/4ft-taper International integration
and FIA's exact NSVB-taper integration; the tightest the International rule
reaches, matching the established MAXLEN=64 evidence. Default (flag unset)
preserves the operational Region-9 8-ft-log board foot bit-for-bit.
