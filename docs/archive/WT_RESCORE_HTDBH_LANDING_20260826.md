# Landing wt-rescore and wt-htdbh onto a review branch (2026-08-26)

Branch: `fix/land-wt-rescore-htdbh-20260826`, based on `release/v1-code-20260811` tip
(c40e2af, "Merge save/worktree-20260725 into release/v1-code-20260811"). This branch
is a review artifact only. Nothing has been merged into `release/v1-code-20260811`,
`main`, or `v1.0.0`, and none of those refs were moved.

## Landed: native HT-DBH (height-diameter) hook

Source: worktree `wt-htdbh`, branch `feat/conus-htdbh-native`, HEAD `7ad08e0`
("Merge branch 'feat/htdbh-serializer' into feat/conus-htdbh-native"), six commits
between the release merge-base and HEAD. Prior audit notes described this as "the
correct default."

What it does: adds `GREGHTDBH`, a native Fortran HT-DBH evaluator supporting six
functional forms (`src-converted/base/greghtdbh.f90`, common block
`src-converted/common/GREGHDMC.f90`), loaded once per run by `GREGLOADHD` and wired
into the NE variant's `HTDBH` MODE 0 dispatch behind the `FVS_GREGHTDBH` environment
toggle. The toggle is unset by default, so this hook changes no existing run's
behavior unless explicitly enabled. Coefficients come from a CONUS-wide Marshall
height-diameter refit, serialized per FIA species code and crosswalked to each of
25 FVS variants through the existing `FIAJSP` species mapping
(`sf_integration_dev/htdbh_native/62m_htdbh_to_variant.R`): 412 species covered at a
mean 96.1% per-variant coverage, with species lacking an acceptable CONUS fit
falling back to native NE-TWIGS HTDBH.

Evidence of completeness, from the original commit messages: a real(8) Fortran
evaluator was checked against an R `model.formula` reference across 412 species by
DBH 1 to 50, with an overall max absolute difference of 5.1e-13 against a 1e-6
threshold, and FVSne built clean with the hook compiled in.

Conflict resolution (commit `4e85db0`): one trivial conflict in
`src-converted/bin/FVSne_sourceList.txt`. Both sides added a different new source
filename at the same line position, release's side already carrying
`../base/greghghg.f90` (an unrelated Greg four-arm height-growth and crown-recession
hook, added independently, confirmed present and distinct by reading its header and
`git log`) and the incoming side adding `../base/greghtdbh.f90`. Verified both files
are real, distinct, already-present source files with no overlapping content, so
both lines were kept in the source list rather than either being dropped.

Also landed (commit `03cfc2c`): two files that were sitting uncommitted in the
`wt-htdbh` worktree on top of its HEAD, judged to be genuine fix content rather than
build output: `calibration/R/30b_build_htdbh_input.R` (the input-prep script feeding
the `62m_htdbh_to_variant.R` serializer already landed) and
`calibration/data/conus_htdbh_build_summary.csv` (the per-species and per-variant
coverage summary the commit messages cite: 412 of 438 species, mean 96.1%
coverage). Excluded: `bin-htdbh/FVSne`, a compiled binary with no source-control
value.

## Not landed: DG size-deceleration ceiling (from wt-rescore)

`wt-rescore` (branch `eval/rescore-fixed`, HEAD `1380970`, described in its own
commit message as "ready to ship, PR #101") carries two additional fixes on top of
the same HT-DBH ancestry, merged in via `origin/feat/dg-size-ceiling` and
`origin/feat/gompit-smalldbh-guard`. Neither is landed on this branch.

The DG ceiling (commit `a217eae`, "DG native hook: add size-deceleration ceiling to
fix long-horizon runaway") adds a per-species `GDGMAX(MAXSP)` maximum-diameter field
to the `GREGMC.f90` common block and a logistic `DECEL` term in `GREGDGV`
(`src-converted/base/gregdghg.f90`) that ramps the annual DG increment to zero as
DBH approaches 85 to 100 percent of `GDGMAX`, ported identically from the PR #90
native-hook track. Its own verification claims are concrete: on 20 NE stands over
300 years, runaway stands (topHT > 80 m, QMD > 150 cm, or BA > 120 m^2/ha) at year
300 dropped from 15 of 20 to 0 of 20; year-300 QMD with the ceiling on had a median
of 64.5 cm and a p90 of 73.2 cm versus a p90 of 141.6 cm without it; short-horizon
DBH-increment percent RMSE moved only 174.97 to 175.95 percent.

This is the **same underlying defect and the same mechanism name** as the parallel
session task tracked as "Fix DG runaway: complete GDGMAX patch + disable default"
(branch `fix/gregdg-deceleration-20260826`), which as of this check had an empty
worktree (zero commits, zero diff against `release/v1-code-20260811`). Attempting to
land `a217eae` here surfaced a **real, substantive conflict**, not just a
duplicate-effort concern: `release/v1-code-20260811`'s current tip already carries a
`GDGMAX(MAXSP)` field in `GREGMC.f90`'s `/GREGMR/` common block (added by the bulk
`495541e` "Save Cardinal worktree" commit, merged into release via
`save/worktree-20260725`), but paired with a different, PC-form climate variable set
(`GDD0, GTD, GPPT_SM, GDD18, GPC1, GPC2`) rather than wt-rescore's simpler
`GEMT, GTD, GELEV` set, and with no corresponding `DECEL` logic yet wired into
`GREGDGV` in `gregdghg.f90` at all (grep for `GDGMAX`/`DECEL`/`DBHMAX` there returns
nothing). In other words, release already has an incomplete, differently-shaped
`GDGMAX` scaffold sitting unused, which matches the word "complete" in the parallel
task's own title. Cherry-picking `a217eae` produces genuine content conflicts in
both `GREGMC.f90` (conflicting common-block layouts) and `gregdghg.f90` (conflicting
variable declarations in `GREGDGV`), not a mechanical rename-style collision.
Resolving this safely requires deciding whether the review branch adopts
wt-rescore's simpler climate-variable layout or re-derives the deceleration logic
against release's already-merged PC-form layout, which is a design decision, not a
merge-conflict resolution, and belongs to whoever is completing the parallel GDGMAX
task rather than to this landing pass. **Recommendation: point the
`fix/gregdg-deceleration-20260826` work at wt-rescore's `a217eae` as a working
reference implementation with real verification numbers, but re-derive it against
release's current `GREGMC.f90` layout rather than cherry-picking as-is.**

## Not landed: gompit small-DBH survival guard (from wt-rescore)

Commits `58006dd` ("feat(gompit): small-DBH survival guard to stop young-cohort
collapse") and `a56b453` ("chore(gompit): log guard params at activation; add
verification harness") add a `DBHV` argument to `GOMPSURV`
(`src-converted/base/gompmort.f90`, called from all 22 variant mortality drivers)
and floor per-cycle survival only for trees with `0 <= DBHV < GDBHMIN`, leaving
established trees (`DBHV >= GDBHMIN`) on the byte-identical pre-guard hazard path
(verified 0 of 1380 divergences across 92 species by DBH by crown-ratio
combinations in the commit's own harness). The stated defect: raw gompit annual
survival can collapse toward 0 for some NE and SN conifers at the crown-ratio and
crown-closure combination of a young dense seedling cohort, killing close to 100
percent of the stand in the first 10-year cycle, falsifying the Bakuzis assessment.

Attempting to cherry-pick this onto the review branch (after HT-DBH was already
landed) produced conflicts in `src-converted/base/gompmort.f90`,
`src-converted/common/GOMPMC.f90`, and `src-converted/vls/morts.f90`. Inspection
shows this is another real design fork against content already merged into
`release/v1-code-20260811`, from the same `495541e` bulk commit: release's current
`GOMPSURV` already has its own fix for a related problem, a **universal** survival
floor (`GSFLOOR`, paired with a `GMORTMULT` hazard multiplier) that clamps every
tree's per-cycle survival at `GSFLOOR**FINTL` whenever `GSFLOOR > 0`, explicitly
documented on release as substituting for the FVS `MORTMULT` keyword "which does not
reach GOMPIT." Release's `GOMPMC.f90` has no `GDBHMIN`-equivalent at all, so its
floor is coarser by design (every tree, always) where wt-rescore's is targeted
(small trees only, established trees untouched). These are two different design
philosophies for the same subroutine's failure mode, not a mechanical collision, so
resolving it here would mean silently picking a mortality-floor design for Aaron
rather than presenting the choice. Left unresolved and unlanded.
**Recommendation: this needs an explicit decision from Aaron, universal GSFLOOR
floor already on release versus wt-rescore's targeted small-DBH guard, or whether
the two are meant to compose (e.g., GSFLOOR as an operator-set global backstop with
wt-rescore's GDBHMIN gate layered underneath it for the seedling-collapse case
specifically).**

## Conflict-check discipline

A throwaway worktree was created at `release/v1-code-20260811` tip (later found to
be dirty from a prior stale worktree registration reused by `git worktree add -f`;
that worktree was fully removed and a fresh one created and verified clean, 0 files
in `git status`, HEAD matching the release tip exactly) and used to dry-run
`git merge --no-commit --no-ff` and `git cherry-pick -n` for each candidate commit
before touching the real landing branch. All conflict markers found were read in
full on both sides and their content verified (file existence, distinctness, and in
the two "not landed" cases, the actual competing common-block layouts and
subroutine logic already on release) before any resolution or landing decision was
made. No conflict marker was deleted without confirming what both sides actually
contained.

## Build artifacts excluded

`bin-rescore/` (compiled `FVSne`, `FVSpn`, `FVSsn` binaries in the wt-rescore
worktree) and `bin-htdbh/FVSne` (compiled binary in the wt-htdbh worktree) were not
landed. Both are build output, not source, and neither worktree's git history
tracks them.

## Cardinal MCP connectivity note

Cardinal MCP connection reliability during this session was substantially worse
than the documented roughly 40-45 percent interactive-call timeout rate, with
stretches of 6 to 9 consecutive `Connection closed` failures on both
`cardinal_test_connection` and `cardinal_run` before recovering. All work was
completed by reissuing on timeout per standing convention; no step was skipped or
guessed around a blocked call.
