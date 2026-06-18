# Full build and regression test, all variants (deployment readiness)
2026-06-18. Branch conus-sf-integration-2026-05-21. Built and tested on OSC Cardinal (gcc 12.3.0).

## Result: deployment-ready

The full all-variant build is clean and the regression suite passes with no failures once a path bug in
the test harness is fixed.

| stage | result |
|---|---|
| Build (build_fvs_libraries.sh src-converted ./lib) | 27 shared libraries (.so) built, all variants, no errors |
| Standalone executables | present and valid (ELF 64-bit) for every variant |
| Regression suite (run_regression_tests.sh, absolute paths) | 42 passed, 0 failed, 2 skipped |

The skipped two are the rFVS library-load tests (rFVS is not installed in the batch environment); they are
an environment dependency, not a code failure. Every standalone simulation test passed, including exact-match
and database-output checks across cr, cs, ie, ne, sn and the others.

## Harness bug found and fixed

The first regression run reported 23 variants as "no executable" and only ran FVSak (which then failed). The
cause was not the build: every executable was present and valid. The harness took a relative FVS_BIN
(./lib) and TEST_DIR, and Part 2 cd's into a per-test working directory for each simulation without
returning, so on every test after the first the relative ./lib/FVSx path resolved against the wrong
directory and was reported missing. Passing absolute paths makes all variants run; the fix in
deployment/scripts/run_regression_tests.sh resolves FVS_BIN and TEST_DIR to absolute paths right after
argument validation, so the harness is correct regardless of how it is invoked. With the fix: 42/42.

## Deployment readiness

The build produces the .so libraries for the ctypes/fvs2py and microFVS deployment paths plus standalone
executables for every variant, and the regression suite confirms numerical correctness against saved
baselines. This is the deployable state for the new release. Next: the release-and-deposit runbook
(20260618_release_deposit_runbook.md) takes it to a tagged v2026.06.1 and a Zenodo new version of the
fvs-modern concept DOI.

Jobs: 11758564 (build + first regression, exposed the path bug), 11759089 (regression re-run, surfaced the
relative-path cause), 11759xxx (regression with absolute paths, 42/42).
