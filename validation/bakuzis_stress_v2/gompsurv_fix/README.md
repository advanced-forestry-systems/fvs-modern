# GOMPSURV sixth-argument fix (October 28, 2026)

`base/gompmort.f90` has carried a six-argument `GOMPSURV(ISPC, CR, CCHV, FINTL, SURV, DBHV)` since commit 495541e (July 25, 2026), when the size-aware `b5 * log(DBH)` term was added. Only `vls/morts.f90` (Northeast, Central States, Lake States) was updated to pass the tree diameter. The other 21 variant mortality routines kept the five-argument call, so under GOMPIT (`FVS_GOMPIT=1`) they read `DBHV` from whatever the sixth argument slot held.

On firebreather (gfortran 10.2.1, Debian) this is not benign. The Acadian library built from the five-argument source segfaults on the first scenario of the 36-cell synthetic matrix under GOMPIT, reproducibly (`ab_acd.log`, `ab_acd2.log`, exit 139), and runs 36 of 36 after the fix (`acd/seed_site_gate_gompit_postfix.csv`). Native mortality is unaffected: the pre-fix and post-fix Acadian native runs are identical in every value (`acd/engine_bakuzis_100yr_native_prefix.csv` against `_postfix.csv`), and the Northeast GOMPIT run, whose routine already passed the diameter, is identical pre and post fix (`ne/`).

Consequence: any GOMPIT result produced through fvs-modern on a variant other than NE, CS and LS between July 25 and October 28, 2026 came from an undefined read and should be rerun. The Paper 2 five-arm near-term run (August 4) is Northeast only and is not affected. The fvs-conus-pipeline projector implements GOMPIT in its own R code and does not call this routine.

The fix is one token in each of 21 files: `CALL GOMPSURV(..., SURVG)` becomes `CALL GOMPSURV(..., SURVG, D)`, where `D` is the tree's DBH already assigned as `D = DBH(I)` in every routine before the call (`vwc/morts.f90` passes `1.0` for the interval and now `D` as well). With the current five-column coefficient file `b5 = 0`, so the term is inert and the fix changes no number where the old code happened to run.

Harness: `validation/bakuzis_stress_v2/engine_bakuzis_v2.py`, 36 scenarios, 100 yr at 5 yr, libraries built with `deployment/scripts/build_fvs_libraries.sh` (the `-ldl` variant) from the same tree with and without the 21-file change.
