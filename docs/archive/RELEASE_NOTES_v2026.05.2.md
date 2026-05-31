# Release notes — v2026.05.2 (2026-04-26)

## What's new

This release ships substantial build-infrastructure improvements
to the FVS-modern project. The FVS-PN, FVS-SN, and FVS-IE shared
libraries now build cleanly and load via Python ctypes for the
first time since the F77-to-F90 source conversion, removing a
foundational blocker for the FIA Bakuzis evaluation on western
and southern variants.

## Build infrastructure improvements

- `deployment/scripts/build_fvs_libraries.sh`: include-order fix
  resolves the `vwc/morts.f90` shared-source conundrum that
  silently broke PN and SN compilation. Variant-specific include
  paths now precede `-Ibase` so that `PRGPRM.f90` and other
  variant-aware headers resolve correctly.
- Repaired seven F77-to-F90 conversion bugs in `src-converted/econ/`
  (ecvol.f90, ecinit.f90, echarv.f90, ecstatus.f90, eccalc.f90,
  ecsetp.f90, ecin.f90). The conversion had mangled C-prefixed
  comments into broken trailing-`&` continuations, leaving the
  econ extension uncompilable.
- Twelve previously-undefined symbols are now resolved in the
  PN, SN, IE shared libraries: `morcon_`, `ecvol_`, `ecinit_`,
  `echarv_`, `ecstatus_`, `setpretendstatus_`, `eccalc_`,
  `ecsetp_`, `eckey_`, plus three downstream.

## New tools

- `scripts/fix_f77_comment_bug.py`: heuristic auto-repair for the
  comment-conversion pattern. Walks any Fortran file looking for
  trailing-`&` followed by stripped-comment fragments and rewrites
  them as proper `!` comments.
- `calibration/python/marshall_to_fia_csv.py`: adapter that
  converts Marshall-format FIA CSVs (Aaron's curated state-level
  exports) to the standard FIA DataMart schema. Already used to
  prepare per-state inventory data for OR, WA, AL, FL, GA, MS,
  SC, TN, ID, MT.

## Manuscript updates

- Section 5.4 limitations rewritten to accurately reflect the new
  FVS-PN/SN/IE state: shared libraries build cleanly and load via
  ctypes; the remaining blocker is a Fortran runtime keyword
  reader EOF that affects all PN/SN/IE keyfiles regardless of
  format or content.
- Section 5.5 future work adds the runtime regression as a sixth
  outstanding item, pointing readers at the diagnostic document.

## Known issue

FVS-PN/SN/IE libraries load via ctypes but hit a Fortran runtime
EOF error in `base/keyrdr.f90` line 47 when reading any keyword
file. This affects both our INVENTORY-mode keyfiles and the
canonical upstream USDA test files. Documented in
`calibration/python/PN_SN_LIBRARY_DIAGNOSIS.md`. Resolution
requires comparison against a USDA reference binary or deeper
Fortran I/O debugging — beyond the scope of this release.

## Operational notes

- Zenodo GitHub webhook for `holoros/fvs-modern` is still pending
  enablement. Once enabled, this tag (v2026.05.2) will mint a DOI.

## Verification

- `bash deployment/scripts/build_fvs_libraries.sh src-converted lib pn sn ie`
  on Cardinal with `FC=gfortran CC=gcc CXX=g++` produces working
  FVSpn.so, FVSsn.so, FVSie.so libraries.
- Python ctypes load test: `ctypes.CDLL("/path/to/FVSpn.so", RTLD_GLOBAL)`
  succeeds without OSError.
- nm symbol verification: `nm -D lib/FVSpn.so | grep -E " T (morcon|ecvol|ecinit)"`
  shows all three exported.

## Backwards compatibility

The build script change improves correctness without breaking
existing variant builds. Eastern variants (NE, ACD) that previously
worked continue to build identically.
