# Known Issues

## Recently Resolved

### PN/SN/IE keyrdr.f90 EOF blocker (resolved 2026-05-09)

The Pacific Northwest, Southern, and Inland Empire shared libraries previously
loaded via Python ctypes but hit a Fortran runtime EOF error in
`base/keyrdr.f90` line 47 when reading any keyword file. Documented in the
v2026.05.2 release notes as the remaining blocker for these variants.

Root cause was identified as five files in `src-converted/rd/` (RDADD, RDARRY,
RDCOM, RDCRY, RDPARM) that are pure Fortran INCLUDE files (variable
declarations and COMMON-block definitions, no SUBROUTINE wrapper) referenced
via `INCLUDE 'RDPARM.f90'` etc. in dozens of other rd/ sources. They have
`.f90` extensions so the build script also compiled them as standalone units,
polluting the link namespace with COMMON-block declarations and corrupting
keyword reader state at runtime.

Fix lives on branch `build-fixes-2026-05-06` (commits fb60191, 685d2c2,
ba10be9):
1. Source-list parser skips the 5 INCLUDE-only rd/ files
2. New `build_stubs()` function generates `libfvs_stubs.so` and
   `libfvs_stubs_final.so` from sources at `src-converted/stubs/`
3. F77-to-F90 conversion bug fixed in `cmdline.f90` and `apisubs.f90`
   (incorrect `, bind(c)` on local declarations of subroutine arguments)
4. VARVER stub added for the Acadian variant

Verification: `env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib`
produces 25 .so files plus 2 stubs; 23 of 25 load via ctypes; PN/SN/IE all
print the FVS variant banner under a keyword-file invocation. The 2 still-
failing variants (bc, on) are Canadian and have separate-scope compilation
gaps that were not in the April-May production set either.

## Regression Test Failures

### iet03 segfault (Inland Empire variant, test 03)

The iet03 regression test previously produced a segmentation fault during
execution, inherited from upstream USDA FVS source (Open-FVS revision 3360).
As of 2026-04-21 the test exits cleanly with `STOP 10` on the updated base
and fire extension sources (the SUMOUT, OPADD, OPCSET, OPGET3, FILOPN
restorations land before the FFE snag initializer fires). The summary output
diverges numerically from the 2025-04-25 baseline, which is expected given
the fire and SDI plumbing changes between those snapshots. The regression
harness currently counts this as a pass (simulation completed, summary lines
differ from stale baseline).

**Status**: Resolved for crash; baseline refresh pending before we claim
exact-match parity. Does not affect core growth/mortality projections or the
Bayesian calibration pipeline. Tracked by GitHub issues #3 and #5.

**Triage notes**: The iet03 keyword file exercises the Fire and Fuels
Extension (FFE) together with SNAGINIT, SIMFIRE, SALVAGE, DEFULMOD, POTFIRE,
FMORTMLT, and the full BurnRept/FuelOut/FuelRept/SOILHEAT reporting stack,
seeded only with a TREEFMT descriptor and no TreeInit tree list. The segfault
is consistent with FFE snag and fuels initialization running before the tree
list is populated. See `src-converted/tests/FVSie/iet03.key`.

**Planned investigation**:

1. Run `FVSie` under `gdb --args FVSie --keyword=iet03.key` to capture the
   faulting frame; expected location is the FFE snag or pot-fire initializer.
2. Compare against upstream Open-FVS r3360 to confirm the fault is inherited
   (it is listed as preexisting but never pinned to a commit).
3. If reproducible upstream, file an issue on the USDA Open-FVS tracker and
   keep the test skipped here with an explicit xfail marker in
   `run_regression_tests.sh`.

**Workaround**: Avoid the specific stand initialization configuration used in
iet03. Normal FVS runs with standard StandInit/TreeInit inputs are unaffected.

## Calibration Pipeline

### Crown ratio model performance

The crown ratio change model (script 05/05b) has the lowest predictive
performance among calibrated components, with R-squared values ranging from
0.13 to 0.38 across variants. This is consistent with the difficulty of
predicting crown dynamics from standard inventory variables. Crown ratio
predictions should be interpreted with appropriate uncertainty bounds.

**Planned improvement**: Explore lidar-derived crown allometrics and nonlinear
model forms in the CONUS variant development (branch: `conus-variant`).

### Ingrowth model is empirical only

The ingrowth component uses observed FIA recruitment rates by variant rather
than a mechanistic regeneration model. This is adequate for short-term (5 to
15 year) projections but may not capture long-term regeneration dynamics or
climate-driven species composition changes.

**Planned improvement**: Develop climate-sensitive ingrowth submodel for
the CONUS variant.

### ClimateSI and SDIMAX raster coverage

The `plot_raster_lookup.csv` file provides climate site index (ClimateNA-based)
and Emmerson SDIMAX values for FIA plot locations. Coverage is:

- ClimateSI: 99.6% of FIA conditions matched
- SDIMAX: 69.0% of FIA conditions matched

Conditions without raster matches fall back to FIA SICOND (for site index)
or default variant-level SDIMAX values.

## Build System

### ACD variant on cardinal.osc.edu (resolved 2026-04-21)

Earlier development cycles reported ACD build instability on Cardinal.
The resolution thread closed on 2026-04-21 with a fresh rebuild under
gcc/12.3.0 on Cardinal login01:

- `src-converted/acd/blkdat.f90` and `src-converted/acd/crown.f90` are
  byte-identical on Cardinal and in the workspace
  (sha256 `59db704b6857284e05281b721093e48647a664e7bcbead1f76aa67febe8f791d`
  and `33415929d3e442d6c448cf0b6867211e7901dbdacdd4e20cbfafe4b75b9efb93`).
- All MAXSP-dimensioned DATA statements in ACD carry 108 values, matching
  the declared MAXSP in `src-converted/acd/common/PRGPRM.f90`.
- ifort 2021.10.0 and ifx 2023.2.3 previously compiled ACD sources cleanly
  on Cardinal. gcc/12.3.0 now does the same: the 2026-04-21 Cardinal build
  produced `FVSacd.so` at 7.7 MB
  (sha256 `357ac26b51a5dc18804d7f764c43b609bbac7b1f46bc3991cf94dadddf9105af`),
  loaded via `ctypes.CDLL` with `RTLD_LAZY`, and exposes all four public
  API entry points (`fvssetcmdline_`, `fvssummary_`, `fvsdimsizes_`,
  `fvstreeattr_`) under `nm -D --defined-only`.
- The earlier gfortran 11.4.1 `-fimplicit-none` smoke-test artifact for
  `iosum` in `src-converted/common/OUTCOM.f90` is not triggered by the
  production build path (`build_fvs_libraries.sh`, which uses implicit
  typing). It remains tracked as an IMPLICIT NONE hardening item for
  the next calendar tag alongside `FMCOM.f90`.

Status: **resolved**. ACD now ships as a fully supported variant in the
REST surface and the shared-library matrix. `diagnose_acd_cardinal.sh` is
retained in `deployment/scripts/` for future triage should the discrepancy
recur on a different Cardinal login node or compiler revision.

### macOS compilation

On Apple Silicon Macs, the Homebrew gfortran installation may require explicit
library path configuration. See `deployment/scripts/setup_macos.sh` for the
recommended setup procedure.

### Windows native builds

FVS-modern is not tested on native Windows. Windows users should use WSL2
(see `deployment/scripts/setup_wsl.sh`) or Docker.
