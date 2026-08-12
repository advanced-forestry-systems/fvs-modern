# Known Issues

## Limitations of the v1 code release (2026-08-11)

This release ships the Northeast model as frozen on 5 August 2026. The defects listed here were
diagnosed, sized, and deliberately left in place rather than repaired, and they are stated as
limitations of the release rather than as open bugs it is waiting on. Anyone quoting a benchmark
number from this codebase should read this section first. Full evidence sits in three project
records held outside this repository, in the CRSF project archive: the model freeze document of
5 August 2026 and its section 4 defect table, the mortality kernel validation of 11 August 2026, and
the unity multiplier diagnosis of 11 August 2026. Contact the maintainer for copies.

**Table 1.** Defects carried by the v1 code release, with the evidence and the disposition taken.

| Defect | Evidence | Disposition in v1 |
|:---|:---|:---|
| Background mortality is a diameter only logistic, halved in source | `vls/morts.f90:508` computes `RI=(1.0/(1.0+EXP(B0+B1*D)))`, then line 512 applies `RI = 0.5 * RI` | Documented, not fixed. The `FVS_MORT_HALVE` hook exposes the constant without changing the default, and the validated interim level of 2.675 is available but not adopted |
| Density gate that 78 percent of plots never reach | `PMSDIL = 55.` at `ne/grinit.f90:245` and `PMSDIU = 85.` at line 259; only 22 percent of the cohort crosses the trigger | Documented. `FVS_MORT_PMSDIL` exposes the threshold; closing it properly needs a continuous density dependent hazard |
| Crown ratio deficit sits in the allocator, not the rate | A crown ratio term in the background logit with a 6.7 fold range moved the most suppressed class by 0.006 | Landed as the environment gated `FVS_MORT_VMCRB` and `FVS_MORT_VMCRREF` hooks, default off, so frozen behaviour is bit identical unless a caller opts in |
| Reineke slopes fail for the native arms | -1.007 for `fvs_base` and -1.197 for `fvs_regional` against -1.605 expected on an acceptance band of -2.2 to -1.2, with R squared above 0.93 throughout | Documented. The 11 August mortality kernel validation shows the level fix overshoots past the steep edge of the band, so no tested setting closes this |
| Stocking guide passes only because of a clamp | `vls/morts.f90:750-786` inflate mortality until residual basal area is within 1 ft² ac⁻¹ of BAMAX; 29.6 percent of `fvs_base` and 39.5 percent of `fvs_regional` plot cycles pinned | Documented. Report the Bakuzis stocking guide pass as conditional and never as evidence of fidelity |
| Greg arms under predict diameter growth | Mean annual diameter growth 0.0952 in yr⁻¹ against 0.1166 observed on the HTCD screened basis and 0.1063 unscreened, so 10 to 18 percent low depending on the screen | Documented. Closing needs a refit or recalibration of the Greg diameter growth coefficients on an Acadian cohort with a working BAL term |
| Native arms hold physically impossible basal areas indefinitely from an overstocked start | Basal area as a fraction of the BAMAX ceiling runs flat against the clamp over much of the 440 plot, 150 year projection in both native arms | Documented. This is the joint symptom of the three mortality rows above, not a separate defect to patch |
| Greg arms annihilate the stand over a long horizon | Fractional stem loss of 0.837, 0.796, 0.792 and 0.789 from poorest to best site, flat and slightly backward across site; relative density falls from 0.396 to about 0.30 so the FVS density machinery never engages and the GOMPIT override alone drives the loss | Documented. GOMPIT is density aware through crown closure but carries no productivity term and needs a site covariate |
| `conus_climate` height coefficient coverage is partial, with silent fallback | Coverage on the benchmark cohort is 18.91 percent of basal area and 21.81 percent of stems; the remainder falls back to the native height model without a warning | Documented, with a correction. The 11.14 percent basal area and 16.49 percent stem figures that appear in the working record and the stress battery memo were measured on a small unrepresentative sample and are superseded by the two values given here |
| `config/config_loader.py` six column offset persists in three writers | MORTMULT was fixed at line 1177; BAMAX at 1144 and 1150, BAIMULT at 1191 and HTGMULT at 1201 still write the keyword padded to 16 columns | Documented. Applying the MORTMULT layout to the remaining writers needs the owning track's agreement because BAIMULT and HTGMULT govern growth |
| StandID aliasing collapses 3,321 control numbers to 50 identifiers | `f"{variant.upper()}_{t1_cn % 10_000_000:07d}"` puts 3,320 plots in a collision group, largest group 147 plots | Documented and provably latent, since each `run_stand` call uses its own `tempfile.mkdtemp` and nothing downstream keys on StandID |


### Calibrated configuration defects diagnosed 11 August 2026

Three further defects were established while auditing the calibrated configuration files, after the
freeze document was written. They are not in Table 1 above because they sit in the calibration and
serialization path rather than in the Fortran model, but they change how a benchmark number from
this release should be labelled.

First, the per species diameter growth and height growth multipliers are absent from every Northeast
calibrated configuration. In `config/calibrated/ne.json`, and identically in `calibrated_sdifix/ne.json`,
`dds_multiplier` and `htg_multiplier` are exactly 1.0 for all 108 species slots, while the three
species code keyed components in the same block carry real variation. The underlying posteriors are
not degenerate: the fit holds 81 species level diameter growth terms with a standard deviation of
0.0976, implying multipliers from 0.7459 to 1.3211, and 27 height increment terms implying 0.8765 to
1.1194. Real fitted structure was discarded rather than correctly reported as absent. The mechanism is
a dense index mapping failure in the serializer, which needs a `<component>_species_index.csv`
crosswalk that was never written; the fit wrote `diameter_growth_map.csv` and `height_increment_map.csv`
instead, which are variable and estimate dumps rather than crosswalks. A unity array is not a harmless
no operation here, because `_array_or_none` accepts it as valid input and thereby suppresses the legacy
fallback multiplier routines.

Second, the emitter that produced the shipped configurations is not under version control. The
provenance schema present in `ne.json` does not match what the current `R/multipliers.R` emits, and no
file in either repository writes the string `mort_n_re`. The shipped configurations therefore cannot be
regenerated byte identically at present, and the mechanism above is inferred rather than replayed,
though the inference is strong because the counts of 81 and 27 are exactly the fitted term counts.

Third, the adoption record contradicts itself. The header comment in `R/multipliers.R` says diameter
growth was adopted for 7 variants and height increment for 6, while the authoritative
`calibration/data/equation_availability_full.csv` marks diameter growth TRUE for all 25 variants and
height increment TRUE for 6. The emitted configurations follow neither. Someone has to decide whether
18 variants are missing calibration they should have, or 7 received calibration they should not.

What this changes for the release is a labelling matter and it is material. The `fvs_regional` gate
numbers already published come from a model calibrated on height diameter, crown ratio and mortality
only, running default diameter growth and default height growth. Any claim about regional calibration
performance on volume or basal area is currently a claim about three of five components. The volume
gate still passes on its own terms. Whether a repair would move it is open, since the Northeast
diameter growth multiplier median is 0.9965 and pooled aggregate bias may barely shift, but the species
composition of that bias would change and the Reineke and Eichhorn diagnostics are sensitive to exactly
that. Repair was deliberately kept out of v1 because it requires a recovered or refitted crosswalk plus
re derivation of every downstream configuration and gate number.

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
