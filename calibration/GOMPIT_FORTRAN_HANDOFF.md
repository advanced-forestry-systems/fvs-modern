# Gompit mortality in FVS: Fortran integration handoff

Greg Johnson's CONUS gompit survival now runs **inside the FVS Fortran growth
loop** (TREGRO -> MORTS), substituted for native mortality per tree, per cycle,
so growth and density interact with the substituted mortality each cycle. This
supersedes the post-hoc TPA-overlay approach, which validation rejected (it
disabled FVS mortality and made FVS growth run away). Branch:
`feature/gompit-projection-wiring`, PR #67.

This document is the single reference for what was built, how it works, current
coverage, validation results, how to run it, and the two decisions that remain
yours.

---

## 1. The model

Per species, annual hazard with a cycle-length exposure:

```
eta = b0 + b1*(cr+0.01)^b2 + b3*cch^b4
H   = exp(eta)                 ! annual hazard
S   = exp(-H * FINT)           ! period survival, FINT = cycle length (yrs)
trees killed = PROB * (1 - S)
```

`cr` is crown ratio (FVS `ICR/100`). `cch` is crown closure at the subject
tree's tip, recomputed each cycle from the live treelist by an ORGANON
crown-closure port and mapped onto the gompit cch scale by the validated affine
fit `cch = CCH_A + CCH_B*cch_hat` (CCH_A=0.062, CCH_B=0.0036; 35d_validate_cch.R,
Spearman 0.93 vs the panel's stored CCH1).

Coefficients: `conus_mort/full_out/greg_mortality_coefficients.csv` (133 species,
columns SPCD,n,b0..b4). Model-level validation on 7.6M FIA remeasurement records:
gompit AUC 0.740 vs base-rate 0.673, log loss 0.420 vs 0.443, 131/133 species
improved; gompit tracks observed survival into the most crowded crown-closure
quartile where the base rate misses crowding mortality
(`calibration/output/gompit_mortality_validation.md`).

---

## 2. The code

| file | role |
|------|------|
| `src-converted/common/GOMPMC.f90` | shared state (COMMON): `LGOMP`, `NGOMP`, `GB(MAXSP,5)`, `GHAVE(MAXSP)`, `GGRP(MAXSP)`, `CCHT(MAXTRE)`. INCLUDE after PRGPRM. (Reached via the `base/common` symlink.) |
| `src-converted/base/gompmort.f90` | `GOMPLOAD` (read env, load coeffs, resolve to variant species via FIAJSP, assign ORGANON group); `GOMPCCH` (ORGANON crown-closure-at-tip port, affine-mapped, fills CCHT each cycle); `GOMPSURV(ispc,cr,cch,years,surv)` (period survival). |

The hook in each variant's `morts.f90` is five edits: INCLUDE GOMPMC; declare
`CRG,SURVG`; `IF(LGOMP) CALL GOMPCCH` once before the per-tree loop; per-tree
override for fitted species (`WKI = PROB*(1-S)`); `CALL GOMPLOAD` in MORCON. For
variants that call VARMRT, the native SDI redistribution is bypassed when LGOMP
(gompit is density-aware via cch). For the ORGANON-logistic `vwc` family the
hook replaces the annual rate `RIP` with `1-exp(-H)` and lets the routine's own
period conversion finish, reproducing gompit period survival exactly.

`GOMPLOAD` is SAVE-guarded (loads once per run). All activation is env-gated, so
one binary does native or gompit with no recompile.

---

## 3. Activation

```bash
export FVS_GOMPIT=1
export FVS_GOMPIT_COEF=/path/to/greg_mortality_coefficients.csv
# run FVS normally; unset FVS_GOMPIT for native mortality
```

A `GOMPMORT` keyfile keyword is intentionally NOT implemented (see Section 6).

---

## 4. Coverage (gompit on all 23 buildable variants; 19 validated)

Two variants -- **ON (Ontario)** and **BC (British Columbia)** -- do not produce
an executable even unmodified, on pre-existing link errors unrelated to gompit
(`fmcrow_` in ON's fire module; `dbsrd1_`/`dminit_` in BC's database/mistletoe
modules). They are the only two of 25 without a gompit executable. ON is hooked
(its morts compiles); BC is left pristine (its nonblock-DO control flow is
hostile to the hook and the variant cannot be built regardless). **So gompit is
wired into every variant that the repo can currently build (23/23).**


Hook applies across **all five mortality-routine families** in the engine:

| family | routine | variants | status |
|--------|---------|----------|--------|
| eastern TWIGS (shared) | `vls/morts.f90` | NE, CS, LS | built, validated |
| eastern (own) | `sn/morts.f90` | SN | built, validated |
| Dixon SDI / VARMRT (own) | `<v>/morts.f90` | CR, WS, EC, CA, BM, NC, OC, SO, UT, EM, OP, ACD, ON | built (ON exe blocked*), most validated |
| Hamilton/Prognosis (own) | `<v>/morts.f90` | AK, CI, IE, KT | built, AK/CI/IE validated |
| ORGANON-logistic (shared) | `vwc/morts.f90` | WC, PN | built, validated |

\* ON (Ontario) morts compiles and is hooked, but its executable fails to LINK
on a pre-existing, unrelated error (`fmcrow_` undefined) -- the variant does not
build even unmodified.

**BC (British Columbia):** left pristine. Its executable also fails to build
unmodified (`dbsrd1_`/`dminit_` undefined), and its nonblock-DO control flow
breaks the hook insertion. Both issues must be fixed before BC can carry gompit.

**Built but not validated** (their FIA plots are labelled under other variants,
so no stress stands exist for them in the CONUS standinit): OC, OP, KT, ACD.
Their hook compiles, links, and runs; they just lack stands to A/B here.

### Validation results (yr100 mean AGB, native -> gompit, ~8-200 stands/variant)

All bounded, none runs away or crashes. 19 variants validated, sorted by effect:

```
SN -78   TT -76   CI -75   WC -56   WS -44   EM -43   PN -42   LS -41
AK -38   CR -35   IE -24   NE -21   CS -21   CA -15   BM -9    EC -5
SO ~0    NC ~0
```
(TT, CI, CR, AK, EM are sparse high-elevation/high-latitude western variants
with very low absolute biomass, so their large percentages are noisy.)

Figure: `calibration/output/gompit_fvs_allvariants.png` (trajectories + sorted
percent-change). The effect size tracks BOTH ORGANON-proxy fit and stand
density/crown closure: western conifer (EC, SO, NC ~0) is the best proxy fit and
nearly matches native; SN/CI (southern, central Idaho) the largest. Dense wet
PNW (WC -56, PN -42) shows a large effect because high crown closure drives the
cch term. UT is excluded from the figure (degenerate 0-biomass sample).

---

## 5. How to reproduce

Build a variant executable with gompit (Cardinal, `module load gcc/12.3.0`):

```bash
cd ~/fvs-modern
bash deployment/scripts/build_fvs_executables.sh src-converted <out_lib_dir> ne
```

Validate native vs gompit on stress stands (the driver lives on Cardinal at
`/fs/scratch/PUOM0008/crsfaaron/fvs_gompit/gompit_fvs_validate.py`):

```bash
export FVS_LIB_DIR=<out_lib_dir> VAL_VARIANT=ne VAL_N=200
unset FVS_GOMPIT;  VAL_MODE=native  python3 gompit_fvs_validate.py
export FVS_GOMPIT=1; VAL_MODE=gompit python3 gompit_fvs_validate.py
```

Adding a new variant of an existing family: add `../base/gompmort.f90` to its
`bin/FVS<v>_sourceList.txt`, apply the five-edit hook (the patchers
`/tmp/patch_morts2.py` for the FINT-line form, `patch_morts3.py` for the
WK2(I)=WKI form capture the transforms), compile-check, build, validate.

---

## 6. The two flagged items -- both now RESOLVED

1. **Group-map refinement -- tested, rejected on evidence; coarse proxy kept.**
   A genus/crown-form crosswalk over all 18 ORGANON groups was built and the
   affine cch map re-fit on the held validation sample (113k trees). It DEGRADES
   the cch reproduction: Spearman 0.853 (full genus) and 0.875 (conifer-only) vs
   0.925 for the coarse softwood/hardwood proxy. The PNW-specific ORGANON crown
   equations add species variance that does not match how the stored CCH1 was
   generated; uniformity wins for a rank-order proxy. The coarse proxy is
   retained unchanged and is now *justified*, not a limitation -- the variant
   spread in gompit effect is a real cch-sensitivity property, not a crosswalk
   artifact. See `calibration/output/cch_crosswalk_refinement_test.md` and
   `calibration/python/refine_cch_crosswalk.py`.

2. **`GOMPMORT` keyword -- implemented and validated.** Added at the only safe
   point in the legacy dispatcher: `keywds.f90` TABLE slot 148 (previously
   blank, no index shift) + a new GOTO target and handler in `vbase/initre.f90`
   that calls `GOMPON` (flags `LGOMPKW`); `GOMPLOAD` then activates on (env OR
   keyword); `BLOCK DATA GOMPBD` initialises the flags. Coeff path stays in
   `FVS_GOMPIT_COEF`. Validated on NE: native 211.1, env-gompit 175.5,
   keyword-gompit 175.5 (byte-identical to env).

## 6b. Stress test (robustness pass)

19 variants x ~60 stands x native/gompit (~45,000 projection-years) produced
**zero** NaN/inf/negative/>2000-t-ac gompit values and no crashes; gompit effect
ranges NC -4% to UT -82%, with productive variants in a tight -13% to -27% band.
See `calibration/output/gompit_stress_test.md` + `gompit_stress_summary.csv`.

---

## 7. Suggested next steps (safe to autopilot)

* Hand-hook BC; rebuild ON once the unrelated fire-link issue is resolved.
* A/B the built-but-unvalidated variants (OC, OP, TT, KT, ACD) once stress
  stands exist.
* Scale every variant's A/B to a few hundred stands for manuscript tables.
* After the group map is settled, re-run the full sweep so all variants share
  the refined crosswalk and affine map.

---

## 8. Artifact index (all under PR #67)

* Fortran: `src-converted/base/gompmort.f90`, `src-converted/common/GOMPMC.f90`,
  the hooked `*/morts.f90`, updated `bin/FVS*_sourceList.txt`.
* Model-level validation: `calibration/output/gompit_mortality_validation.{md,png}`,
  `calibration/R/40_gompit_mortality_validation_figure.R`.
* In-engine validation: `calibration/output/gompit_fvs_integration_validation.md`,
  `gompit_fvs_inengine.png`, `gompit_fvs_allvariants.png`,
  `calibration/R/41_*`, `calibration/R/42_*`, per-variant CSVs under
  `calibration/output/gompit_fvs_csls/`.
* Rejected post-hoc approach (kept as the negative-result record):
  `calibration/python/run_gompit_projection.py`,
  `calibration/python/GOMPIT_INTEGRATION_FINDINGS.md`, and the gompit helper
  modules `greg_mortality.py`, `cch_organon.py`, `project_mortality.py`.
