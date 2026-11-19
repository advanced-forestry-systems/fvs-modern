# Mortality models

Every variant library carries two mortality paths: the variant's native mortality routine and the GOMPIT crown-closure survival model. GOMPIT is off by default and is switched on at run time without recompiling. The configuration catalog names these two paths `native` and `gompit` (see [Configuration files](configuration-files.md#catalog-tables)).

## Native variant mortality

Native mortality lives in each variant's `morts.f90`, called from `TREGRO` when cycling for growth. Some variants share a file through the source lists: NE, for example, compiles `src-converted/vls/morts.f90`, and `vwc/morts.f90` serves both PN and WC. The routine computes a per-tree periodic mortality, reduces each record's trees per acre, and has an entry point `MORCON` that loads one-time constants.

In `vls/morts.f90` the per-tree rate combines a background rate and a density-driven rate. The header describes the density part: SDI-based mortality applies while the stand diameter metric is below 10 inches, after which BAMAX-based mortality takes over, and BAMAX defaults to the value implied by maximum SDI at 10 inches. After the tree loop, density mortality is redistributed by the variant routine `VARMRT`, and a loop adjusts mortality until residual basal area is within 1 ft² ac⁻¹ of `BAMAX`.

The NE background rate is a Hamilton logistic in diameter, `RI = 1/(1 + EXP(B0 + B1*D))`, followed by a legacy halving. That halving and an optional crown-ratio term are now runtime parameters with defaults that reproduce legacy behavior:

| Environment variable | Default | Effect | Location |
|---|---|---|---|
| `FVS_MORT_HALVE` | 0.5 | Multiplier on `RI` (the former hardcoded halving) | `vls/morts.f90` |
| `FVS_MORT_CRB` | 0.0 | Coefficient on `(CR - CRREF)` added to the background logit | `vls/morts.f90` |
| `FVS_MORT_CRREF` | 45.0 | Crown ratio centring constant, percent | `vls/morts.f90` |
| `FVS_MORT_PMSDIL` | 55 | Lower SDI percentage threshold `PMSDIL` | `ne/grinit.f90` |
| `FVS_MORT_VMCRB`, `FVS_MORT_VMCRREF`, `FVS_MORT_RELCAP`, `FVS_MORT_VMDIAG` | | Hooks in the NE density allocator | `ne/varmrt.f90` |

The MORTMULT keyword is processed at the top of the NE routine (activity code 94) into per-species `XMMULT`, `XMDIA1` and `XMDIA2`, and applied as `XMORT` to trees with `D .GE. D1 .AND. D .LT. D2`.

`KNOWN_ISSUES.md` lists the known limitations of the native NE path as shipped in the v1 code release, including the halved background rate, a density gate (`PMSDIL = 55`, `PMSDIU = 85` in `ne/grinit.f90`) that most plots never reach, and Reineke slopes outside the acceptance band.

## GOMPIT survival

GOMPIT replaces native mortality per tree, per cycle, inside the growth loop, for species that have a fitted coefficient row. Species without a row keep native background mortality. The code is in `src-converted/base/gompmort.f90`; shared state is the include `src-converted/common/GOMPMC.f90`.

| Routine | Role |
|---|---|
| `GOMPLOAD` | Called from `MORCON`. Resolves activation once per run, reads the coefficient file, maps rows onto FVS species, reads the two tuning variables. |
| `GOMPCCH` | Recomputes crown closure at each tree's tip for the current cycle from the live tree list, using an ORGANON SWO crown-width port, and maps it with `CCH = 0.062 + 0.0036 * cch_hat`. |
| `GOMPSURV(ISPC, CR, CCHV, FINTL, SURV, DBHV)` | Survival of one tree over a cycle of `FINTL` years: annual survival raised to `FINTL`. |
| `GOMPON` | Hook for the `GOMPMORT` keyword; sets `LGOMPKW`. |
| `BLOCK DATA GOMPBD` | Initialises `LGOMP`, `LGOMPKW`, `GHAVE`, `NGOMP`, `GGRP`, `GSFLOOR` (0.0) and `GMORTMULT` (1.0). |

### Six-argument GOMPSURV

`GOMPSURV` has taken six arguments since commit 495541e (2026-07-25). The sixth, `DBHV`, is the tree diameter and enters as a size term `B5*LOG(MAX(DBHV, 0.1))`. In PR 119 the 21 variant `morts.f90` call sites were updated to pass `D`. At this commit 22 `morts.f90` files call it with six arguments: `acd`, `adk`, `ak`, `bm`, `ca`, `canada/on`, `ci`, `cr`, `ec`, `em`, `ie`, `kt`, `nc`, `oc`, `op`, `sn`, `so`, `tt`, `ut`, `vls`, `vwc` and `ws`. When adding or merging a variant `morts.f90`, check that its `GOMPSURV` call passes all six arguments.

### How GOMPIT is switched on

Activation is resolved in `GOMPLOAD` and needs two things:

1. **A switch.** Either the environment variable `FVS_GOMPIT` is set to a value whose first character is not `0`, `n`, `N`, `f` or `F`, or the keyword file contains `GOMPMORT` (handled in `vbase/initre.f90`, which calls `GOMPON`; `GOMPMORT` is registered in `base/keywds.f90`).
2. **A coefficient file.** `FVS_GOMPIT_COEF` must name a readable file. If it is empty, or the file cannot be opened, GOMPIT stays off and a message goes to the main output listing (`GOMPMORT: ... gompit mortality NOT enabled.`).

On success the listing reports `GOMPMORT enabled: <n> fitted species read, <m> matched to this variant.` and the cap settings. The keyword records the on or off choice in the keyfile; the file path always comes from `FVS_GOMPIT_COEF`.

```bash
export FVS_GOMPIT=1
export FVS_GOMPIT_COEF=/path/to/fvs-modern/config/greg_mortality_coefficients.csv
```

`fvs2py.catalog.select_config` returns exactly these two variables in `choice.env` for gompit arms.

Two optional variables tune GOMPIT:

| Variable | Default | Effect |
|---|---|---|
| `FVS_GOMP_MORTMULT` | 1.0 | Multiplies `exp(eta)` inside the annual survival `1 - exp(-exp(eta) * MORTMULT)`, so annual mortality becomes `m**MORTMULT`: values above 1 lower mortality and values below 1 raise it (0.5 gives the square root of `m`, not half). Values of 0 or less reset to 1.0. |
| `FVS_GOMP_SFLOOR` | 0.0 (off) | Floor on annual survival for all trees; period survival is floored at `SFLOOR**FINTL`. Values outside 0 to 1 reset to 0. |

`GOMPMC.f90` records that these two substitute for the MORTMULT keyword, which does not reach GOMPIT. In `vls/morts.f90` the GOMPIT branch sets the trees killed directly, and `VARMRT` is skipped when GOMPIT is active.

### Coefficient file

`config/greg_mortality_coefficients.csv` is the default file named by `fvs2py.catalog` (`GOMPIT_COEF_DEFAULT`). Its header is `SPCD,n,b0,b1,b2,b3,b4` and it has 92 species rows. `GOMPLOAD` reads the header line, then each row as `SPCD, n, b0, b1, b2, b3, b4, b5` and falls back to the seven-field form with `b5 = 0` when the size term is absent, so the shipped file is size-blind. It reads at most 600 rows. Rows are matched to FVS species through `FIAJSP`; each FVS species also gets an ORGANON group (1 for FIA codes below 300, otherwise 16) that `GOMPCCH` uses for crown geometry.

`config/greg_mortality_coefficients_cch2.csv` has the same header. The `config/greg_mort_coefficients_<driver>.csv` files (`bgi`, `cspi`, `elev`, `esi`, `none`) use a different header (`SPCD,n,deaths,b0,b_lnDBH,b_CR,b_CCH,b_driver`); their columns do not line up with the `b0` to `b5` positions `GOMPLOAD` reads, so do not point `FVS_GOMPIT_COEF` at them.

### The survival computation in code

The executable lines of `GOMPSURV` are:

```fortran
ETA = B0 + B1*(CRC+0.01)**B2 + B3*CTERM + B5*LOG(MAX(DBHV,0.1))
IF (ETA.GT.30.0)  ETA = 30.0
IF (ETA.LT.-30.0) ETA = -30.0
HZ = 1.0 - EXP(-EXP(ETA) * GMORTMULT)   ! annual survival (named HZ historically); GMORTMULT multiplies exp(eta)
SURV = MAX(0.0,MIN(1.0,HZ)) ** FINTL   ! compound annual survival over the cycle
```

`CRC` is crown ratio clamped to 0.0001 to 1, and `CTERM` is `CCHV**B4` for positive crown closure and 0 otherwise. Callers pass the cycle length as `FINTL` and kill `PROB * (1 - SURV)` trees. `vwc/morts.f90` instead calls with `FINTL = 1.0`, sets the native annual rate to `1 - SURV`, and lets its own period conversion apply.

!!! success "Fixed in PR 143 (October 30, 2026)"
    The comments in `gompmort.f90` and `GOMPMC.f90` now describe what the code computes: annual survival `1 - exp(-exp(eta) * GMORTMULT)` raised to the cycle length, with `FVS_GOMP_MORTMULT` acting as an exponent on annual mortality. Before this fix the header described period survival `exp(-exp(eta) * T)`, and the comments said a multiplier of 0.5 halves mortality, when in fact it raises it. No executable line changed.

!!! warning "Framing differs between the engine and the fitting code"
    The fitting script `calibration/R/35_fit_greg_mortality_cch.R` and `calibration/python/greg_mortality.py` still write period survival as `exp(-exp(eta) * T)`. The engine computes annual survival `1 - exp(-exp(eta) * GMORTMULT)` raised to the cycle length. These are different functions of `eta`, and the repository does not state which framing the shipped coefficients are on. Treat the engine code as authoritative for what it does, and confirm the framing before fitting new coefficients or comparing engine output with the Python implementation.

## Choosing between them

Use the catalog rather than choosing by hand. `config/catalog/catalog_recommendation_by_variant.csv` names, per variant, which of `default_native`, `default_gompit`, `calibrated_native` and `calibrated_gompit` passed the long-horizon gates. See [Validation](validation.md).
