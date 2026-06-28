# GOMPMC v2: native engine spec for our richer gompit survival (options 2 and 3)

Design for extending fvs-modern's native gompit mortality (base/gompmort.f90, GOMPMC common) from Greg's
5-parameter form to our deployed survival_unified_v2_crz form, so options 2 (species-dependent) and 3
(species-free) adopt natively the same way option 1 (Greg) now does. Grounded in the validated GOMPSURV path
(this turn: native formula reproduces the R reference to 1.6e-7) and the survival_unified_v2_crz model that
already carries z_sp.

## 1. The form to support

Deployed survival eta (per tree). NOTE the sign: the hazard is exp(-eta), so higher eta means higher
survival (survival_unified_v2_crz.stan lines 12-13: annual_hazard = exp(-eta); P(survive T) = exp(-exp(-eta)*T)):

  eta = b0 + trait_effect[sp] + z_sp[sp]
      + z_L1[L1] + z_L2[L2] + z_L3[L3] + z_FT[fortype]
      + b1*dbh + b2*dbh^2
      + b3*cr_z + b3b*cr_z^2
      + b4*ln_csi
      + b5*bal_metric
      + b6*sqrt_ba_rd
      + b7*cch_z + b7b*cch_z^2

  trait_effect[sp] = W[sp,] . gamma   (8 traits)
  cr_z  = (cr  - cr_mean ) / cr_sd
  cch_z = (cch - cch_mean) / cch_sd

Option 3 (species-free, b1) uses trait_effect only; option 2 (species-dependent, b2) uses trait_effect +
z_sp. Everything else is identical.

## 2. The key simplification: collapse to per-species and per-stand scalars offline

Most terms do not need to live in the engine as separate effects. Precompute, offline, two scalars:

  a_sp[sp]    = trait_effect[sp] + z_sp[sp]      (option 2)   or   trait_effect[sp]   (option 3)
  a_stand     = z_L1[L1] + z_L2[L2] + z_L3[L3] + z_FT[fortype]

a_sp is fixed per FVS species and ships in the coefficient file (exactly like Greg's b0..b4 do today).
a_stand is resolved ONCE per stand at setup from the stand's EPA L1/L2/L3 ecoregion codes and FORTYPCD, then
held constant for the run. b0 folds into a_sp (add b0 to every a_sp) so the engine carries no separate
intercept. The slopes b1..b7b are GLOBAL fixed effects (identical across species), so they are stored once,
not per species.

In-engine evaluation then reduces to:

  eta = a_sp[isp] + a_stand
      + BV(1)*dbh + BV(2)*dbh^2
      + BV(3)*cr_z + BV(4)*cr_z^2
      + BV(5)*ln_csi
      + BV(6)*bal_metric
      + BV(7)*sqrt_ba_rd
      + BV(8)*cch_z + BV(9)*cch_z^2

with cr_z, cch_z standardized using four constants (cr_mean, cr_sd, cch_mean, cch_sd). This is the same shape
as the current GOMPSURV (a per-species intercept plus a handful of slopes), so the engine change is modest.

## 3. GOMPMC v2 common-block layout

Extend GOMPMC.f90 (keep v1 fields for Greg's option; add v2 in parallel, selected by a mode flag):

  INTEGER  NGOMP, GMODE              ! GMODE: 1 = Greg b0..b4 (v1), 2 = surv_crz (v2)
  REAL     GB(MAXSP,5)               ! v1 (unchanged)
  REAL     ASP(MAXSP)                ! v2 per-species scalar a_sp (includes b0)
  LOGICAL  GHAVE(MAXSP)
  REAL     BV(9)                     ! v2 global slopes b1,b2,b3,b3b,b4,b5,b6,b7,b7b
  REAL     STDZ(4)                   ! cr_mean, cr_sd, cch_mean, cch_sd
  REAL     ASTAND                    ! v2 per-stand addend, set in GOMPLOAD/stand setup
  REAL     CCHT(MAXTRE)              ! per-tree cch on the gompit scale (unchanged)

## 4. Coefficient files (emitted by config_loader / an R emitter, per option)

Three small files, mirroring how 62h emits Greg's CSV:
1. a_sp file: header + rows SPCD, a_sp        (option 2 uses trait+z_sp; option 3 uses trait only)
2. global file: BV(1..9) and STDZ(1..4)       (one row; identical across stands)
3. stand-effect lookup: ecoregion L1/L2/L3 and FORTYPCD -> z components, summed to a_stand at stand setup
   (the projector already carries these RE tables: surv_crz *_re_L1/L2/L3/FT.csv, and z_sp now via the
   re_sp.csv export added this turn).

Because a_sp absorbs trait_effect + z_sp + b0, the engine never sees traits, gamma, or the W matrix; all of
that is reduced offline. Switching options 2 vs 3 is just swapping the a_sp file, exactly like the config JSON
blocks differ today.

## 5. Per-tree inputs the engine must supply (availability)

  dbh           tree record (have)
  cr            tree record (have; same as v1)
  cch           GOMPCCH each cycle (have; v1 already does this)
  ln_csi        stand climate site index -> log; from the stand record / climate attach (stand-level, set once)
  bal_metric    BAL_SW + BAL_HW basal-area-in-larger (FVS computes BAL; map softwood/hardwood split)
  sqrt_ba_rd    sqrt( BA*0.2296 * RD ), RD = SDI / SDImax; FVS has BA and SDI; SDImax from the stand
  dbh^2, cr_z, cch_z, cr_z^2, cch_z^2   derived in-engine from the above + STDZ

Two need wiring beyond v1: ln_csi (a stand climate scalar) and the BAL softwood/hardwood split for bal_metric.
Both are stand- or tree-level quantities FVS already tracks; the work is mapping them, not computing new state.

## 6. GOMPSURV v2 (sketch)

  SUBROUTINE GOMPSURV2(ISPC, DBH, CR, CCHV, BALSW, BALHW, BA, SDI, SDIMAX, LNCSI, FINTL, SURV)
    eta = ASP(ISPC) + ASTAND
    crz  = (CR  - STDZ(1)) / STDZ(2)
    cchz = (max(CCHV,0) - STDZ(3)) / STDZ(4)
    rd   = SDI / max(SDIMAX,1e-3)
    sbar = sqrt( max(BA*0.2296,0) * max(rd,0) )
    eta = eta + BV(1)*DBH + BV(2)*DBH*DBH + BV(3)*crz + BV(4)*crz*crz
              + BV(5)*LNCSI + BV(6)*(BALSW+BALHW) + BV(7)*sbar + BV(8)*cchz + BV(9)*cchz*cchz
    clamp eta to [-30,30]; HZ = exp(-eta); SURV = exp(-HZ*FINTL)   ! NOTE exp(-eta): high eta -> high survival

Sign matters: surv_crz uses hazard = exp(-eta), the OPPOSITE of v1 gompmort.f90 as currently coded (which uses
exp(+eta) and is itself wrong for Greg's gompit; see the stress-test finding and task #24). The morts.f90 hook
(IF(LGOMP.AND.GHAVE(ISPC)) ... GOMPSURV ...) calls GOMPSURV2 when GMODE==2, passing the extra tree/stand
arguments it already has in scope. Add the same benign/stress biological-plausibility gate as v1 before enabling.

## 7. Validation plan (mirrors this turn's v1 check)

1. Unit: compile GOMPSURV2 standalone, feed the a_sp + global files, compare eta and S to the R surv_crz
   prediction for a tree grid. Target agreement < 1e-4 (single vs double), as achieved for v1.
2. Stand: once a working NE keyword fixture exists (the in-repo net01.key is a broken echo; FVS01 INVALID
   KEYWORD, a maintainer fix), run the A/B (native default vs GMODE=2) and confirm density regulation.

## 8. Effort and sequencing

- Offline emitters (a_sp, global, stand-effect): small R, reuse the projector's coefficient tables + the new
  re_sp.csv. Low risk, do first.
- GOMPMC.f90 + GOMPSURV2 + the morts.f90 branch: a focused Fortran change patterned on the working v1 trio;
  the maintainer or a careful patch. Medium.
- ln_csi and BAL softwood/hardwood wiring: the only genuinely new engine plumbing. Medium.

Net: options 2 and 3 reach native mortality by the same route option 1 just did, with the richer eta reduced
to a per-species scalar plus nine global slopes plus a per-stand addend. The science is unchanged; this is the
in-engine wiring of the already-fitted model.
