# Greg Johnson's deployed CONUS DG / HG / mortality forms (recovered 2026-06-28)

This resolves the #28 blocker. The deployed (climate-based) equation algebra was NOT in a C++ file as
previously assumed. It is implemented, authoritatively, in the Greg projector on Cardinal, and documented in
the fvs_remodeling PDFs and quarto scripts (and on Google Drive "fvs remodeling").

## Sources
- AUTHORITATIVE (code): /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_greg.R
  (copied to this folder). The dg_annual / hg_annual / gomp_surv_annual functions are the deployed forms.
- Documentation: ~/fvs_remodeling/pdfs/{diameter_growth_equations_for_conus.pdf,
  Height_Growth_Equations_for_CONUS.pdf, Mortality_Equations_for_CONUS.pdf} and the matching .qmd fit scripts
  in ~/fvs_remodeling/scripts/{diameter_growth,height_growth,mortality}. Mirrored on the user's Google Drive
  "fvs remodeling" folder.
- Coefficients: ~/fvs_remodeling/rds/{dg_parms.RDS (84 spp, B0..B6), hg_parms.RDS (96 spp, B0..B8),
  mort_parm_base_rate_cr_cch.RDS (b0..b4)}.

## 1. Diameter growth (annual increment, inches/yr) -- dg_parms B0..B6
    z  = B0 + B1*log( (dbh+1)^2 / (cr*ht+1)^B3 ) + B2 * bal^B4 / log(dbh+2.7) + B5*elev + B6*EMT
    dg = exp( clamp(z, -30, 5) ),  dg >= 0
Inputs: dbh (in), cr (fraction), ht (ft), bal (ft2/ac basal area in larger trees), elev (ft), EMT (deg C).
Note: B3 and B4 are exponents inside the form (not linear slopes); B0..B6 are the 7 stored per-species params.

## 2. Height growth (annual increment, ft/yr) -- hg_parms B0=max_height, b1..b8
    dht = mx*b1*b2 * cr^b3 * exp( -b1*ht - b4*ccfl - b8*cch^0.5 - b5*elev + b6*sqrt(TD) + b7*EMT )
              * (1 - exp(-b1*ht))^(b2-1)
where mx = B0 = species max height. Inputs: ht (ft), cr (fraction), ccfl (ft2/ac crown competition factor in
larger trees), cch (crown closure at tip, fraction 0-1), elev (ft), TD (deg C, = MWMT - MCMT), EMT (deg C).
This is a Chapman-Richards-style height increment modulated by competition (ccfl, cch) and climate (TD, EMT).

## 3. Mortality (annual SURVIVAL probability) -- mort_parm b0..b4
    eta = b0 + b1*(cr+0.01)^b2 + b3 * (cch>0 ? cch^b4 : 0)
    P_surv_annual = 1 - exp( -exp( clamp(eta, -30, 30) ) )
Inputs: cr (fraction), cch (fraction). Fallback for unfit species: softwood/hardwood median by SPCD<300.

### Independent confirmation of the #75 fix
Greg's own deployed projector computes survival as `ps <- 1 - exp(-exp(eta))` (gomp_surv_annual). This is
EXACTLY the form I corrected gompmort.f90 to in PR #76 (it had the inverted exp(-exp(eta)*FINT)). So the engine
fix is now confirmed three independent ways: the empirical stress test, the 62g documented equation, and Greg's
production projector code. High confidence the fix is correct.

## 4. Climate and competition plumbing the engine needs (for native DG/HG hooks)
- elev (ft): FVS has it (stand ELEV).
- bal, ccfl: FVS computes basal-area-in-larger and crown-competition-factor-in-larger.
- cch: already recomputed natively each cycle by GOMPCCH (the mortality hook), affine-mapped to the gompit
  scale; HG uses cch on a 0-1 fraction, so check the scale alignment.
- EMT (extreme minimum temperature) and TD (= MWMT - MCMT): per-STAND climate scalars from ClimateNA 1991-2020
  normals, looked up by stand lat/lon. This is the only genuinely new engine input; supply per stand via an
  aux table (greg_emt_td_lookup.rds is the projector's per-STAND_CN lookup) or a keyword.

## 5. What this unblocks (#28)
The DG and HG native hooks can now be authored on the GOMPMORT template (common block + LOAD-from-CSV +
evaluator + one-line growth-loop hook in dgf.f90 / htgf.f90), because the exact deployed algebra and the
per-species coefficient tables are both in hand. The CR update in the projector is the fvs-conus kernel
(non-Greg), so Greg's option 1 leaves CR to the native/our kernel (already flagged).

Next steps:
1. Emit dg_parms and hg_parms as flat per-species CSVs (SPCD + B0..B6 / B0..B8), analog to the validated
   62h mortality CSV.
2. Author GOMPDG (dgf.f90) and GOMPHG (htgf.f90) substitution routines from the forms above, env/keyword
   activated, reading those CSVs, plus a per-stand EMT/TD/elev supply.
3. Unit-validate each against the projector's dg_annual / hg_annual on a tree grid (the same approach that
   validated the mortality CSV to 1.6e-7), then a stand A/B on net01.key.
