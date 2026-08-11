# Adopting the CONUS equations natively into fvs-modern: assessment and plan (2026-06-28)

Question: what can we do on fvs-modern to prepare adopting the new CONUS equations into the FVS framework?
Answer: more than expected, and without waiting on the blocked Python stop-point hook (engine bug 2). The
engine already has a working native-substitution pattern. This turn I generated the one missing file that
makes Greg's mortality run natively today, and below is the component-by-component path for the rest.

## 1. The key finding: a native substitution mechanism already exists (and bypasses bug 2)

fvs-modern is the standard FVS source tree: each geographic variant under src-converted/<variant>/ carries the
native Fortran growth subroutines (dgf diameter growth, htgf height growth, crown, morts mortality, regent
regeneration/ingrowth), and there is a shared base/. The NE variant (src-converted/ne/, where conus_greg
routes) is fully present.

Inside base/gompmort.f90 there is already a complete, integrated native substitution of a CONUS gompit
mortality, the GOMPMORT system:
- GOMPMC.f90 (common): per-species coefficient table GB(MAXSP,5)=b0..b4, GHAVE flag (gompit owns a species
  only if it carries a fitted row), CCHT per-tree crown closure.
- GOMPLOAD: reads the coefficient CSV (header line + rows SPCD n b0 b1 b2 b3 b4), maps FIA SPCD onto FVS
  species via FIAJSP.
- GOMPSURV(ISPC,CR,CCH,FINT,SURV): eta = b0 + b1*(cr+0.01)^b2 + b3*cch^b4; H=exp(eta); S=exp(-H*FINT).
- GOMPCCH: recomputes cch each cycle from the live treelist via an ORGANON crown-closure port, affine-mapped
  to the gompit cch scale (CCH_A=0.062, CCH_B=0.0036, Spearman 0.93).
- Hooked into morts.f90 in every eastern variant: `IF(LGOMP) CALL GOMPCCH` and
  `IF(LGOMP.AND.GHAVE(ISPC)) ... GOMPSURV ...` substitute the gompit survival INSIDE the growth loop, so
  density and growth interact with the substituted mortality each cycle.
- Activated with NO recompile: env var FVS_GOMPIT=1 plus FVS_GOMPIT_COEF=<csv>, or a GOMPMORT keyword.

This is the production adoption pattern: env-var/keyword-activated Fortran substitution subroutines that read a
flat coefficient file and evaluate our form inside the growth loop. It does NOT use the dlopen stop-point hook
(sf_shadow_injector.py) that is blocked by bug 2. Bug 2 is therefore demoted from critical path to "nice to
have generic alternative." The native route is the real one and it already works for mortality.

## 2. What I did this turn: Greg mortality goes native

The GOMPMORT machinery was complete but the coefficient file it reads did not exist on disk. I wrote
calibration/R/62h_greg_mortality_csv.R, which emits greg_mortality_coefficients.csv from Greg's
mort_parm_base_rate_cr_cch.RDS in the exact GOMPLOAD format (header + SPCD,n,b0..b4), de-duplicated to one
best-fit row per species by minimum nll. Result: 92 unique FIA species,
~/fvs-modern/config/greg_mortality_coefficients.csv.

Greg's CONUS mortality (option 1, the cch-gompit form) is now natively runnable:
    export FVS_GOMPIT=1
    export FVS_GOMPIT_COEF=/users/PUOM0008/crsfaaron/fvs-modern/config/greg_mortality_coefficients.csv
This is the first of the new equations to become a first-class native FVS equation rather than a projector
overlay. Next concrete step here: an A/B run (native default mortality vs FVS_GOMPIT) on an NE stand set to
confirm the substitution moves density as expected.

## 3. Component-by-component adoption status

Mechanism legend: NATIVE-DONE (substitution exists + integrated), TEMPLATE (mirror the GOMPMORT pattern in
the named subroutine), KEYWORD (config_loader multiplier path), EXTEND (native routine exists but our form is
richer than what it currently evaluates).

| component | our form | native subroutine | adoption path | status |
|---|---|---|---|---|
| Mortality (Greg, opt 1) | cch-gompit b0..b4 | base/gompmort.f90 GOMPSURV, hook in morts.f90 | supply coeff CSV | NATIVE-DONE; CSV emitted this turn |
| Mortality (ours, opt 2/3) | gompit b0..b7b + W.gamma + ecoregion/FT/species REs + std cr_z,cch_z | base/gompmort.f90 | EXTEND GB to N cols + add trait/RE/standardization terms to GOMPSURV | spec below; coeffs ready (z_sp now exported) |
| Diameter growth | Kuehne (opt2) / Greg B0..B7 (opt1) / species-free (opt3) | <variant>/dgf.f90 (+ dgdriv.f90) | TEMPLATE: GOMPDG-style hook reading a DG coeff CSV | not started |
| Height growth | ORGANON (opt2) / Greg B0..B8 (opt1) | <variant>/htgf.f90 | TEMPLATE | not started |
| HCB / crown | species-free v2split | <variant>/crown.f90, cratet.f90 | TEMPLATE | not started |
| HT-DBH | Chapman-Richards | <variant>/dubscr.f90 (bark/dbh-ht) | TEMPLATE; new asymptotic branch | pending the 50k refit verdict |
| Ingrowth | hurdle count + multinomial composition | <variant>/regent.f90 (essubh.f90) | TEMPLATE | not started |
| Stand SDImax / BA cap (opt2/3 Garcia) | quantile SDImax + monomolecular BA | SDIMAX keyword + a carrying-capacity hook | KEYWORD (fix field order, re-enable) + TEMPLATE | SDIMAX disabled 2026-06-16, needs fix |

Options map cleanly onto this: option 1 (Greg) needs native DG + HG hooks plus the mortality CSV (done) and no
Garcia layer; options 2 and 3 share the same native routines and differ only in which coefficient file is
loaded (species-dependent b2 vs species-free b1), exactly like the config JSON blocks already differ.

## 4. Prep we can do now (no fits, no bug 2)

1. DONE: greg_mortality_coefficients.csv emitter, file generated. Run the NE A/B to validate.
2. Coefficient-CSV emitters for the remaining components and all three options, in the flat per-species format
   the native hooks read. config_loader already decomposes categories_conus / _sf / _greg into runtime form;
   add a thin "export native CSV" method beside generate_keywords. One emitter per component (DG, HG, HCB,
   HT-DBH, ingrowth) times three options, but they are mechanical once the schema is set.
3. Author the TEMPLATE hooks by copying the GOMPMORT pattern: a common block, a LOAD routine (env var +
   keyword), an evaluator subroutine, and a one-line substitution in the host subroutine (dgf/htgf/crown/
   regent), each guarded by GHAVE-style per-species coverage. The mortality trio is the worked example.
4. EXTEND GOMPSURV for our richer gompit (options 2/3): widen GB beyond 5 columns to carry b1..b7b plus the
   per-species trait effect (precompute W.gamma per species so the engine stores one number), the
   ecoregion/forest-type/species REs (also collapsible to one per-species + per-stand-ecoregion addend), and
   the cr_z / cch_z standardization constants. Because trait and RE effects can be pre-reduced to per-species
   and per-ecoregion scalars offline, the in-engine form stays close to the current one. Spec this as the
   GOMPMC v2 layout before coding.
5. Re-enable SDIMAX with the corrected field order (the 2026-06-16 WO-1 disable was a wrong-field-order bug),
   then add the monomolecular BA carrying-capacity term for the Garcia layer used by options 2/3.
6. Keep engine bug 2 as a maintainer ticket for the generic dlopen stop-point path, but it is no longer
   gating: the native substitution route above is the production path.

## 5. Recommended order
1. Validate native Greg mortality (NE A/B with FVS_GOMPIT) now that the CSV exists.
2. Write the GOMPMC v2 spec (extended coefficient layout) for our richer gompit, since the survival z_sp export
   just landed and the coefficients are ready.
3. Add the DG and HG native hooks (GOMPMORT pattern) so option 1 (Greg) is fully native end to end; that is the
   smallest closed loop and a clean A/B against the projector.
4. Generalize the coefficient-CSV export in config_loader for all components and options.
5. Re-enable and fix SDIMAX, then add the BA carrying-capacity hook for options 2/3.

## 6. Artifacts this turn
- ~/fvs-conus/R/62h_greg_mortality_csv.R (emitter)
- ~/fvs-modern/config/greg_mortality_coefficients.csv (92 species, GOMPLOAD format)
Both on Cardinal; should be committed (the emitter to fvs-conus-components, the CSV to fvs-modern alongside
the conus_greg config). Not yet committed.
