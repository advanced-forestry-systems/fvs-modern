# fvs-conus: forest-type + ecoregion structure, and cross-repo refinements
2026-06-18. Branch conus-sf-integration-2026-05-21. Implements the request that the fvs-conus equations
depend on forest type and ecoregion, and records the refinement recommendations from this session's full
assessment of fvs-modern and fvs-conus.

## 1. Forest type and ecoregion in the species-free components

Finding: ecoregion is already in. The species-free DG model (v6_data_driven) already carries EPA ecoregion
at three nested levels as non-centered random intercepts, z_L1, z_L2, z_L3, plus a level-1 site slope
(z_L1_csi). Forest type was NOT present. The genuine addition is a FIA forest-type-group random intercept.

Implemented and verified: dg_kuehne2022_speciesfree_v8_forest_eco.stan adds a forest-type-group level
z_FT mirroring the ecoregion pattern, so the equation now depends on both forest type and ecoregion. The
model compiles cleanly (cmdstanr, gcc 12.3.0); the new parameters z_FT_raw and sigma_FT are present. The
five additions over v6:

- data: int N_FT; array[N_obs] int<lower=1,upper=N_FT> FT_idx;
- parameters: vector[N_FT] z_FT_raw; real<lower=0> sigma_FT;
- transformed parameters: vector[N_FT] z_FT = sigma_FT * z_FT_raw;
- model: z_FT_raw ~ std_normal(); sigma_FT ~ normal(0, 0.3);
- linear predictor (model and generated quantities): eta += z_FT[FT_idx]

The data already carries the needed column. conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds (1,174,749
rows) has FORTYPCD (135 codes), plus EPA_L1_CODE, EPA_L2_CODE, EPA_L3_CODE already used by the ecoregion
levels. FORTYPCD at 135 levels is too granular for a stable random effect; bin to the FIA forest-type
GROUP (roughly FORTYPCD rounded down to the nearest 10, about 28 groups), which is the right grain for a
partially pooled intercept.

## 2. Driver wiring (the remaining step to fit)

Extend the species-free fit driver (R/32_fit_dg_kuehne_speciesfree.R) at the index-build and data-list
steps, exactly mirroring the L1/L2/L3 pattern:

    dat[, FT_grp := (as.integer(FORTYPCD) %/% 10L) * 10L]   # FIA forest-type group
    FT_levels <- sort(unique(dat$FT_grp))
    dat[, FT_idx := match(FT_grp, FT_levels)]
    # in stan_data: N_FT = length(FT_levels), FT_idx = dat$FT_idx

Caveat to resolve first: the pulled species-free driver builds its data list with ln_csi, ba_x_rd,
bal_x_rd, while v6/v8 expect ln_sicond, ccfl1, is_plantation, ln_elev, sdi_complexity, rd_additive. The
driver and the Stan model versions must be matched before fitting; use the driver that actually produced
the v6 data-driven fit (its data list names match v6), then add the three FT lines above and point
STAN_FILE at v8. A subsample pilot (SUBSAMPLE = 20000) should sample in a few minutes and confirm sigma_FT
is identified before the full 100k production run.

## 3. Roll the FT level out to the other components

The same z_FT block applies unchanged to the HG, HCB, CR, mortality, and HT-DBH species-free models (they
share the ecoregion-level structure). Add z_FT to each, refit, and read sigma_FT to see where forest type
carries signal. Expectation from the species stress test: forest type should matter most where stand
context drives the rate, plausibly mortality, crown recession, and HCB.

## 4. Tie to the species stress test (this session)

The stress test showed the per-species term is well identified for mortality, diameter growth, and crown
recession, and weak for height growth, height-diameter, and HCB. Forest type and ecoregion are the natural
substitute structure for the components where species identity is weak: a forest-type-group plus ecoregion
intercept can absorb the among-stand variation that a noisy species term was straining to fit, and it
transfers to the roughly 232 FIA species without a fitted curve. The clean experiment is a paired
held-out ELPD comparison per component: species-specific, versus species-free with forest type and
ecoregion. That single comparison decides the production structure and supersedes the earlier
species-free-only confirmatory check.

## 5. Cross-repo refinement recommendations (from the full assessment)

fvs-conus:
- Make forest-type-group plus ecoregion the standard grouping for every species-free component (Section 3),
  then run the paired held-out ELPD to confirm it matches or beats the species-specific fits.
- Density-dependent recruitment needs a site-resolved form so Pacific Northwest TPH transfers out of
  sample; refit with the SDI-headroom rate scaled by site productivity (BGI or CSPI) and re-run the
  held-out validation.
- Put fvs-conus under version control with a data-aware gitignore; it is a 4 GB working directory with no
  git history, which is the main reproducibility gap for the manuscript and the Zenodo deposit.
- Feed the CONUS-consistent consensus MCW (mcw_conus_consensus.csv) into the fvs-conus crown component so
  crown width is harmonized across variants.

fvs-modern:
- Finish Route A by loading stands via fvsAddTrees in memory (bypasses the in-process database read that
  causes the extree.f90 segfault), which unblocks the true in-engine arms C and D and a fair four-arm.
- Adopt the disturbance-aware benchmark as the standard evaluation: stratify by FIA TRTCD and DSTRBCD, or
  simulate recorded removals, before reporting bias. This is the one change that most affects the headline.
- Prune the branch sprawl and merge conus-sf-integration into main; cut a tagged release once the
  in-engine four-arm lands.
- Keep the framing honest: disturbance-aware benchmark plus prototype adjustment layer, size levers
  validated out of sample, recruitment a prototype.

Shared:
- The single highest-value build remains the true single-framework four-arm (Route A via fvsAddTrees); it
  is the shared headline figure and resolves the projector-versus-engine sign issue for good.
- Defer the Zenodo deposit until the validation is out-of-sample-honest, then deposit as a new version of
  the fvs_perseus_conus concept DOI and backfill the DOI into both manuscripts.

## Artifacts

- stan/dg_kuehne2022_speciesfree_v8_forest_eco.stan (compiles; forest type + ecoregion).
- This note (wiring spec + refinements). The production fit is a driver-wiring step, not a model gap.
