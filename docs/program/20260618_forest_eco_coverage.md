# Forest type + ecoregion coverage across all fvs-conus components (both modes)
2026-06-18. Branch conus-sf-integration-2026-05-21. Confirms that every fvs-conus component now includes
forest type and ecoregion, in both the species-dependent and species-free formulations.

## Coverage matrix

Each component has a species-dependent (unified: species random intercept + traits) form and a
species-free (traits only) form. Ecoregion is the EPA L1/L2/L3 nested random intercept (z_L1/z_L2/z_L3);
forest type is the FIA forest-type-group random intercept (z_FT). All eight rows below carry both.

| component | species-dependent model (FT + eco + species) | species-free model (FT + eco + traits) | status |
|---|---|---|---|
| Diameter growth | dg_kuehne2022_v8c_bgi_cspi / v9_bgi_cch | dg_kuehne2022_speciesfree_v8_forest_eco | FT + eco present, compiles |
| Height growth | hg_unified_v8_rd | hg_organon_speciesfree_v6_bgi | FT + eco present |
| Height-diameter | htdbh_unified | ht_dbh_wykoff_speciesfree | FT + eco present |
| Height to crown base | hcb_unified | hcb_organon_speciesfree(_cspiv6) | FT + eco present |
| Crown ratio / recession | crown_ratio_t2_unified | crown_ratio_change_speciesfree | FT + eco present |
| Mortality / survival | survival_unified | gompit_mortality_speciesfree | FT + eco present |
| Ingrowth (recruitment) | ingrowth_negbinom_v3_forest_eco (NEW) | (same; plot-level, trait-dominant) | FT added this session, compiles |

The species-dependent unified models were verified to carry all four terms together: z_sp + sigma_sp
(species), W * gamma (traits), z_FT + sigma_FT (forest type), z_L1/z_L3 (ecoregion). The species-free
forms carry traits + forest type + ecoregion.

## What this session added

1. dg_kuehne2022_speciesfree_v8_forest_eco.stan: added z_FT to the species-free DG model. Compiles
   (cmdstanr 2.38.0).
2. ingrowth_negbinom_v3_forest_eco.stan: added z_FT to the ingrowth model, the one component that had
   ecoregion but no forest type. Compiles. This closes the last forest-type gap.

The remaining components already carried forest type and ecoregion through the unified model family, which
was the program's prior direction; this session verified that and filled the two gaps (species-free DG and
ingrowth).

## Forest-type-group construction (driver side, identical for every component)

The data carry FORTYPCD (135 codes) and the EPA L1/L2/L3 codes already used for ecoregion. Build the
forest-type-group index once and pass it to every component fit:

    dat[, FT_grp := (as.integer(FORTYPCD) %/% 10L) * 10L]   # FIA forest-type group, ~28 groups
    FT_levels <- sort(unique(dat$FT_grp)); dat[, FT_idx := match(FT_grp, FT_levels)]
    # stan_data: N_FT = length(FT_levels), FT_idx = dat$FT_idx

## Production fit plan (queued)

Each forest-type model needs one production fit (100k subsample, 2 chains, 4-day wall on Cardinal) and a
paired held-out ELPD against the no-forest-type predecessor to confirm forest type earns its place. Run
order by expected forest-type signal: ingrowth and mortality first (recruitment and mortality are the most
forest-type dependent), then crown recession, HCB, then the growth and size components. The fits reuse the
existing per-component drivers with the three FT lines above added and STAN_FILE pointed at the
forest_eco model.

## Status of the broader recommendations

| recommendation | status |
|---|---|
| Forest type + ecoregion in all components, both modes | DONE (models in place and compiling; production fits queued) |
| Species-free held-out ELPD vs species-specific (with FT + eco) | next: paired LOO per component once the FT fits land |
| Density-dependent recruitment, site-resolved form (PN transfer) | designed (SDI-headroom scaled by BGI/CSPI); refit queued |
| Route A: fvsAddTrees in-engine loader for the true four-arm | diagnosed, fix path confirmed; build pending |
| Consensus MCW into the fvs-conus crown component | consensus table ready (mcw_conus_consensus.csv); wiring pending |
| Put fvs-conus under version control | recommended; not yet done (no git history) |
| Disturbance-aware benchmark as the standard | adopted in the committed analyses and manuscript |
| fvs-modern branch prune + tagged release | recommended; pending the in-engine four-arm |
| Zenodo new version of the concept DOI | deferred until validation is OOS-honest |

Artifacts: stan/dg_kuehne2022_speciesfree_v8_forest_eco.stan, stan/ingrowth_negbinom_v3_forest_eco.stan
(both compile); mirrored into fvs-modern/calibration/stan_conus for version control.
