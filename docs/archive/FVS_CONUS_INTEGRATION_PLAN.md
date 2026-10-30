# Plan for incorporating fvs-conus parameters into fvs-modern

*Drafted 2026-05-09 from a Cardinal-side audit of `~/fvs-modern/calibration/output/conus/`. The plan assumes the CONUS Phase 4 fits remain in active development and final model selection is pending.*

## Current state on Cardinal

The fvs-conus calibration lives at `~/fvs-modern/calibration/output/conus/` (56 GB total). All six target components have main posterior fits in place plus modifier fits at three regularization levels.

| Component | Main model | Latest fit timestamp | Largest artifact |
|---|---|---|---|
| Diameter growth | dg_kuehne, dg_organon | 2026-05-09 (residuals) | 13 GB `dg_kuehne_cspi_traits1_fit.rds` |
| Height growth | hg_organon_fixedK | 2026-04-26 | 4.9 GB `hg_organon_fixedK_cspi_traits1_fit.rds` |
| Height-DBH | htdbh_wykoff_lognormal | 2026-04-29 | 7.3 GB fit |
| Height-crown base | hcb_organon | 2026-05-03 | 8.5 GB fit |
| Mortality | mort_logit_simple, mort_gompit | 2026-05-02 | 14 GB fit |
| Crown recession | cr_recession | 2026-04-29 | 8.3 GB fit |

Each component also has modifier fits at lambda 5, 10, 20 (regularization sweep). The infrastructure is structurally complete; what remains "under construction" is final model selection per component (Kuehne vs Organon for DG, gompit vs logit_simple for mortality, lambda choice for modifiers, formula variant decisions).

Already-extracted parameter summaries at the top of `calibration/output/conus/`:

- `dg_organon_fixed_draws.csv` (4,000 rows × 16 cols): per-draw fixed-effect parameters (mu_b0, b1..b10, K1, K2, sigma_sp, sigma_eco, sigma)
- `dg_organon_species_intercepts.csv` (104 rows): SPCD-keyed random intercepts
- `dg_organon_ecodiv_intercepts.csv` (151 rows): ecodivision-keyed random intercepts

## Schema gap to close

The existing `config/calibrated/{variant}.json` schema is shallow:

```json
"categories": {
  "growth": {"B1": [...], "B2": [...]},
  "height_diameter": {"HD_A": [...], "HD_B": [...], "HD_C": [...]},
  "crown": {"BCR1": [...], "BCR2": [...], "BCR3": [...], "BCR4": [...]},
  "site_index": {"SDICON": ...},
  "self_thinning": {"CALIBRATED_SLOPE": ..., "REINEKE_STANDARD": ...}
}
```

Two coefficients per species for diameter growth, three for height-DBH, four for crown. The CONUS posteriors are richer: 16 fixed-effect parameters plus species and ecodivision random intercepts plus posterior draws plus modifier coefficients.

Three options for closing the gap, from least to most invasive:

**Option A: collapse to existing schema.** Compute MAP point estimates per species per variant from the CONUS posterior, fit them onto the existing B1/B2 keys. Loses uncertainty quantification and the ecodivision dimension. Cheapest to integrate but the worst science.

**Option B: alternate parameter set.** Add a parallel `categories_conus` block in each `{variant}.json` keyed by component name with the full CONUS parameter shape (fixed effects + species RE table + ecodiv RE table + draws ID for posterior sampling). Existing readers continue to use `categories`; a new flag in `config_loader.py` selects which. Preserves backwards compatibility. Recommended.

**Option C: separate config tree.** Create `config/conus/{variant}.json` with the new schema, leave `config/calibrated/{variant}.json` as the existing per-variant fits. Cleanest separation but doubles the config surface.

Option B fits best because the runtime layer (`config_loader.py`) already has versioned parameter switching (`default | calibrated | custom`); adding `conus` is one more switch. The existing per-variant fits stay as a fallback, and the CONUS fits land as a richer alternative.

## Proposed integration sequence

### Phase 1: extraction infrastructure (do once now)

1. Generalize the DG-only summary extractor so it handles all six components. Add `calibration/R/61_extract_conus_summaries.R` that reads each component's `*_fit.rds` and writes:
   - `{component}_fixed_draws.csv` — full posterior of fixed effects
   - `{component}_species_intercepts.csv` — SPCD-keyed RE
   - `{component}_ecodiv_intercepts.csv` — ecodivision-keyed RE
   - `{component}_modifier_summary.csv` — per-lambda modifier means
2. Generalize the modifier comparison logic so the chosen lambda for each component is captured in metadata (not picked by hand each time).
3. Run on one component to verify (DG is already done; HG is the next-most-mature).

### Phase 2: schema extension (do when first component is ready)

1. Define the `categories_conus` block schema in `config_loader.py` and `config/uncertainty.py`.
2. Add a `parameter_source` enum field to track which calibration produced each block (`per_variant` for legacy, `conus_phase4` for the new fits).
3. Write `calibration/R/62_conus_to_variant_json.R` that, given a component and a variant:
   - Looks up the variant's species set (`species_definitions.FIAJSP`)
   - Selects matching SPCD intercepts from the CONUS species_intercepts CSV
   - Weights ecodivision intercepts by the variant's geographic coverage (or includes all and lets runtime select)
   - Writes a `categories_conus.{component}` block in `config/calibrated/{variant}.json`

### Phase 3: per-component land-as-they-finish (the rolling integration)

When each component's analysis lands, it goes through:

1. Run `61_extract_conus_summaries.R --component {name}` on Cardinal to dump CSVs.
2. Run `62_conus_to_variant_json.R --component {name}` per variant to land the parameters into per-variant configs.
3. Run a per-variant validation (basal area RMSE, R-squared) against the FIA benchmark to confirm the CONUS-derived parameters perform at least as well as the per-variant fits.
4. Update `MEMORY.md` and `CALIBRATION.md` to mark the component's status `production`.
5. Add the component to `components_updated` in the `categories_conus.calibration` metadata block.

### Phase 4: runtime switch (do when at least 4 of 6 components are integrated)

1. Wire `config_loader.py version="conus"` to read `categories_conus` instead of `categories`.
2. Add a fallback rule: components missing from `categories_conus` fall back to `categories`. This lets a partial integration run with mixed sources.
3. Add a `version="hybrid"` mode that defaults to `conus` per component but allows per-component override (e.g., DG from CONUS, mortality still from per-variant).
4. Smoke-test `version="conus"` against the FIA benchmark across the 19 testable variants.

### Phase 5: deprecation (do when 6 of 6 components are integrated and validated)

1. Mark the per-variant fits as `legacy` in metadata.
2. Update `MANUSCRIPT_REVISION_PLAN.md` to integrate CONUS-vs-per-variant comparison results.
3. Decide whether to keep the per-variant configs as a backup or remove them.

## Integration triggers

The natural triggers for each phase are:

- **Trigger Phase 1 immediately**: extraction code is reusable and won't be wasted even if final model selection changes.
- **Trigger Phase 2 when**: at least one component is finalized (HG looks closest given its 2026-04-26 fit timestamp and "HG complete" note in MEMORY.md).
- **Trigger Phase 3 per component**: when its model selection is final and modifier lambda is chosen.
- **Trigger Phase 4 when**: 4 of 6 components have landed and you trust the runtime fallback rule.
- **Trigger Phase 5 when**: 6 of 6 + validation passes.

## What to do this week

Three concrete actions while waiting for analysis to finish:

1. **Run Phase 1 step 1** for HG (since HG is closest to done). Produces `hg_organon_fixedK_fixed_draws.csv`, `hg_organon_fixedK_species_intercepts.csv`, `hg_organon_fixedK_ecodiv_intercepts.csv`. About 30 minutes of R work on Cardinal — one Rscript invocation per existing fit.

2. **Draft the schema extension** in `config_loader.py` and `config/uncertainty.py` against the HG output. This makes the schema concrete and shakes out any field-name issues before applying to the other five components.

3. **Write `61_extract_conus_summaries.R` and `62_conus_to_variant_json.R`** as templates with HG-specific implementations and TODO stubs for the other five. Each component fills in its own conversion logic (parameter naming, link function, etc.) in its own block of the script.

The total work outside Cardinal compute time is about 1 to 2 days of R/Python plus a half-day for the runtime switch in `config_loader.py`. The Cardinal compute cost is whatever extraction takes per component (modest — reading existing RDS files and writing CSVs).

## Risk register

| Risk | Mitigation |
|---|---|
| Schema lock-in: pick wrong shape for `categories_conus` and force later migration | Draft the schema with HG only, validate against one variant's needs before extending to all 6 components |
| Performance: 4,000 posterior draws × 25 variants × 6 components could blow up runtime memory | Store draws per component once globally; variants reference by index |
| Variant-CONUS coverage gaps: some species in a variant may not have CONUS coefficients (out-of-range) | Define an explicit fallback rule (use nearest species, or fall back to per-variant fit) |
| Ecodivision weighting: each variant covers multiple ecodivisions in different proportions | Compute per-variant ecodivision weight tables once at integration time, store in `categories_conus.ecodiv_weights` |
| Backwards incompatibility: existing scripts that read `categories.growth` directly will break | Keep `categories` populated as the legacy view; only ADD `categories_conus` |
| Modifier lambda selection: 3 lambdas per component is a hyperparameter you'll need to commit to | Pick the one with best held-out FIA RMSE; store the choice in component metadata |

## Outputs of this plan

When Phase 5 completes, the deliverables on origin/main will be:

- New R scripts: `61_extract_conus_summaries.R`, `62_conus_to_variant_json.R`
- Extended schema: `categories_conus` block in all 25 `config/calibrated/{variant}.json`
- Updated runtime: `config_loader.py` supports `version="conus"` and `version="hybrid"`
- Updated docs: CALIBRATION.md, MEMORY.md, MANUSCRIPT_REVISION_PLAN.md
- A v2026.06 calendar tag carrying the CONUS integration

Most of this can land incrementally. None of it requires waiting for all six components to finish at once.
