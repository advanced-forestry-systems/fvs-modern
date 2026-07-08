# What is actually implemented in fvs-modern (audit, 2026-06-15)

Audited the holoros/fvs-modern repo to settle what is in production versus prototype.

## Two config layers, both on main and on the open PR #70 branch

- `config/<variant>.json`: the native FVS coefficients extracted per variant.
- `config/calibrated/<variant>.json`: native coefficients plus a `calibration_multipliers`
  block plus a `self_thinning` slope. This is what the engine calibration reads.
- PR #70 (the species-free integration scaffold) is OPEN, not merged. The calibrated configs,
  however, exist on both main and the branch with the same content, so the calibration framework
  is on main.

## Component calibration factors: framework complete, values partial

All 25 variants have a `calibration_multipliers` block with per-species arrays for five components,
emitted by `config_loader.py` as FVS keywords (MORTMULT, BAIMULT, HTGMULT) or applied via
`set_species_attr`. Coverage by component is uneven, and this is the key point:

| component | keyword | calibrated? |
|---|---|---|
| Mortality | MORTMULT | yes, all 25 variants (real per-species values) |
| Crown ratio | (cr) | yes, all 25 |
| Height-diameter | (htdbh) | yes, all 25 |
| Diameter growth | BAIMULT (dds) | partial: only ~7 variants (acd, ca, cs, kt, ls, nc, on); NE, PN, SN, CR, UT and most others are 1.0 (uncalibrated placeholder) |
| Height growth | HTGMULT (htg) | none: 1.0 for every variant |

So mortality, crown, and height-diameter are fully calibrated across all variants; diameter growth is
calibrated for only about a quarter of variants; height growth is not calibrated anywhere. The 1.0
arrays are placeholders awaiting the per-variant Bayesian fits to be serialized
(`calibration/R/06_posterior_to_json.R`). Some mortality multipliers hit the clip bounds (0.10, 10.0),
worth a look.

## Max SDI: revised values are NOT in production

The production loader emits SDIMAX from `_find_sdi_param(categories)`, i.e. the native (species-weighted
style) per-species values extracted into the config. `config_loader.py` has no reference to brms,
localized, or TreeMap. The revised localized max SDI lives only as:

- `calibration/sdimax/localized_sdimax.py` (and `.R`): a standalone module that computes a per-stand
  localized max SDI from the TreeMap raster / brms table and emits the SDIMAX keyword.
- the benchmark/analysis path (`var_scale_diag.py`), which injects it via `extra_keywords` for testing.

It has not replaced the SDIMAX in the production configs.

## Role of the revised max SDI

It is the density limit that governs long-term self-thinning and density-dependent mortality. Validated
to predict observed self-thinning ~85 percent better than species-weighting, and inside the engine a
level-calibrated localized max SDI matches or beats native in 19 of 20 variants (gains modest,
concentrated in high-bias variants). The important caveat for implementation: mortality is already
calibrated per variant (MORTMULT populated everywhere), and our joint-fit analysis showed the max SDI
level a variant needs inside FVS largely compensates for that mortality calibration. So the revised max
SDI cannot be dropped in under an already-tuned mortality without revalidating; the two must be set
consistently (the unified joint fit, or a co-derived per-variant level alongside MORTMULT).

## Bottom line for the mental model

- Step 1 (FIA calibration modifiers) is real and per-variant for mortality, crown, and height-diameter;
  partial for diameter growth; absent for height growth. Not "fully implemented across all variants."
- Step 2 (revised max SDI) is validated analysis plus a ready module, but is NOT in the production
  configs; adopting it interacts with the already-calibrated mortality.
- Step 3 (species-free equations) are banked bundles plus an injection prototype on the open PR #70;
  not merged, not in production.
