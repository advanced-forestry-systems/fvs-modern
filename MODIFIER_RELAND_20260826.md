# CONUS modifier layer re-land, 26 August 2026

## Defect

Commit `a5e1e11` (3 July 2026, "Land CONUS modifier layer (6 components) into 25 variant
configs") ran the production landing script `62m_modifier_to_variant_json.R` against summary
CSVs that were stale or otherwise not the paths `MODIFIER_LANDING_MANIFEST.md` names as
authoritative. The defect shipped in `v1.0.0` (tag object pointing to commit `c40e2af`).

A second, larger problem was found during this re-land and is documented separately below:
in the current `v1.0.0` tree, 19 of the 25 variant configs carry **no** `categories_conus_sf`
or `categories_conus_sf_modifier` block at all. Only six variants (`acd`, `ca`, `cs`, `ls`,
`nc`, `on`) retained any CONUS species-free block by the time `release/v1-code-20260811`
reached `c40e2af`. The block was present in all 25 configs immediately after `a5e1e11`; later
commits that regenerated individual variant configs (for example `1d0a124`, "Recover ne
diameter growth species crosswalk") rebuilt the JSON from an older base and dropped every
`categories_conus_*` key, including the modifier block, without anyone noticing. This is a
distinct defect from the wrong-CSV problem and is flagged for Aaron's attention; it is not
fixed by this branch, because fixing it means re-landing the base species-free block (`62c`)
for those 19 variants, which is a separate write outside this task's scope.

## Authoritative source CSVs used (per `MODIFIER_LANDING_MANIFEST.md`)

All paths under `/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/` on Cardinal:

| component | form | summary file(s) |
|---|---|---|
| height crown base (HCB) | traitmed | `hcb/modifier_traitmed/hcb_speciesfree_traitmed_lambda10_global_summary.csv` + `hcb_speciesfree_traitmed_lambda10_gamma_summary.csv`; base gamma `sf_integration/hcb_v2split_sf_gamma.csv` |
| crown recession (CR) | common | `cr/modifier_common_prod/cr_speciesfree_modifier_lambda10_summary.csv` |
| diameter growth (DG) | common | `dg/modifier_common_v8_for_loo/dg_kuehne_v8_modifier_lambda10_summary.csv` |
| height growth (HG) | common | `hg/modifier_common_v5_prod/hg_organon_v5_modifier_lambda10_summary.csv` |
| mortality | common | `mort/modifier_common_prod/mort_speciesfree_plotlevel_modifier_lambda10_summary.csv` |
| ingrowth | common | `ingrowth/modifier_common_prod/ingrowth_v4_plotlevel_modifier_lambda10_summary.csv` |

All six files were confirmed present and read in full on Cardinal before re-landing.

## Landing method

`scripts/62m_modifier_to_variant_json.R` (v2, 2026-06-05, `md5 5d4164a804e498eca7a65618bf39e300`)
from `~/fvs-conus/scripts/` on Cardinal, run against a scratch copy of the 25 `c40e2af` variant
configs, one invocation per component, `--dry_run=FALSE`. R 4.5.2 via
`module purge; module load gcc/12.3.0 R/4.5.2`. Output verified valid JSON and diffed key-by-key
against the working tree to confirm the writer touched only the `categories_conus_sf_modifier`
key and nothing else in any config.

## Before/after coefficient table (representative variant, alphas identical across all 25
variants since the modifier layer is species-free/national, not variant-specific)

Comparison is against `cs.json`, one of the six variants that still had a (wrong) modifier
block at `c40e2af`; the other 19 variants had no block at all (see above), so for them this
change is "absent -> correct" rather than "wrong -> correct."

### Crown recession (CR) — structural defect, now fixed

Shipped form was **trait-mediated**; the corrected, paired-LOO-validated decision is
**common**. This is the confirmed structural error.

| coefficient | shipped (wrong, traitmed) | corrected (common) |
|---|---|---|
| alpha_0 | 0.001038 | 0.001990 |
| alpha_plant | 0.0000214 | -0.004677 |
| alpha_fire | -0.002371 | -0.002189 |
| alpha_insect | -0.000676 | -0.001599 |
| alpha_disease | -0.000782 | -0.002562 |
| alpha_wind | -0.000374 | 0.003411 |
| alpha_harvest | 0.000476 | 0.001945 |
| alpha_cutting | 0.001156 | 0.000889 |
| alpha_siteprep | -0.015229 | -0.007681 |

The old values are not just off, they come from a different-form fit entirely (a
species-trait-weighted model instead of the species-independent common-alpha model), so no
single coefficient in the old table is even estimating the same quantity as its replacement.

### Diameter growth (DG) — numerically wrong, now fixed

| coefficient | shipped (wrong) | corrected | % change |
|---|---|---|---|
| alpha_0 | 0.039565 | 0.050696 | +28.1% |
| alpha_plant | 0.083824 | 0.086573 | +3.3% |
| alpha_fire | -0.037890 | -0.041787 | -10.3% |
| alpha_insect | -0.011425 | -0.019865 | -73.9% |
| alpha_disease | -0.005779 | -0.009277 | -60.5% |
| alpha_wind | 0.055977 | 0.049405 | -11.7% |
| alpha_harvest | 0.015999 | 0.015709 | -1.8% |
| alpha_cutting | 0.219832 | 0.226485 | +3.0% |
| alpha_siteprep | 0.116326 | 0.117158 | +0.7% |

alpha_0 moved 28%, consistent with the audit's flagged "22 percent wrong" (the audit's figure
likely used a different reference statistic; direction and order of magnitude both confirm the
same underlying defect).

### Height growth (HG) — numerically wrong, previously unresolved, now confirmed and fixed

The task named height growth as unresolved. It is now confirmed wrong and corrected: several
disturbance alphas moved by more than 50%, and alpha_fire changed sign-adjacent magnitude by
nearly 3x.

| coefficient | shipped (wrong) | corrected | % change |
|---|---|---|---|
| alpha_0 | -0.012554 | -0.009048 | +27.9% |
| alpha_plant | 0.112054 | 0.168298 | +50.2% |
| alpha_fire | -0.047378 | -0.186230 | -293.1% |
| alpha_insect | -0.079361 | -0.121367 | -52.9% |
| alpha_disease | -0.073522 | 0.001809 | sign flip |
| alpha_wind | -0.105967 | -0.045421 | +57.1% |
| alpha_harvest | 0.066095 | 0.020029 | -69.7% |
| alpha_cutting | -0.023346 | -0.049580 | -112.4% |
| alpha_siteprep | -0.032764 | 0.018835 | sign flip |

### Mortality — previously unresolved, now confirmed CORRECT (no change)

The task named mortality as unresolved. Direct comparison shows the six variants that had a
modifier block already carried the exact authoritative mortality coefficients (0.0% change on
every alpha). Mortality was not actually defective in those files; it was already correct.

### Ingrowth — confirmed CORRECT (no change) in the six variants that had it

Also an exact match, 0.0% change on every alpha, including alpha_plant (-0.200783 in both).
The audit's "sign flip" finding was not reproduced in the currently shipped `cs.json`
ingrowth block; if a sign flip exists it must be in a variant-specific artifact this
re-landing did not touch, since the modifier layer is otherwise species-free/national and
identical across variants. Flagging this discrepancy for Aaron rather than asserting a fix
for a defect that direct inspection did not confirm on the shipped file.

### Height crown base (HCB, trait-mediated) — numerically wrong, now fixed

| coefficient | shipped (wrong) | corrected | % change |
|---|---|---|---|
| alpha_0 | 0.004154 | 0.011433 | +175.2% |
| alpha_plant | -0.162633 | -0.164432 | -1.1% |
| alpha_fire | 0.254582 | 0.277894 | +9.2% |
| alpha_insect | 0.022462 | 0.028893 | +28.6% |
| alpha_disease | -0.004439 | 0.000923 | sign flip |
| alpha_wind | 0.118346 | 0.113462 | -4.1% |
| alpha_harvest | 0.012710 | 0.012929 | +1.7% |
| alpha_cutting | 0.047051 | 0.055731 | +18.4% |
| alpha_siteprep | -0.100786 | -0.089684 | +11.0% |

Form (trait-mediated, trait_mediated_types = plant, insect, harvest, cutting, P_trait = 8) was
already correct; only the coefficient values were wrong.

## Coverage across all 25 variants

Verified programmatically: after this re-land, all 25 variant configs (`acd, ak, bc, bm, ca,
ci, cr, cs, ec, em, ie, kt, ls, nc, ne, oc, on, op, pn, sn, so, tt, ut, wc, ws`) carry a
complete `categories_conus_sf_modifier` block for all six components (`crown_recession`,
`diameter_growth`, `height_growth`, `mortality`, `ingrowth`, `height_crown_base`). Before this
branch, only 6 of 25 had any block, and that block was wrong on crown recession (structural),
diameter growth, height growth, and height crown base.

Confirmed by direct JSON diff that the landing script touched only the
`categories_conus_sf_modifier` key in every file; no other config content changed.

## Caveat: base block still missing for 19 variants

For the 19 variants that lacked `categories_conus_sf` entirely at `c40e2af` (`ak, bc, bm, ci,
cr, ec, em, ie, kt, ne, oc, op, pn, sn, so, tt, ut, wc, ws`), the corrected modifier block this
branch adds has no base species-free equation to multiply against, and per the manifest's own
sequencing note landing the modifier alone onto a config lacking the base block is inert at
runtime. This branch still lands the corrected modifier coefficients for those variants (per
the task's explicit instruction to land all 25), but the fix is functionally inert for those
19 until `62c` (the base species-free landing script) is also re-run for them and the engine
modifier hook is confirmed to read the block. Recommend this become the next tracked item.

## Branch and provenance

Branch `fix/modifier-relanding-20260826`, based on `release/v1-code-20260811` tip (`c40e2af`,
same commit as the `v1.0.0` tag). `release/v1-code-20260811`, `main`, and the `v1.0.0` tag were
not modified, moved, or merged into. This branch was not merged anywhere; it is staged for
Aaron's review.

## Addendum: base species-free block (`categories_conus_sf`) restored for 19 variants

Follow-up to the paragraph above. Two corrections to the record and the fix itself:

**Correction 1 -- `a5e1e11` never wrote `categories_conus_sf`.** `git show a5e1e11 -- config/calibrated/ak.json`
shows the commit's only change is appending `categories_conus_sf_modifier` after the existing
`_emit_sdimax` key; `categories_conus_sf` is absent both before and after `a5e1e11` for `ak.json`
(checked directly against `a5e1e11^`). `a5e1e11` landed the modifier layer only, via
`62m_modifier_to_variant_json.R`; it is not the base-block landing commit and is not usable as
structural ground truth for `categories_conus_sf`.

**Correction 2 -- the true base-block landing history.** `categories_conus_sf` was first landed by
`c1af06d` ("Land categories_conus_sf (Leg B species-free) across 25 variants, 6 components (#83)",
2 Jul 2026), which ran `62b_speciesfree_to_variant_json.R` (present on Cardinal at
`calibration/R/62b_speciesfree_to_variant_json.R` and `sf_integration_dev/62b_speciesfree_to_variant_json.R`,
md5 `924e2c0a75b43e0a1082b4c5f1e90d32` for both copies). This is the script the manifest and prior
audit call "62c" -- Cardinal's `sf_integration_dev/62c_arms_to_variant_json.R` is a different,
later (2 Jul, same day) script that lands the unrelated `categories_conus_mod` / `categories_conus_organon`
arms, not the species-free base block. `586b3cc` ("Gate hybrid source map on per-species reliability,
shrinkage w>=0.5") re-ran the same landing shortly after and is reflected in every surviving block's
`metadata.default_policy` string. The block then rode unmodified through `4d276d8`, `8a5c3a5`, `a5e1e11`
(modifier only), `54bb3b0`, `c2f30c0`, `10c2f82`, `e46f60a`, `2880573`, and `adc372c` ("Reconcile
three-arm integration into main (union)", the commit `release/v1-code-20260811` and `fix/modifier-relanding-20260826`
both descend from). For all 19 affected variants, `adc372c` is the last commit in that variant's own
file history where `categories_conus_sf` is present, confirmed individually per variant by walking
each file's `git log` and testing every commit for the key's presence. Immediately after `adc372c`,
each variant's config was rebuilt by its own per-variant "recover / regenerate" commit (`dd088af` for
`ak`, and an analogous commit per remaining variant) that reconstructed the JSON from an older base
and silently dropped every `categories_conus_*` key gained since, including `categories_conus_sf`. That
regeneration pattern, not `a5e1e11`, is the actual defect that produced the 19-variant gap.

**Fix.** For each of the 19 variants (`ak`, `bc`, `bm`, `ci`, `cr`, `ec`, `em`, `ie`, `kt`, `ne`, `oc`,
`op`, `pn`, `sn`, `so`, `tt`, `ut`, `wc`, `ws`), `categories_conus_sf` was extracted verbatim from that
variant's own `adc372c` blob (`git show adc372c:config/calibrated/<v>.json`) and spliced into the
current (post `c71564f`) config as a pure text insertion: parse to locate the file's trailing `}`,
insert `"categories_conus_sf": {...}` as a new key immediately before it, re-serialize nothing else.
Verified by structural diff that every other top-level key, including the just-corrected
`categories_conus_sf_modifier`, is byte-identical before and after (`git diff --stat` shows insertions
only, zero deletions, across all 19 files). The restored block carries the same schema and
`metadata.default_policy` (shrinkage-gated hybrid source map) as the six variants that never lost it,
confirming `acd`/`ca`/`cs`/`ls`/`nc`/`on` and the 19 restored variants now share one lineage for this
block. Species-specific content (trait tables, per-species hybrid source maps, RE tables) is unique
per variant, pulled from that variant's own `adc372c` state, not copied across variants.

All 25 configs now carry structurally complete `categories_conus_sf` (diameter_growth, height_growth,
height_diameter, height_crown_base, crown_recession, mortality, metadata) and
`categories_conus_sf_modifier` (diameter_growth, height_growth, mortality, ingrowth, height_crown_base,
crown_recession) blocks, so the corrected modifier coefficients from the section above are no longer
orphaned in any variant. This restoration does not touch modifier coefficient values; it only restores
the base scaffolding those coefficients depend on. The still-open item from the section above (add the
engine modifier hook that reads `categories_conus_sf_modifier` and applies it against
`categories_conus_sf`, and validate held-out disturbed/treated plots) remains open and is out of scope
for this commit.

`release/v1-code-20260811`, `main`, and the `v1.0.0` tag are untouched. This restoration lands as a new
commit on `fix/modifier-relanding-20260826`, on top of `c71564f`, per the same staged-for-review
convention (not merged).
