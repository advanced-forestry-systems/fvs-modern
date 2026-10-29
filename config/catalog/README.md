# config/catalog

`catalog_recommendation_by_variant.csv` and `catalog_gates_by_arm.csv` are the fvs-modern configuration catalog, read by `fvs2py.catalog.select_config(variant)`. They were regenerated on October 29, 2026 (second version, after the output red team) by `validation/bakuzis_stress_v2/catalog/` (see its README for the stand sample, the thresholds, the pre-registered rule and the stress-test alternative). Regenerate them only through that pipeline; do not edit by hand. `status` and `recommended` follow the primary rule; `admissible_alt` and `rule_sensitive` show the alternative.
