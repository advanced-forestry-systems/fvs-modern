"""
Patch DG v8, HG v5, ingrowth v4 drivers to auto-detect traits_v3 and use the
9-trait set (climate_exposure + low_adaptive_cap replacing vulnerability_score
+ adding LAC). Behavior with v1/v2 traits unchanged.

Detection rule: if the loaded traits table has both 'climate_exposure' and
'low_adaptive_cap' columns, use the v3 trait set; otherwise use the v2 set.
This is robust to any path convention.
"""

paths = [
    "/users/PUOM0008/crsfaaron/fvs-modern/calibration/R/32_fit_dg_kuehne_v8.R",
    "/users/PUOM0008/crsfaaron/fvs-modern/calibration/R/32_fit_hg_speciesfree_v5.R",
    "/users/PUOM0008/crsfaaron/fvs-modern/calibration/R/35_fit_ingrowth_negbinom_v4.R",
    "/users/PUOM0008/crsfaaron/fvs-modern/calibration/R/35b_fit_ingrowth_hurdle.R",
    "/users/PUOM0008/crsfaaron/fvs-modern/calibration/R/36_fit_ingrowth_species_composition.R",
]

OLD = (
    'trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",\n'
    '                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",\n'
    '                "vulnerability_score", "sensitivity")'
)

NEW = (
    '# Auto-detect traits_v3 (decomposed Potter VCC: CE + S + LAC) vs v2 (composite vuln_score + S)\n'
    'use_v3_traits <- all(c("climate_exposure", "low_adaptive_cap") %in% names(traits))\n'
    'if (use_v3_traits) {\n'
    '  cat("[traits] detected v3 layout: using decomposed Potter components (CE+S+LAC)\\n")\n'
    '  trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",\n'
    '                  "leaf_longevity_months", "max_ht_m", "max_dbh_cm",\n'
    '                  "climate_exposure", "sensitivity", "low_adaptive_cap")\n'
    '} else {\n'
    '  trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",\n'
    '                  "leaf_longevity_months", "max_ht_m", "max_dbh_cm",\n'
    '                  "vulnerability_score", "sensitivity")\n'
    '}'
)

for p in paths:
    s = open(p).read()
    if OLD in s:
        s = s.replace(OLD, NEW, 1)
        open(p, "w").write(s)
        print(f"patched: {p}")
    elif "use_v3_traits <- all" in s:
        print(f"already patched: {p}")
    else:
        # Try a more flexible needle: indent could differ
        print(f"NEEDLE NOT FOUND in {p}")
