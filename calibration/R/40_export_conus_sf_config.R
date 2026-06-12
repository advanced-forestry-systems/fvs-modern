#!/usr/bin/env Rscript
# 40_export_conus_sf_config.R
# Companion to 06_posterior_to_json.R for the unified CONUS species-free variant.
# Reads the six fvs-conus component bundles + Layer 2 modifier JSONs + CSPI metadata
# and writes config/calibrated/conus_sf.json (schema: calibration/schema/conus_sf_config.schema.json)
# plus thinned posterior draws for uncertainty propagation.
#
# Emits BOTH legs per component:
#   species_independent : trait term W %*% gamma  (Leg B, every species)
#   species_dependent   : fitted per-species b0   (Leg A, well-sampled species)
#   blend               : shrinkage weight w_s = n_s/(n_s+kappa)
#
# STATUS: scaffold. Fill the bundle paths (TODO) once Aaron confirms locations.

suppressMessages({ library(jsonlite); library(data.table); library(arrow) })

# ---- config -----------------------------------------------------------------
BUNDLE_DIR <- Sys.getenv("CONUS_SF_BUNDLE_DIR", "output/conus/sf_integration")  # TODO confirm
MODIF_DIR  <- Sys.getenv("CONUS_SF_MODIFIER_DIR", "fvs-modern/config/calibrated") # TODO confirm
OUT_JSON   <- "config/calibrated/conus_sf.json"
DRAW_DIR   <- "config/calibrated/draws"
N_DRAWS    <- as.integer(Sys.getenv("CONUS_SF_N_DRAWS", "200"))

# component -> bundle prefix, site term, modifier form, blend kappa (stress-test defaults)
COMPONENTS <- list(
  diameter_growth   = list(prefix="dg_sf",         site=list(source="BGI"),          modifier="common",         kappa=1500),
  height_growth     = list(prefix="hg_v8rd_sf_sf", site=list(source="BGI"),          modifier="common",         kappa=1500),
  height_diameter   = list(prefix="htdbh_v2split_sf", site=list(source="CSPI", v="v4"), modifier="none",        kappa=1500),
  height_crown_base = list(prefix="hcb_v2split_sf",   site=list(source="CSPI", v="v4"), modifier="trait_mediated", kappa=1500),
  survival          = list(prefix="surv_crz_sf",   site=list(source="CSPI", v="v4"), modifier="common",         kappa=1500),
  crown_recession   = list(prefix="cr_t2_sf",      site=list(source="CSPI", v="v4"), modifier="trait_mediated", kappa=400)
)

read_bundle <- function(prefix) {
  # TODO: read <prefix>_sf_fixed/_sf_gamma/_sf_species/_sf_re_*/_sf_manifest
  stop("TODO: implement bundle reader for ", prefix)
}

build_component <- function(name, spec) {
  b <- read_bundle(spec$prefix)
  list(
    form  = spec$prefix,
    shared = b$fixed,                                  # shared slopes
    species_independent = list(                        # Leg B
      trait_cols = b$trait_cols, gamma = b$gamma,
      scale_mean = b$scale_mean, scale_sd = b$scale_sd),
    species_dependent = list(                          # Leg A
      spcd = b$spcd, b0 = b$b0, n_obs = b$n_obs),
    blend = list(rule = "shrinkage", kappa = spec$kappa, min_n_legA = 500),
    site_term = spec$site,
    modifier  = if (spec$modifier == "none") list(form="none") else read_modifier(name, spec$modifier),
    posterior = list(n_draws = N_DRAWS, path = file.path("draws", paste0(name, "_draws.parquet")))
  )
}

read_modifier <- function(name, form) stop("TODO: read 62m modifier JSON for ", name)
write_draws   <- function(name, draws) stop("TODO: thin to N_DRAWS and write parquet for ", name)

main <- function() {
  dir.create(dirname(OUT_JSON), recursive = TRUE, showWarnings = FALSE)
  dir.create(DRAW_DIR,          recursive = TRUE, showWarnings = FALSE)
  comps <- lapply(names(COMPONENTS), function(n) build_component(n, COMPONENTS[[n]]))
  names(comps) <- names(COMPONENTS)
  cfg <- list(variant = "conus_sf", units = "metric", components = comps)
  write_json(cfg, OUT_JSON, auto_unbox = TRUE, pretty = TRUE, digits = 10)
  message("wrote ", OUT_JSON)
}
if (sys.nframe() == 0) main()
