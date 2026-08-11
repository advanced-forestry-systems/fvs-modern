## Localized maximum SDI for CONUS growth-and-yield models (R companion)
##
## Replaces the species-weighted maximum SDI (biased ~+28%, near-zero plot-level skill)
## with a per-stand, data-derived maximum from an FIA-based surface. Returns a per-stand
## maximum SDI (trees/ha) any growth-and-yield engine can use as its density limit, and
## emits the FVS SDIMAX keyword block. Model-agnostic.
##
## Source fidelity order: (1) TreeMap 2022 30 m SDImax raster (Zenodo 10.5281/zenodo.19509367),
## (2) brms FIA plot table, (3) forest-type + geography fallback.

TPA_PER_HA <- 2.4710538  # trees/ha / this = trees/acre (FVS internal units)

## brms_csv cols: STATECD,UNITCD,COUNTYCD,PLOT,ID(plot_key),SDImax.mean,SDImax.median
load_brms_sdimax <- function(brms_csv) {
  b <- utils::read.csv(brms_csv, stringsAsFactors = FALSE)
  names(b) <- trimws(gsub('"', "", names(b)))
  stats::setNames(as.numeric(b$SDImax.mean), as.character(b$ID))
}

sdimax_for_plot <- function(plot_key, brms_lookup) {
  v <- brms_lookup[[as.character(plot_key)]]
  if (is.null(v) || !is.finite(v)) NA_real_ else v
}

## TreeMap SDImax raster lookup (needs terra); raster_path -> SDImax band
sdimax_for_coord <- function(lon, lat, raster_path) {
  if (!requireNamespace("terra", quietly = TRUE) || !file.exists(raster_path)) return(NA_real_)
  r  <- terra::rast(raster_path)
  pt <- terra::vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
  pt <- terra::project(pt, terra::crs(r))
  v  <- terra::extract(r, pt)[, 2]
  ifelse(is.finite(v) & v > 0, v, NA_real_)
}

## Per-stand SDIMAX keyword block: every species set to the localized stand value, so the
## BA-weighted stand maximum equals the localized value regardless of composition.
fvs_sdimax_keywords <- function(max_sdi_tph, maxsp, to_acre = TRUE) {
  val <- if (to_acre) max_sdi_tph / TPA_PER_HA else max_sdi_tph
  paste(sprintf("SDIMAX  %10d%10.1f", seq_len(maxsp), val), collapse = "\n")
}

## One-call resolver: raster preferred, then plot table; NA if nothing resolves.
localized_sdimax_keywords <- function(plot_key = NULL, lon = NULL, lat = NULL, maxsp = 120,
                                      brms_lookup = NULL, raster_path = NULL) {
  v <- NA_real_
  if (!is.null(lon) && !is.null(lat) && !is.null(raster_path)) v <- sdimax_for_coord(lon, lat, raster_path)
  if (!is.finite(v) && !is.null(plot_key) && !is.null(brms_lookup)) v <- sdimax_for_plot(plot_key, brms_lookup)
  if (!is.finite(v)) return(NA_character_)
  fvs_sdimax_keywords(v, maxsp)
}

## Non-circular validation: deviance explained of OBSERVED self-thinning by relative density
## under each max-SDI candidate. The better maximum is the one whose RD predicts observed
## density loss best. Needs mgcv.
selfthin_skill <- function(df, sdi_col = "SDI1",
                           maxsd_cols = c(brms = "SDImax_brms", fvs = "fvs_sdimax"),
                           tph1 = "TPH1", tph2 = "TPH2", years = "YEARS") {
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  df$mort <- -log(df[[tph2]] / df[[tph1]]) / df[[years]]
  df <- df[is.finite(df$mort), ]
  sapply(maxsd_cols, function(mc) {
    rd <- df[[sdi_col]] / df[[mc]]; ok <- is.finite(rd) & is.finite(df$mort)
    summary(mgcv::bam(df$mort[ok] ~ s(rd[ok], k = 8)))$dev.expl
  })
}
