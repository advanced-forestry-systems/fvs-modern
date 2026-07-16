#!/usr/bin/env Rscript
# build_sf_coeffs.R  (extended 2026-06-22 for v8_forest_eco: b9..b12, gamma_site,
# species_site_slope, sigma_FT, sigma_L1_csi). Pulls posterior-mean coefficients
# and residual sigma from a species-free fit summary.csv (no multi-GB fit load).
#
# DG v8_forest_eco linear predictor (see stan/dg_kuehne2022_speciesfree_v8_forest_eco.stan):
#   b_site = b6 + z_L1_csi[L1] + species_site_slope[sp]
#   eta = b0 + trait_effect[sp] + z_L1 + z_L2 + z_L3 + z_FT
#       + b1 ln(DBH) + b2 DBH + b3 ln_cr_adj + b4 ln_bal_sw_adj + b5 bal_hw
#       + b_site*ln_sicond + b9 ln_sicond^2 + b7 ccfl1 + b8 is_plantation
#       + b10 ln_elev + b11 sdi_complexity + b12 (ln_sicond*rd_additive)
#   dg_obs_a ~ lognormal(eta, sigma/sqrt(YEARS))
#   PREDICTION CONVENTION (median exp(eta) vs lognormal-mean exp(eta+sigma_eff^2/2))
#   resolved empirically in the projector; sigma=1.765 makes it a ~4.7x lever at YEARS=1.
#
# CLI: --summary=PATH --meta=PATH --component=dg|hg|mort --out=PATH
suppressPackageStartupMessages(library(data.table))
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
args <- commandArgs(trailingOnly=TRUE)
ga <- function(n,d=NULL){m=grep(paste0("^--",n,"="),args,value=TRUE);if(!length(m))return(d);sub(paste0("^--",n,"="),"",m[1])}
SUMMARY <- ga("summary"); META <- ga("meta"); COMP <- ga("component","dg"); OUT <- ga("out")
stopifnot(file.exists(SUMMARY))
S <- fread(SUMMARY); setnames(S, names(S)[1], "var")
gv <- function(p){ x <- S[var==p, mean]; if(length(x)) x[1] else NA_real_ }
meta <- if(!is.null(META) && file.exists(META)) readRDS(META) else NULL
sp_levels <- meta$sp %||% meta$prep_meta$sp %||% meta$sp_levels %||% NA

# number of b coefficients present (b0..bK)
bmax <- max(c(0, as.integer(sub("^b","", grep("^b[0-9]+$", S$var, value=TRUE)))), na.rm=TRUE)
bundle <- list(
  component=COMP, summary_path=SUMMARY,
  b = setNames(sapply(0:bmax, function(i) gv(paste0("b",i))), paste0("b",0:bmax)),
  gamma            = S[grepl("^gamma\\[",var), mean],
  gamma_site       = S[grepl("^gamma_site\\[",var), mean],
  trait_effect     = S[grepl("^trait_effect\\[",var), mean],
  species_site_slope = S[grepl("^species_site_slope\\[",var), mean],
  sigma   = gv("sigma"),
  sigma_L1= gv("sigma_L1"), sigma_L2= gv("sigma_L2"), sigma_L3= gv("sigma_L3"),
  sigma_FT= gv("sigma_FT"), sigma_L1_csi= gv("sigma_L1_csi"),
  sp_levels = sp_levels,
  L1=meta$L1, L2=meta$L2, L3=meta$L3, FT=meta$FT, trait_cols=meta$trait_cols,
  cspi_shift = meta$cspi_shift %||% meta$prep_meta$cspi_shift %||% 1.0,
  note = "Ecoregion/FT REs (z_L1/L2/L3/FT, z_L1_csi) NOT in summary; population mean sets them 0. Prediction convention resolved in projector.")
if(!is.null(OUT)) saveRDS(bundle, OUT)
cat("component:", COMP, "| b0..b",bmax,": ", paste(round(bundle$b,4),collapse=", "),
    "\nsigma:", round(bundle$sigma,4),
    "| n gamma:", length(bundle$gamma), "| n gamma_site:", length(bundle$gamma_site),
    "| n trait_effect:", length(bundle$trait_effect),
    "| n species_site_slope:", length(bundle$species_site_slope),
    "| n species:", length(sp_levels), "\n", sep="")
if(!is.null(OUT)) cat("saved:", OUT, "\n")
