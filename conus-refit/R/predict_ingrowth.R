#!/usr/bin/env Rscript
# predict_ingrowth.R
# FIX 2 recruitment predictor for the projector: replaces the empirical per-variant ingrowth lookup
# with the fitted models. Two parts:
#   COUNT  (how many recruits): negative binomial v4 (ingrowth_negbinom_v2.stan). Population-average
#          (EPA random effects z_L = 0). Linear predictor transcribed from the Stan model header:
#            ln(lambda_per_year) = b0 + W_dom . gamma
#               + b1 ln(BA_m2ha + 1) + b2 ln(BAL_m2ha + 1) + b3 rd_sdimax
#               + b4 ln(CSI) + b5 ln(HT40_m + 1) + b6 clim_pca1
#          Recruits over dt years = dt * lambda_per_year. Trait standardization is deterministic from
#          the full species_traits file (NOT the fit subsample), so the count predictor is exact and
#          self-contained from the summary CSV alone.
#   SHARES (which species): trait-driven multinomial (ingrowth_species_composition_v2.stan):
#            eta[s] = alpha_0[s] + (X_std . gamma_cov') . W_sp[s] ; shares = softmax(eta)
#          where X_std is the standardized (ln_ba, ln_bal, RD, ln_csi) vector. The composition X
#          standardization is subsample-dependent, so it must come from the fit meta's
#          scale_mean/scale_sd (saved by the patched 36 from the next refit). recover_shares() errors
#          clearly if the meta lacks them.
# recruits_by_species = round(count_total * shares).
suppressPackageStartupMessages({ library(data.table) })

## ---- COUNT: negative binomial recruit rate (self-contained) ----
load_count_model <- function(summary_csv, meta_rds, traits_rds) {
  s <- fread(summary_csv); getm <- function(v) as.numeric(s[variable == v, mean][1])
  b <- sapply(0:6, function(i) getm(paste0("b", i)))            # b0..b6
  gv <- s[grepl("^gamma\\[", variable), mean]                  # trait gammas
  meta <- readRDS(meta_rds); trait_cols <- meta$trait_cols     # exact training trait set (9 cols)
  tr <- as.data.table(readRDS(traits_rds))
  trait_cols <- intersect(trait_cols, names(tr))
  stopifnot(length(trait_cols) == length(gv))                  # gammas must match trait count
  tmean <- sapply(tr[, trait_cols, with = FALSE], function(x) mean(x, na.rm = TRUE))
  tsd   <- sapply(tr[, trait_cols, with = FALSE], function(x) sd(x, na.rm = TRUE))
  Wtab  <- tr[, c("SPCD", trait_cols), with = FALSE]
  list(b = b, gamma = as.numeric(gv), trait_cols = trait_cols,
       tmean = tmean, tsd = tsd, Wtab = Wtab)
}

# expected recruits per YEAR for one stand (population-average; z_L = 0)
predict_recruit_rate <- function(m, BA_m2ha, BAL_m2ha, rd_sdimax, CSI, HT40_m, clim_pca1, dom_spcd) {
  w <- m$Wtab[SPCD == dom_spcd, m$trait_cols, with = FALSE]
  if (!nrow(w)) wz <- rep(0, length(m$trait_cols)) else
    wz <- (as.numeric(w[1]) - m$tmean) / m$tsd
  wz[!is.finite(wz)] <- 0
  eta <- m$b[1] + sum(wz * m$gamma) +
    m$b[2]*log(BA_m2ha + 1) + m$b[3]*log(BAL_m2ha + 1) + m$b[4]*rd_sdimax +
    m$b[5]*log(max(CSI, 0.01)) + m$b[6]*log(HT40_m + 1) + m$b[7]*clim_pca1
  exp(eta)                                                     # recruits per year
}

## ---- SHARES: trait-driven multinomial composition ----
# meta must carry scale_mean/scale_sd (patched 36 saves them); gamma_cov from the complete summary.
load_composition_model <- function(meta_rds, gamma_cov_summary_csv, traits_rds) {
  m <- readRDS(meta_rds)
  if (is.null(m$scale_mean) || is.null(m$scale_sd))
    stop("composition meta lacks scale_mean/scale_sd; re-fit 36 (patched) or reconstruct the seed-42 subsample standardization.")
  s <- fread(gamma_cov_summary_csv)
  P_cov <- length(m$cov_cols); P_tr <- length(m$trait_cols)
  gc <- matrix(NA_real_, P_tr, P_cov)
  for (p in seq_len(P_tr)) for (c in seq_len(P_cov))
    gc[p, c] <- as.numeric(s[variable == sprintf("gamma_cov[%d,%d]", p, c), mean][1])
  a0 <- sapply(seq_along(m$sp_levels), function(k) as.numeric(s[variable == sprintf("alpha_0[%d]", k), mean][1]))
  tr <- as.data.table(readRDS(traits_rds))
  Wsp <- as.matrix(tr[match(m$sp_levels, SPCD), m$trait_cols, with = FALSE])
  for (j in seq_len(ncol(Wsp))) { na <- is.na(Wsp[,j]); if (any(na)) Wsp[na,j] <- median(Wsp[!na,j], na.rm=TRUE)
    Wsp[,j] <- (Wsp[,j] - mean(Wsp[,j])) / sd(Wsp[,j]) }
  list(alpha_0 = a0, gamma_cov = gc, Wsp = Wsp, sp_levels = m$sp_levels,
       cov_cols = m$cov_cols, scale_mean = m$scale_mean, scale_sd = m$scale_sd)
}

# species share vector for one stand; X order must match meta$cov_cols (default ln_ba, ln_bal, RD, ln_csi)
predict_shares <- function(cm, ln_ba, ln_bal, RD, ln_csi) {
  X <- c(ln_ba = ln_ba, ln_bal = ln_bal, RD = RD, ln_csi = ln_csi)[cm$cov_cols]
  Xs <- (X - cm$scale_mean[cm$cov_cols]) / cm$scale_sd[cm$cov_cols]
  Xs[!is.finite(Xs)] <- 0
  eta <- cm$alpha_0 + as.numeric(cm$Wsp %*% as.numeric(t(cm$gamma_cov) %*% Xs))
  ex <- exp(eta - max(eta)); setNames(ex / sum(ex), cm$sp_levels)
}

if (sys.nframe() == 0) cat("predict_ingrowth.R loaded: load_count_model/predict_recruit_rate (count, self-contained);",
                           "load_composition_model/predict_shares (needs patched-36 meta).\n")
