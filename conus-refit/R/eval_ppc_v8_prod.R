##=============================================================================
## eval_ppc_v8_prod.R
##
## Posterior predictive checks for DG_Kuehne v8 production fit.
## Generates:
##   1. Marginal histogram observed vs replicated
##   2. Test statistic comparisons (mean, sd, q05, q95)
##   3. Per-species mean PPC (do species fit individually?)
##   4. Per-ecoregion mean PPC
##   5. Residual vs predicted (should be flat)
##   6. Residual vs DBH (should be flat)
##   7. Residual vs BGI (should be flat) — key check for v8 hockey-stick
##
## Output: ppc_v8_prod.png with 6 panels + ppc_v8_prod_table.csv
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(ggplot2)
  library(bit64)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

FIT_FILE   <- get_arg("fit", "calibration/output/conus/dg_kue/v8/dg_kuehne_v8_100k_prod_fit.rds")
META_FILE  <- get_arg("meta", sub("_fit\\.rds$", "_meta.rds", FIT_FILE))
PAIRS_FILE <- get_arg("pairs", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")
OUT_DIR    <- get_arg("out", "calibration/output/conus/dg_kue/v8/ppc")
N_DRAWS    <- as.integer(get_arg("n_draws", "200"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== eval_ppc_v8_prod.R ==\n")
cat("  fit:", FIT_FILE, "\n  out:", OUT_DIR, "\n  n_draws:", N_DRAWS, "\n\n")

fit  <- readRDS(FIT_FILE)
meta <- readRDS(META_FILE)
pairs <- as.data.table(readRDS(PAIRS_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))

sp_levels <- meta$sp_levels
L1_levels <- meta$L1_levels
L2_levels <- meta$L2_levels
L3_levels <- meta$L3_levels
FT_levels <- meta$FT_levels

# Reproduce data prep
pairs[, dg_obs_a := (DBH2 - DBH1) / YEARS]
pairs[, sqrt_years := sqrt(YEARS)]
pairs[, ln_dbh := log(DBH1)]
pairs[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
pairs[, ln_bal_sw_adj := log(BAL_SW1 + 0.01)]
pairs[, rd_additive := sdi_additive1 / SDImax_brms]
pairs[, sdi_complexity := sdi_additive1 / pmax(SDI1, 1.0)]

filt <- is.finite(pairs$DBH1) & pairs$DBH1 >= 2.54 & is.finite(pairs$DBH2) &
        is.finite(pairs$CR1) & pairs$CR1 > 0 & pairs$CR1 <= 1.0 &
        is.finite(pairs$YEARS) & pairs$YEARS >= 1 & pairs$YEARS <= 20 &
        !is.na(pairs$EPA_L1_CODE) & pairs$EPA_L1_CODE != "" &
        pairs$TREESTATUS1 == 1 & pairs$TREESTATUS2 == 1 &
        is.finite(pairs$BAL_SW1) & pairs$BAL_SW1 >= 0 &
        is.finite(pairs$BAL_HW1) & pairs$BAL_HW1 >= 0 &
        is.finite(pairs$rd_additive) & pairs$rd_additive > 0 & pairs$rd_additive < 3.0 &
        is.finite(pairs$sdi_complexity) & pairs$sdi_complexity > 0 & pairs$sdi_complexity < 10 &
        is.finite(pairs$BA1) & pairs$BA1 >= 0 &
        is.finite(pairs$bgi) &
        !is.na(pairs$FORTYPCD_cond1) & pairs$FORTYPCD_cond1 > 0 &
        pairs$dg_obs_a > 0.01 & pairs$dg_obs_a < 5.0
d <- pairs[filt]
d <- d[SPCD %in% sp_levels]
d[, sp_idx := match(SPCD, sp_levels)]
d[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
d[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
d[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
d[, FT_idx := match(as.integer(FORTYPCD_cond1), FT_levels)]
d <- d[!is.na(sp_idx) & !is.na(L1_idx) & !is.na(L2_idx) &
        !is.na(L3_idx) & !is.na(FT_idx)]
cat("After filters:", nrow(d), "rows; matching base fit:", meta$n_obs, "\n")

# Subsample to match the original 100K subsample seed
set.seed(42)
idx <- sort(sample.int(nrow(d), min(100000, nrow(d))))
d <- d[idx]
cat("Subsampled to:", nrow(d), "rows\n")

# Trait matrix
trait_cols <- meta$trait_cols
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
softwood_by_sp <- traits_sub$softwood
softwood_by_sp[is.na(softwood_by_sp)] <- 0
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}
softwood_per_tree <- softwood_by_sp[d$sp_idx]
softwood_per_tree_c <- softwood_per_tree - mean(softwood_per_tree)

# Get N_DRAWS posterior draws of key parameters
draws <- fit$draws(
  variables = c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9a","b9b",
                "b11","b12","b13","b14","b15",
                "sigma","trait_effect","species_site_slope","z_sp",
                "z_L1","z_L2","z_L3","z_FT","z_L1_bgi",
                paste0("gamma_site[", seq_len(ncol(W)), "]")),
  format = "draws_matrix"
)

# Thin to N_DRAWS
all_draws <- nrow(draws)
keep_idx <- round(seq(1, all_draws, length.out = N_DRAWS))
draws <- draws[keep_idx, ]
cat("Using", nrow(draws), "draws for PPC\n")

# Knots from meta
knot1 <- meta$bgi_knots[1]
knot2 <- meta$bgi_knots[2]
bgi <- d$bgi
clim_b1 <- bgi
clim_b2 <- pmax(bgi - knot1, 0)
clim_b3 <- pmax(bgi - knot2, 0)

# Pre-compute fixed predictor matrix for each obs
n_obs <- nrow(d)
sp_i <- d$sp_idx
L1_i <- d$L1_idx
L2_i <- d$L2_idx
L3_i <- d$L3_idx
FT_i <- d$FT_idx

# Generate y_rep for each draw
y_obs <- d$dg_obs_a
cat("Generating", N_DRAWS, "replicated datasets...\n")

# Storage: marginal stats per draw
draw_stats <- data.table(
  draw = 1:N_DRAWS,
  mean_yrep = NA_real_,
  sd_yrep = NA_real_,
  q05_yrep = NA_real_,
  q50_yrep = NA_real_,
  q95_yrep = NA_real_
)

# For per-species: store residuals (obs - pred mean) per species
sp_resid_mean <- numeric(length(sp_levels))
sp_n <- integer(length(sp_levels))

# For residual vs predictor: compute mean prediction across draws
eta_sum <- numeric(n_obs)

for (di in 1:N_DRAWS) {
  drow <- draws[di, ]
  # Extract scalars
  b0 <- drow[, "b0"]; b1 <- drow[, "b1"]; b2 <- drow[, "b2"]
  b3 <- drow[, "b3"]; b4 <- drow[, "b4"]; b5 <- drow[, "b5"]
  b6 <- drow[, "b6"]; b7 <- drow[, "b7"]; b8 <- drow[, "b8"]
  b9a <- drow[, "b9a"]; b9b <- drow[, "b9b"]
  b11 <- drow[, "b11"]
  b12 <- drow[, "b12"]; b13 <- drow[, "b13"]
  b14 <- drow[, "b14"]; b15 <- drow[, "b15"]
  sigma <- drow[, "sigma"]
  gamma_site <- as.numeric(drow[, grep("^gamma_site\\[", colnames(drow))])
  # Extract vectors
  trait_effect <- as.numeric(drow[, grep("^trait_effect\\[", colnames(drow))])
  species_site_slope <- as.numeric(drow[, grep("^species_site_slope\\[", colnames(drow))])
  z_sp <- as.numeric(drow[, grep("^z_sp\\[", colnames(drow))])
  z_L1 <- as.numeric(drow[, grep("^z_L1\\[", colnames(drow))])
  z_L2 <- as.numeric(drow[, grep("^z_L2\\[", colnames(drow))])
  z_L3 <- as.numeric(drow[, grep("^z_L3\\[", colnames(drow))])
  z_FT <- as.numeric(drow[, grep("^z_FT\\[", colnames(drow))])
  z_L1_bgi <- as.numeric(drow[, grep("^z_L1_bgi\\[", colnames(drow))])

  b_site_i <- b6 + z_L1_bgi[L1_i] + species_site_slope[sp_i]
  eta <- b0 + trait_effect[sp_i] + z_sp[sp_i] +
         z_L1[L1_i] + z_L2[L2_i] + z_L3[L3_i] + z_FT[FT_i] +
         b1 * d$ln_dbh + b2 * d$DBH1 + b3 * d$ln_cr_adj +
         b4 * d$ln_bal_sw_adj + b5 * d$BAL_HW1 +
         b_site_i * clim_b1 + b9a * clim_b2 + b9b * clim_b3 +
         b7 * (d$BA1 * 0.2296 * d$rd_additive) +
         b8 * (d$BAL_SW1 * d$rd_additive) +
         b11 * d$sdi_complexity +
         b12 * (bgi * d$rd_additive) +
         b13 * (bgi * d$ln_dbh) +
         b14 * (bgi * softwood_per_tree_c) +
         b15 * (bgi * d$ln_cr_adj)
  eta_safe <- pmin(pmax(eta, -30), 20)
  sigma_i <- pmin(pmax(sigma / d$sqrt_years, 1e-4), 50)
  y_rep <- rlnorm(n_obs, meanlog = eta_safe, sdlog = sigma_i)

  draw_stats[di, mean_yrep := mean(y_rep, na.rm = TRUE)]
  draw_stats[di, sd_yrep := sd(y_rep, na.rm = TRUE)]
  draw_stats[di, q05_yrep := quantile(y_rep, 0.05, na.rm = TRUE)]
  draw_stats[di, q50_yrep := quantile(y_rep, 0.50, na.rm = TRUE)]
  draw_stats[di, q95_yrep := quantile(y_rep, 0.95, na.rm = TRUE)]

  eta_sum <- eta_sum + eta_safe
  if (di %% 25 == 0) cat("  draw", di, "/", N_DRAWS, "\n")
}
eta_mean <- eta_sum / N_DRAWS
pred_mean <- exp(eta_mean)  # lognormal mean approx

# Observed test stats
obs_stats <- data.table(
  stat = c("mean","sd","q05","q50","q95"),
  obs = c(mean(y_obs), sd(y_obs),
          quantile(y_obs, 0.05), quantile(y_obs, 0.50), quantile(y_obs, 0.95)),
  yrep_mean = c(mean(draw_stats$mean_yrep), mean(draw_stats$sd_yrep),
                mean(draw_stats$q05_yrep), mean(draw_stats$q50_yrep),
                mean(draw_stats$q95_yrep)),
  yrep_q05 = c(quantile(draw_stats$mean_yrep, 0.05),
               quantile(draw_stats$sd_yrep, 0.05),
               quantile(draw_stats$q05_yrep, 0.05),
               quantile(draw_stats$q50_yrep, 0.05),
               quantile(draw_stats$q95_yrep, 0.05)),
  yrep_q95 = c(quantile(draw_stats$mean_yrep, 0.95),
               quantile(draw_stats$sd_yrep, 0.95),
               quantile(draw_stats$q05_yrep, 0.95),
               quantile(draw_stats$q50_yrep, 0.95),
               quantile(draw_stats$q95_yrep, 0.95))
)
obs_stats[, ppc_ok := obs >= yrep_q05 & obs <= yrep_q95]
cat("\n=== Test statistic PPC ===\n")
print(obs_stats)
fwrite(obs_stats, file.path(OUT_DIR, "ppc_test_stats.csv"))

# Residual: obs - exp(eta_mean)
resid_marginal <- y_obs - pred_mean

# Save residual diagnostics summary
resid_summary <- data.table(
  group = c("overall"),
  n = length(resid_marginal),
  mean = mean(resid_marginal),
  sd = sd(resid_marginal),
  median = median(resid_marginal),
  q01 = quantile(resid_marginal, 0.01),
  q99 = quantile(resid_marginal, 0.99)
)

# Per-species residual
sp_resid <- data.table(
  sp_idx = d$sp_idx,
  resid = resid_marginal
)[, .(n = .N, mean = mean(resid), sd = sd(resid)), by = sp_idx]
sp_resid[, SPCD := sp_levels[sp_idx]]
setorder(sp_resid, -n)
fwrite(sp_resid, file.path(OUT_DIR, "ppc_residual_by_species.csv"))

cat("\n=== Top 10 species by N (residuals) ===\n")
print(head(sp_resid, 10))

# Plot: 4-panel PPC summary
df_plot <- data.table(
  observed = y_obs,
  predicted = pred_mean,
  bgi = bgi,
  dbh = d$DBH1,
  resid = resid_marginal,
  sp_idx = d$sp_idx
)

p1 <- ggplot(df_plot, aes(x = observed)) +
  geom_histogram(bins = 60, fill = "#3B7DC4", alpha = 0.5) +
  geom_histogram(aes(x = predicted), bins = 60, fill = "#D45050", alpha = 0.4) +
  scale_x_continuous(limits = c(0, 2)) +
  labs(title = "Marginal distribution: obs (blue) vs pred (red)",
       x = "dg_obs_a (cm/yr)", y = "count") +
  theme_minimal()

p2 <- ggplot(df_plot, aes(x = predicted, y = resid)) +
  geom_hex(bins = 50) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "gam", se = FALSE, color = "yellow", linewidth = 0.7) +
  labs(title = "Residual vs predicted", x = "predicted (cm/yr)", y = "residual") +
  theme_minimal() + theme(legend.position = "none")

p3 <- ggplot(df_plot, aes(x = bgi, y = resid)) +
  geom_hex(bins = 50) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "gam", se = FALSE, color = "yellow", linewidth = 0.7) +
  labs(title = "Residual vs BGI (the v8 piecewise variable)",
       x = "BGI", y = "residual") +
  theme_minimal() + theme(legend.position = "none")

p4 <- ggplot(df_plot, aes(x = dbh, y = resid)) +
  geom_hex(bins = 50) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "Residual vs DBH", x = "DBH (cm)", y = "residual") +
  theme_minimal() + theme(legend.position = "none")

# Use patchwork if available, else save separately
out_png <- file.path(OUT_DIR, "ppc_v8_prod.png")
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = "DG_Kuehne v8 production PPC")
  ggsave(out_png, p_combined, width = 14, height = 10, dpi = 130)
} else {
  ggsave(out_png, p1, width = 7, height = 5, dpi = 130)
  ggsave(sub("\\.png$", "_p2.png", out_png), p2, width = 7, height = 5)
  ggsave(sub("\\.png$", "_p3.png", out_png), p3, width = 7, height = 5)
  ggsave(sub("\\.png$", "_p4.png", out_png), p4, width = 7, height = 5)
}

cat("\nSaved PPC figure:", out_png, "\n")
cat("Done.\n")
