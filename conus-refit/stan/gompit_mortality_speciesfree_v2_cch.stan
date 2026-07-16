// ============================================================================
// gompit_mortality_speciesfree_v2_cch.stan
//
// v2 of the species-free mortality model, motivated by Greg Johnson's CONUS
// mortality residual analysis (Johnson, Marshall, Weiskittel 2026-05-20).
//
// Two additions over v1 (gompit_mortality_speciesfree.stan):
//   1. Crown closure at tree tip (cch), entered as a quadratic (cch + cch^2).
//      Greg's base-rate residuals show strong U-shaped structure in cch that
//      his proposed model captures with a cch^beta term. We were missing this
//      variable entirely. A quadratic captures the U-shape robustly without an
//      estimated power exponent (which samples poorly in hierarchical HMC).
//   2. Quadratic crown ratio (cr + cr^2). Greg's residuals show a steep,
//      nonlinear rise in mortality at low crown ratio that the v1 linear cr
//      term cannot capture. The quadratic adds curvature.
//
// Everything else matches v1: cloglog (gompit) link with exposure offset,
// trait_effect = W * gamma, nested L1/L2/L3 ecoregion REs, forest type RE.
//
// Sign convention (unchanged from v1): higher eta = higher mortality.
//   P_mort_annual = 1 - exp(-exp(eta))
//   P_surv_T = exp(-exp(eta) * T)
//
// Linear predictor:
//   eta = b0 + trait_effect[sp]
//       + z_L1 + z_L2 + z_L3 + z_FT
//       + b1 DBH + b2 DBH^2
//       + b3 CR + b3b CR^2          <- NEW quadratic crown ratio
//       + b4 ln_csi
//       + b5 BAL + b6 sqrt_ba_rd
//       + b7 cch + b7b cch^2        <- NEW crown closure at tree tip
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  array[N_obs] int<lower=0, upper=1> alive;
  vector<lower=0>[N_obs] T_years;

  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;
  vector[N_obs] cr_init;
  vector[N_obs] ln_csi;
  vector[N_obs] bal_metric;
  vector[N_obs] sqrt_ba_rd;
  vector[N_obs] cch;              // NEW: crown closure at tree tip (fraction of acre)

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  // Standardize cch to mean 0 sd 1 for sampling stability; keep raw for the
  // quadratic so the curvature is on the standardized scale.
  real cch_mean = mean(cch);
  real cch_sd   = sd(cch) > 0 ? sd(cch) : 1.0;
  vector[N_obs] cch_z   = (cch - cch_mean) / cch_sd;
  vector[N_obs] cch_z_sq = cch_z .* cch_z;

  vector[N_obs] cr_sq = cr_init .* cr_init;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5; real b6;
  real b3b;                       // NEW: crown ratio quadratic
  real b7; real b7b;              // NEW: cch linear and quadratic

  vector[P_trait] gamma;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;

  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;
}
transformed parameters {
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
  vector[N_FT] z_FT = sigma_FT * z_FT_raw;
}
model {
  b0 ~ normal(-5.0, 2.0);
  b1 ~ normal( 0.0, 0.05);
  b2 ~ normal( 0.0, 0.001);
  b3 ~ normal(-1.0, 1.0);     // CR linear: vigorous trees die less
  b3b ~ normal( 0.0, 0.5);    // CR quadratic: curvature at low CR
  b4 ~ normal(-0.3, 0.3);
  b5 ~ normal( 0.01, 0.05);
  b6 ~ normal( 0.1, 0.3);
  b7 ~ normal( 0.0, 0.5);     // cch linear
  b7b ~ normal( 0.0, 0.5);    // cch quadratic (the U-shape Greg found)

  gamma ~ normal(0, 0.5);

  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);

  vector[N_obs] eta =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * cr_init + b3b * cr_sq
    + b4 * ln_csi
    + b5 * bal_metric + b6 * sqrt_ba_rd
    + b7 * cch_z + b7b * cch_z_sq;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -20.0), 5.0);

  for (i in 1:N_obs) {
    real log_p_surv = -exp(eta_safe[i]) * T_years[i];
    real log_p_mort = log1m_exp(log_p_surv);
    target += alive[i] == 1 ? log_p_surv : log_p_mort;
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] p_mort_annual;
  {
    vector[N_obs] eta =
        b0 + trait_effect[sp_idx]
      + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
      + b1 * dbh + b2 * dbh_sq
      + b3 * cr_init + b3b * cr_sq
      + b4 * ln_csi
      + b5 * bal_metric + b6 * sqrt_ba_rd
      + b7 * cch_z + b7b * cch_z_sq;
    for (i in 1:N_obs) {
      real eta_s = fmin(fmax(eta[i], -20.0), 5.0);
      real log_p_surv = -exp(eta_s) * T_years[i];
      real log_p_mort = log1m_exp(log_p_surv);
      log_lik[i] = alive[i] == 1 ? log_p_surv : log_p_mort;
      p_mort_annual[i] = 1.0 - exp(-exp(eta_s));
    }
  }
}
