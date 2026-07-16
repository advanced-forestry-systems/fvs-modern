// ============================================================================
// survival_unified_v2_crz.stan
//
// Convergence fix for survival_unified.stan. The original left crown ratio (CR)
// on its raw [0,1] scale, so CR and CR^2 were ~0.97 collinear. That produced a
// b3/b3b ridge the sampler could not traverse (100k production fit: b3 = +8.0,
// b3b = -5.3, both many prior-SDs off; sigma_sp rhat 1.58, ESS 7). This variant
// standardizes CR exactly as the model already standardizes CCH, breaking the
// ridge. The model is otherwise identical.
//
// Link (cloglog on annual survival, exposure offset T):
//   annual_hazard      = exp(-eta)
//   P(survive T years) = exp(-exp(-eta) * T)
//   alive ~ Bernoulli(P_surv_T)
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
  vector[N_obs] cch;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  real cch_mean = mean(cch);
  real cch_sd   = sd(cch) > 0 ? sd(cch) : 1.0;
  vector[N_obs] cch_z    = (cch - cch_mean) / cch_sd;
  vector[N_obs] cch_z_sq = cch_z .* cch_z;

  // CR standardized (NEW): breaks the CR / CR^2 collinearity ridge.
  real cr_mean = mean(cr_init);
  real cr_sd   = sd(cr_init) > 0 ? sd(cr_init) : 1.0;
  vector[N_obs] cr_z    = (cr_init - cr_mean) / cr_sd;
  vector[N_obs] cr_z_sq = cr_z .* cr_z;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5; real b6;
  real b3b;
  real b7; real b7b;

  vector[P_trait] gamma;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;
  vector[N_sp] z_sp_raw;

  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;
  real<lower=0> sigma_sp;
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
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
}
model {
  // Priors on the SURVIVAL scale. b3/b3b are now on standardized CR.
  b0 ~ normal( 5.0, 2.0);
  b1 ~ normal( 0.0, 0.05);
  b2 ~ normal( 0.0, 0.001);
  b3 ~ normal( 0.5, 1.0);     // CR (standardized): vigorous trees survive more
  b3b ~ normal( 0.0, 0.5);    // CR curvature (standardized)
  b4 ~ normal( 0.3, 0.3);
  b5 ~ normal(-0.01, 0.05);
  b6 ~ normal(-0.1, 0.3);
  b7 ~ normal( 0.0, 0.5);
  b7b ~ normal( 0.0, 0.5);

  gamma ~ normal(0, 0.5);

  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  z_sp_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);
  sigma_sp ~ normal(0, 0.5);

  vector[N_obs] eta =
      b0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * cr_z + b3b * cr_z_sq
    + b4 * ln_csi
    + b5 * bal_metric + b6 * sqrt_ba_rd
    + b7 * cch_z + b7b * cch_z_sq;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -5.0), 20.0);

  for (i in 1:N_obs) {
    real log_p_surv = -exp(-eta_safe[i]) * T_years[i];
    real log_p_mort = log1m_exp(log_p_surv);
    target += alive[i] == 1 ? log_p_surv : log_p_mort;
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] p_surv_annual;
  {
    vector[N_obs] eta =
        b0 + trait_effect[sp_idx] + z_sp[sp_idx]
      + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
      + b1 * dbh + b2 * dbh_sq
      + b3 * cr_z + b3b * cr_z_sq
      + b4 * ln_csi
      + b5 * bal_metric + b6 * sqrt_ba_rd
      + b7 * cch_z + b7b * cch_z_sq;
    for (i in 1:N_obs) {
      real eta_s = fmin(fmax(eta[i], -5.0), 20.0);
      real log_p_surv = -exp(-eta_s) * T_years[i];
      real log_p_mort = log1m_exp(log_p_surv);
      log_lik[i] = alive[i] == 1 ? log_p_surv : log_p_mort;
      p_surv_annual[i] = exp(-exp(-eta_s));
    }
  }
}
