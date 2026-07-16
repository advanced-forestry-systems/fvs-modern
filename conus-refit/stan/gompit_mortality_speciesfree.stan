// ============================================================================
// gompit_mortality_speciesfree.stan
//
// Species-free B1 architecture for tree mortality. Cloglog (gompit) with
// exposure offset for variable measurement intervals. Mirrors the
// DG_Kuehne v4_full pattern:
//   - trait_effect = W * gamma (no z_sp random intercept)
//   - Forest type random effect (z_FT)
//   - Nested EPA L1 / L2 / L3 ecoregion REs
//
// Link function: cloglog (complementary log log / gompit)
//   cloglog(P_mort_annual) = eta
//   P_mort_annual = 1 - exp(-exp(eta))
//
// With exposure offset for measurement interval T:
//   P_surv_T = exp(-exp(eta) * T)
//   P_mort_T = 1 - exp(-exp(eta) * T)
//   alive_T ~ Bernoulli(P_surv_T)
//
// Linear predictor:
//   eta = b0 + trait_effect[sp]
//       + z_L1 + z_L2 + z_L3 + z_FT
//       + b1 DBH + b2 DBH^2
//       + b3 CR + b4 ln_csi
//       + b5 BAL + b6 sqrt_ba_rd
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  array[N_obs] int<lower=0, upper=1> alive;   // survived (1) or died (0)
  vector<lower=0>[N_obs] T_years;              // measurement interval

  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;
  vector[N_obs] cr_init;
  vector[N_obs] ln_csi;
  vector[N_obs] bal_metric;
  vector[N_obs] sqrt_ba_rd;       // sqrt(BA * RD)

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5; real b6;

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
  b0 ~ normal(-5.0, 2.0);     // baseline annual mortality on cloglog scale
  b1 ~ normal( 0.0, 0.05);    // DBH
  b2 ~ normal( 0.0, 0.001);   // DBH^2
  b3 ~ normal(-1.0, 1.0);     // CR: vigorous trees die less
  b4 ~ normal(-0.3, 0.3);     // ln_csi: better sites lower mortality
  b5 ~ normal( 0.01, 0.05);   // BAL: more competition higher mortality
  b6 ~ normal( 0.1, 0.3);     // sqrt(BA * RD): density-dependent mortality

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
    + b3 * cr_init + b4 * ln_csi
    + b5 * bal_metric + b6 * sqrt_ba_rd;

  // Clamp eta for numerical stability
  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -20.0), 5.0);

  // log P_surv = -exp(eta) * T
  // alive ~ Bernoulli(exp(-exp(eta) * T))
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
      + b3 * cr_init + b4 * ln_csi
      + b5 * bal_metric + b6 * sqrt_ba_rd;
    for (i in 1:N_obs) {
      real eta_s = fmin(fmax(eta[i], -20.0), 5.0);
      real log_p_surv = -exp(eta_s) * T_years[i];
      real log_p_mort = log1m_exp(log_p_surv);
      log_lik[i] = alive[i] == 1 ? log_p_surv : log_p_mort;
      p_mort_annual[i] = 1.0 - exp(-exp(eta_s));
    }
  }
}
