// ============================================================================
// dg_organon.stan  (v2 — traits as fixed effects)
//
// FVS-CONUS diameter growth, ORGANON SWO style with estimated K1, K2, K4.
//
// v2: Trait architecture refactored. W * gamma now enters the linear
//     predictor as a FIXED EFFECT (trait_effect), separate from the species
//     random intercept. z_sp is purely residual (zero-mean) with a tighter
//     sigma_sp prior. Gamma coefficients are directly interpretable as the
//     effect of each trait on species-level growth potential.
//
// Linear predictor (annual diameter growth, cm/yr):
//   ln(dDBH_a) = a0 + trait_effect[sp] + z_sp[sp] + z_L1 + z_L2 + z_L3
//              + a1 ln(DBH + K1) + a2 DBH^K2 + a3 ln((CR+0.2)/1.2)
//              + a4 ln(CSPI_shift) + a5 BAL/ln(DBH+K4)
//              + a6 sqrt(BA) + a7 BA x RD + a8 BAL x RD
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;

  vector[N_obs] dg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] dbh;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] ln_cspi_shift;
  vector[N_obs] bal;
  vector[N_obs] sqrt_ba;
  vector[N_obs] ba_x_rd;
  vector[N_obs] bal_x_rd;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real a0;
  real a1;
  real a2;
  real a3;
  real a4;
  real a5;
  real a6;
  real a7;
  real a8;

  real<lower=0.1>  K1;
  real<lower=0.05> K2;
  real<lower=0.5>  K4;

  vector[P_trait] gamma;
  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

  real<lower=0> sigma;
}
transformed parameters {
  // Traits as fixed effects (separate from species RE)
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
  // Species RE: purely residual, zero-mean
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
}
model {
  a0 ~ normal(-1.5, 2.0);
  a1 ~ normal( 0.4, 0.5);
  a2 ~ normal(-0.02, 0.05);
  a3 ~ normal( 0.8, 0.5);
  a4 ~ normal( 0.3, 0.5);
  a5 ~ normal(-0.005, 0.02);
  a6 ~ normal(-0.03, 0.05);
  a7 ~ normal( 0.0, 0.02);
  a8 ~ normal( 0.0, 0.02);

  K1 ~ normal(1.0, 1.0);
  K2 ~ normal(0.8, 0.3);
  K4 ~ normal(2.7, 1.0);

  gamma ~ normal(0, 0.5);
  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.15);        // tight: traits explain most species variation
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);

  sigma ~ normal(0, 0.5);

  vector[N_obs] ln_dbh_k1 = log(dbh + K1);
  vector[N_obs] dbh_k2;
  for (i in 1:N_obs) dbh_k2[i] = pow(dbh[i], K2);
  vector[N_obs] comp = bal ./ log(dbh + K4);

  vector[N_obs] eta =
      a0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a1 * ln_dbh_k1
    + a2 * dbh_k2
    + a3 * ln_cr_adj
    + a4 * ln_cspi_shift
    + a5 * comp
    + a6 * sqrt_ba
    + a7 * ba_x_rd
    + a8 * bal_x_rd;

  vector[N_obs] mu_safe;
  for (i in 1:N_obs) mu_safe[i] = exp(fmin(eta[i], 20.0));

  dg_obs_a ~ normal(mu_safe, sigma ./ sqrt_years);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_a;
  {
    vector[N_obs] ln_dbh_k1 = log(dbh + K1);
    vector[N_obs] dbh_k2;
    for (i in 1:N_obs) dbh_k2[i] = pow(dbh[i], K2);
    vector[N_obs] comp = bal ./ log(dbh + K4);
    {
      vector[N_obs] eta_gq =
          a0
        + trait_effect[sp_idx] + z_sp[sp_idx]
        + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
        + a1 * ln_dbh_k1 + a2 * dbh_k2 + a3 * ln_cr_adj
        + a4 * ln_cspi_shift + a5 * comp + a6 * sqrt_ba
        + a7 * ba_x_rd + a8 * bal_x_rd;
      for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
    }
  }
  for (i in 1:N_obs) {
    log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
  }
}
