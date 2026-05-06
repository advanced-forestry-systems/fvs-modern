// ============================================================================
// dg_kuehne2022.stan  (v3 — traits as fixed effects)
//
// FVS-CONUS diameter growth, Kuehne et al. 2022 form.
//
// v2: eta cap before exp() to prevent overflow.
// v3: Trait architecture refactored. W * gamma now enters the linear
//     predictor as a FIXED EFFECT (trait_effect), separate from the species
//     random intercept. z_sp is purely residual (zero-mean) with a tighter
//     sigma_sp prior.
//
// Linear predictor (annual diameter growth, cm/yr):
//   ln(dDBH_a) = b0 + trait_effect[sp] + z_sp[sp] + z_L1 + z_L2 + z_L3
//              + b1 ln(DBH) + b2 DBH + b3 ln((CR+0.2)/1.2)
//              + b4 ln(BAL_SW+0.01) + b5 BAL_HW + b6 ln(CSI)
//              + b7 BA x RD + b8 BAL x RD
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

  vector[N_obs] ln_dbh;
  vector[N_obs] dbh;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] ln_bal_sw_adj;
  vector[N_obs] bal_hw;
  vector[N_obs] ln_csi;
  vector[N_obs] ba_x_rd;
  vector[N_obs] bal_x_rd;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real b0;
  real b1;
  real b2;
  real b3;
  real b4;
  real b5;
  real b6;
  real b7;
  real b8;

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
  b0 ~ normal(-1.0, 2.0);
  b1 ~ normal( 0.3, 0.5);
  b2 ~ normal(-0.02, 0.05);
  b3 ~ normal( 0.5, 0.5);
  b4 ~ normal(-0.05, 0.1);
  b5 ~ normal(-0.005, 0.02);
  b6 ~ normal( 0.5, 0.5);
  b7 ~ normal( 0.0, 0.02);
  b8 ~ normal( 0.0, 0.02);

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

  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh
    + b2 * dbh
    + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj
    + b5 * bal_hw
    + b6 * ln_csi
    + b7 * ba_x_rd
    + b8 * bal_x_rd;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -30.0), 20.0);
  // lognormal variant: mirrors HG, avoids the variance-scaled normal collapse.
  // dg_obs_a can include some small negative values (measurement noise on
  // near-zero growth); skip them so lognormal sees only positive responses.
  for (i in 1:N_obs) {
    if (dg_obs_a[i] > 0.001) {
      target += lognormal_lpdf(dg_obs_a[i] | eta_safe[i], fmin(fmax(sigma / sqrt_years[i], 1e-4), 50.0));
    }
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] eta_gq =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw + b6 * ln_csi
    + b7 * ba_x_rd + b8 * bal_x_rd;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) {
    log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
  }
}
