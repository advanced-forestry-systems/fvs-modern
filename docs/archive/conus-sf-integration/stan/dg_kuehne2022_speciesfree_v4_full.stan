// ============================================================================
// dg_kuehne2022_speciesfree_v4_full.stan
//
// Architecture variant 4 (FULL): the proposed "best" combination.
//   - Quadratic site term: b6 + b9 * ln_csi^2 (captures saturation)
//   - L1-varying site slope: z_L1_csi[L1] deviation from national mean
//   - Trait-modulated site slope: W * gamma_site (species responsiveness)
//   - Site x BAL interaction: b10 * (ln_csi * ln_bal_sw_adj)
//
// Effective site slope per observation:
//   b_site_i = b6 + z_L1_csi[L1_i] + species_site_slope[sp_i]
//                 + b10 * ln_bal_sw_adj_i
//
// Full linear predictor:
//   ln(dDBH_a) = b0 + trait_effect[sp] + z_L1 + z_L2 + z_L3
//              + b1 ln(DBH) + b2 DBH + b3 ln((CR+0.2)/1.2)
//              + b4 ln(BAL_SW+0.01) + b5 BAL_HW
//              + b_site_i * ln(CSI) + b9 ln(CSI)^2
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
transformed data {
  vector[N_obs] ln_csi_sq = ln_csi .* ln_csi;
}
parameters {
  real b0; real b1; real b2; real b3; real b4; real b5; real b6; real b7; real b8;
  real b9;                           // quadratic site
  real b10;                          // site x BAL interaction
  vector[P_trait] gamma;
  vector[P_trait] gamma_site;        // trait-modulated site slope
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_L1] z_L1_csi_raw;         // L1-varying site slope deviations
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_L1_csi;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N_sp] trait_effect;
  vector[N_sp] species_site_slope;
  if (P_trait > 0) {
    trait_effect       = W * gamma;
    species_site_slope = W * gamma_site;
  } else {
    trait_effect       = rep_vector(0.0, N_sp);
    species_site_slope = rep_vector(0.0, N_sp);
  }
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
  vector[N_L1] z_L1_csi = sigma_L1_csi * z_L1_csi_raw;
}
model {
  b0 ~ normal(-1.0, 2.0); b1 ~ normal(0.3, 0.5); b2 ~ normal(-0.02, 0.05);
  b3 ~ normal(0.5, 0.5);  b4 ~ normal(-0.05, 0.1); b5 ~ normal(-0.005, 0.02);
  b6 ~ normal(0.5, 0.5);             // national mean linear site coef (positive prior)
  b7 ~ normal(0.0, 0.02); b8 ~ normal(0.0, 0.02);
  b9 ~ normal(-0.1, 0.3);            // quadratic site
  b10 ~ normal(0.0, 0.3);            // site x BAL interaction
  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.3);
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_L1_csi_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_L1_csi ~ normal(0, 0.3);
  sigma ~ normal(0, 0.5);

  vector[N_obs] b_site =
      b6
    + z_L1_csi[L1_idx]
    + species_site_slope[sp_idx]
    + b10 * ln_bal_sw_adj;

  vector[N_obs] eta =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site .* ln_csi
    + b9 * ln_csi_sq
    + b7 * ba_x_rd + b8 * bal_x_rd;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -30.0), 20.0);
  for (i in 1:N_obs) {
    if (dg_obs_a[i] > 0.001) {
      target += lognormal_lpdf(dg_obs_a[i] | eta_safe[i], fmin(fmax(sigma / sqrt_years[i], 1e-4), 50.0));
    }
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] b_site_gq =
      b6
    + z_L1_csi[L1_idx]
    + species_site_slope[sp_idx]
    + b10 * ln_bal_sw_adj;
  vector[N_obs] eta_gq =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site_gq .* ln_csi
    + b9 * ln_csi_sq
    + b7 * ba_x_rd + b8 * bal_x_rd;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
