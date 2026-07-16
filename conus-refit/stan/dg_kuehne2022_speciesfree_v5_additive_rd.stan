// ============================================================================
// dg_kuehne2022_speciesfree_v5_additive_rd.stan
//
// Architecture variant 5: refines v4_full per biometric guidance (May 15):
//   - Use rd_additive (sdi_additive1 / SDImax_brms) instead of Reineke RD
//   - Add SDI complexity ratio (sdi_additive1 / SDI_Reineke) as covariate
//     measuring uneven-aged / structural heterogeneity
//   - Replace BA*RD and BAL*RD precomputed terms with explicit continuous
//     interactions BA*rd_additive and BAL*rd_additive
//   - Add site x rd_additive interaction: site response strength depends on
//     stand density (productive sites tolerate higher density)
//
// All RD-based modifiers are continuous (not bound categorical) per Aaron's
// preference for continuous covariate multipliers over bounded indicators.
//
// Linear predictor:
//   eta = b0 + trait_effect[sp] + z_L1 + z_L2 + z_L3
//       + b1 ln(DBH) + b2 DBH + b3 ln_cr_adj
//       + b4 ln_bal_sw_adj + b5 BAL_HW
//       + b_site * ln_csi + b9 ln_csi^2
//       + b7 BA * rd_additive
//       + b8 BAL_SW * rd_additive
//       + b11 sdi_complexity
//       + b12 ln_csi * rd_additive    <-- continuous site x density modifier
//
// where:
//   b_site_i = b6 + z_L1_csi[L1_i] + W[sp_i,] * gamma_site
//   rd_additive = sdi_additive1 / SDImax_brms (continuous, ~0 to ~1.5)
//   sdi_complexity = sdi_additive1 / SDI_Reineke (continuous, ~0.4 to ~3.6)
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

  vector[N_obs] ba_metric;        // BA in m2/ha
  vector[N_obs] bal_sw_metric;    // BAL_SW in m2/ha
  vector[N_obs] rd_additive;      // sdi_additive1 / SDImax_brms
  vector[N_obs] sdi_complexity;   // sdi_additive1 / SDI_Reineke

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  vector[N_obs] ln_csi_sq = ln_csi .* ln_csi;
  vector[N_obs] ba_x_rdadd  = ba_metric     .* rd_additive;
  vector[N_obs] balsw_x_rdadd = bal_sw_metric .* rd_additive;
  vector[N_obs] lncsi_x_rdadd = ln_csi      .* rd_additive;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5;
  real b6;           // mean linear site coef
  real b7;           // BA x rd_additive
  real b8;           // BAL_SW x rd_additive
  real b9;           // quadratic site
  real b11;          // sdi_complexity
  real b12;          // ln_csi x rd_additive

  vector[P_trait] gamma;
  vector[P_trait] gamma_site;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_L1] z_L1_csi_raw;
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
  b0  ~ normal(-1.0, 2.0);
  b1  ~ normal( 0.3, 0.5);
  b2  ~ normal(-0.02, 0.05);
  b3  ~ normal( 0.5, 0.5);
  b4  ~ normal(-0.05, 0.1);
  b5  ~ normal(-0.005, 0.02);
  b6  ~ normal( 0.5, 0.5);
  b7  ~ normal( 0.0, 0.02);
  b8  ~ normal( 0.0, 0.02);
  b9  ~ normal(-0.1, 0.3);
  b11 ~ normal( 0.0, 0.3);       // sdi_complexity (no strong prior)
  b12 ~ normal( 0.0, 0.3);       // site x rd_additive interaction

  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.3);

  z_L1_raw     ~ std_normal();
  z_L2_raw     ~ std_normal();
  z_L3_raw     ~ std_normal();
  z_L1_csi_raw ~ std_normal();
  sigma_L1     ~ normal(0, 0.5);
  sigma_L2     ~ normal(0, 0.3);
  sigma_L3     ~ normal(0, 0.3);
  sigma_L1_csi ~ normal(0, 0.3);
  sigma        ~ normal(0, 0.5);

  vector[N_obs] b_site = b6 + z_L1_csi[L1_idx] + species_site_slope[sp_idx];

  vector[N_obs] eta =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site .* ln_csi + b9 * ln_csi_sq
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * lncsi_x_rdadd;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -30.0), 20.0);
  for (i in 1:N_obs) {
    if (dg_obs_a[i] > 0.001) {
      target += lognormal_lpdf(dg_obs_a[i] | eta_safe[i],
                fmin(fmax(sigma / sqrt_years[i], 1e-4), 50.0));
    }
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] b_site_gq = b6 + z_L1_csi[L1_idx] + species_site_slope[sp_idx];
  vector[N_obs] eta_gq =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site_gq .* ln_csi + b9 * ln_csi_sq
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * lncsi_x_rdadd;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
