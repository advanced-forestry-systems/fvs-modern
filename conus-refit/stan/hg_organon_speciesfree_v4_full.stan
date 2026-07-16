// ============================================================================
// hg_organon_speciesfree_v4_full.stan
//
// Height growth (ORGANON form) species-free B1 with full architecture:
//   - Quadratic site term (a9 ln_csi^2)
//   - L1-varying site slope (z_L1_csi)
//   - Trait-modulated site slope (species_site_slope = W * gamma_site)
//   - Site x BAL interaction (a10 * ln_csi * bal_log)
//
// Mirrors the DG_Kuehne v4_full winning architecture identified by LOO ELPD
// comparison on 2026-05-15. The HG covariates differ from DG_Kue (single BAL,
// no BAL_SW/HW split), but the site-productivity architecture is identical.
//
// Linear predictor:
//   ln(HTG_a) = a0 + trait_effect[sp] + z_L1 + z_L2 + z_L3
//             + a1 ln(DBH) + a2 ln(HT) + a3 ln((CR+0.2)/1.2)
//             + b_site * ln(CSI) + a9 ln(CSI)^2
//             + a5 BAL + a6 BA + a7 SLOPE + a8 cos(ASP)
// where b_site_i = a4 + z_L1_csi[L1_i] + W[sp_i,] * gamma_site
//                       + a10 * bal_log_i
//
// Linear-site coefficient prior favors positive site response on growth.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=0> P_trait;

  vector[N_obs] hg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] ln_dbh;
  vector[N_obs] ln_ht;          // ln of total height at t1 (m)
  vector[N_obs] ln_cr_adj;
  vector[N_obs] bal_log;        // log(BAL + 5)
  vector[N_obs] ln_csi;         // log of climate site index
  vector[N_obs] ba_metric;      // basal area, m2/ha
  vector[N_obs] slope_pct;
  vector[N_obs] cos_aspect;

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
  real a0; real a1; real a2; real a3;
  real a4;                        // mean linear site slope
  real a5; real a6; real a7; real a8;
  real a9;                        // quadratic site
  real a10;                       // site x BAL interaction
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
  a0 ~ normal(-2.0, 2.0);
  a1 ~ normal(0.5, 0.5);          // ln(DBH) - positive juvenile growth
  a2 ~ normal(-0.3, 0.3);         // ln(HT) - taller trees grow less
  a3 ~ normal(0.5, 0.5);          // CR
  a4 ~ normal(0.3, 0.3);          // linear site slope (positive prior)
  a5 ~ normal(-0.1, 0.3);         // BAL
  a6 ~ normal(0.0, 0.3);          // BA
  a7 ~ normal(0.0, 0.05);         // slope
  a8 ~ normal(0.0, 0.05);         // cos asp
  a9 ~ normal(-0.05, 0.2);        // quadratic site
  a10 ~ normal(0.0, 0.2);         // site x BAL interaction

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
      a4
    + z_L1_csi[L1_idx]
    + species_site_slope[sp_idx]
    + a10 * bal_log;

  vector[N_obs] eta =
      a0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a1 * ln_dbh + a2 * ln_ht + a3 * ln_cr_adj
    + b_site .* ln_csi + a9 * ln_csi_sq
    + a5 * bal_log + a6 * ba_metric
    + a7 * slope_pct + a8 * cos_aspect;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -10.0), 5.0);
  for (i in 1:N_obs) {
    if (hg_obs_a[i] > 0.001) {
      target += lognormal_lpdf(hg_obs_a[i] | eta_safe[i], fmin(fmax(sigma / sqrt_years[i], 1e-4), 5.0));
    }
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] b_site_gq =
      a4
    + z_L1_csi[L1_idx]
    + species_site_slope[sp_idx]
    + a10 * bal_log;
  vector[N_obs] eta_gq =
      a0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a1 * ln_dbh + a2 * ln_ht + a3 * ln_cr_adj
    + b_site_gq .* ln_csi + a9 * ln_csi_sq
    + a5 * bal_log + a6 * ba_metric
    + a7 * slope_pct + a8 * cos_aspect;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 5.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(hg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
