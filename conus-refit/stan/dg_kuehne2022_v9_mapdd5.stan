// ============================================================================
// dg_kuehne2022_v9_mapdd5.stan
//
// Architecture variant 9 (2026-05-15): replace BGI with mapdd5 as climate
// site variable. Driven by climate exploration on 110K joined rows
// (job 9601250, 1:28 wall):
//   - mapdd5 univariate R^2 = 0.046 (vs bgi = 0.032)
//   - mapdd5 conditional R^2 added = 0.052 (vs bgi ~ 0.04)
//   - mapdd5 = mean annual precipitation x growing degree days >5C
//     Captures joint moisture + warmth, both positive-sign, biologically
//     interpretable.
//
// Architecture identical to v8 (rich nonlinearity + interactions) but
// with mapdd5 in place of bgi everywhere.
//
// Linear predictor:
//   eta = b0 + trait_effect[sp] + z_sp[sp]
//       + z_L1 + z_L2 + z_L3 + z_FT
//       + b1 ln(DBH) + b2 DBH + b3 ln_cr_adj
//       + b4 ln_bal_sw_adj + b5 bal_hw
//       + b_site * mapdd5_b1 + b9a mapdd5_b2 + b9b mapdd5_b3   <-- 3-piece basis
//       + b7  BA  * rd_additive
//       + b8  BAL_SW * rd_additive
//       + b11 sdi_complexity
//       + b12 mapdd5 * rd_additive
//       + b13 mapdd5 * ln_dbh
//       + b14 mapdd5 * softwood
//       + b15 mapdd5 * ln_cr_adj
//
// where:
//   b_site_i = b6 + z_L1_site[L1_i] + W[sp_i,] * gamma_site
//   mapdd5 has been z-scored in the R driver, so coefficients are on
//   the standardized scale.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector[N_obs] dg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] ln_dbh;
  vector[N_obs] dbh;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] ln_bal_sw_adj;
  vector[N_obs] bal_hw;

  vector[N_obs] mapdd5;            // z-scored mapdd5
  vector[N_obs] ba_metric;
  vector[N_obs] bal_sw_metric;
  vector[N_obs] rd_additive;
  vector[N_obs] sdi_complexity;
  vector[N_obs] softwood;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;

  real mapdd5_knot1;
  real mapdd5_knot2;
}
transformed data {
  vector[N_obs] mapdd5_b1 = mapdd5;
  vector[N_obs] mapdd5_b2;
  vector[N_obs] mapdd5_b3;
  for (i in 1:N_obs) {
    mapdd5_b2[i] = fmax(mapdd5[i] - mapdd5_knot1, 0.0);
    mapdd5_b3[i] = fmax(mapdd5[i] - mapdd5_knot2, 0.0);
  }

  vector[N_obs] ba_x_rdadd     = ba_metric     .* rd_additive;
  vector[N_obs] balsw_x_rdadd  = bal_sw_metric .* rd_additive;
  vector[N_obs] mapdd5_x_rdadd = mapdd5        .* rd_additive;
  vector[N_obs] mapdd5_x_lndbh = mapdd5        .* ln_dbh;
  vector[N_obs] mapdd5_x_softwood = mapdd5     .* softwood;
  vector[N_obs] mapdd5_x_lncr  = mapdd5        .* ln_cr_adj;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5;
  real b6;
  real b7; real b8;
  real b9a; real b9b;
  real b11;
  real b12; real b13; real b14; real b15;

  vector[P_trait] gamma;
  vector[P_trait] gamma_site;

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;
  vector[N_L1] z_L1_site_raw;

  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;
  real<lower=0> sigma_L1_site;

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
  vector[N_sp] z_sp      = sigma_sp      * z_sp_raw;
  vector[N_L1] z_L1      = sigma_L1      * z_L1_raw;
  vector[N_L2] z_L2      = sigma_L2      * z_L2_raw;
  vector[N_L3] z_L3      = sigma_L3      * z_L3_raw;
  vector[N_FT] z_FT      = sigma_FT      * z_FT_raw;
  vector[N_L1] z_L1_site = sigma_L1_site * z_L1_site_raw;
}
model {
  b0  ~ normal(-1.0, 2.0);
  b1  ~ normal( 0.3, 0.5);
  b2  ~ normal(-0.02, 0.05);
  b3  ~ normal( 0.5, 0.5);
  b4  ~ normal(-0.05, 0.1);
  b5  ~ normal(-0.005, 0.02);
  b6  ~ normal( 0.1, 0.2);     // standardized mapdd5; weak positive prior
  b7  ~ normal( 0.0, 0.02);
  b8  ~ normal( 0.0, 0.02);
  b9a ~ normal( 0.0, 0.2);
  b9b ~ normal( 0.0, 0.2);
  b11 ~ normal( 0.0, 0.3);
  b12 ~ normal( 0.0, 0.2);
  b13 ~ normal( 0.0, 0.1);
  b14 ~ normal( 0.0, 0.2);
  b15 ~ normal( 0.0, 0.2);

  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.1);

  z_sp_raw      ~ std_normal();
  z_L1_raw      ~ std_normal();
  z_L2_raw      ~ std_normal();
  z_L3_raw      ~ std_normal();
  z_FT_raw      ~ std_normal();
  z_L1_site_raw ~ std_normal();

  sigma_sp      ~ normal(0, 0.15);
  sigma_L1      ~ normal(0, 0.5);
  sigma_L2      ~ normal(0, 0.3);
  sigma_L3      ~ normal(0, 0.3);
  sigma_FT      ~ normal(0, 0.3);
  sigma_L1_site ~ normal(0, 0.1);

  sigma ~ normal(0, 0.5);

  vector[N_obs] b_site = b6 + z_L1_site[L1_idx] + species_site_slope[sp_idx];

  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site .* mapdd5_b1
    + b9a * mapdd5_b2 + b9b * mapdd5_b3
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * mapdd5_x_rdadd
    + b13 * mapdd5_x_lndbh
    + b14 * mapdd5_x_softwood
    + b15 * mapdd5_x_lncr;

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
  vector[N_obs] b_site_gq = b6 + z_L1_site[L1_idx] + species_site_slope[sp_idx];
  vector[N_obs] eta_gq =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site_gq .* mapdd5_b1
    + b9a * mapdd5_b2 + b9b * mapdd5_b3
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * mapdd5_x_rdadd
    + b13 * mapdd5_x_lndbh
    + b14 * mapdd5_x_softwood
    + b15 * mapdd5_x_lncr;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
