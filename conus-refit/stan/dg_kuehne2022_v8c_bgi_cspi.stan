// ============================================================================
// dg_kuehne2022_v8c_bgi_cspi.stan
//
// Architecture variant 8c (2026-05-25): v8 BGI nonlinearity PLUS an additive
// CSPI v4 productivity term.
//
// Motivation: Aaron asked whether a robust plot-level site productivity
// metric (CSPI v4, the random-forest predicted site index with 44 covariates)
// can improve diameter growth predictions beyond what BGI plus the L1/L2/L3
// ecoregion REs already capture. The htdbh_unified static height-DBH model
// found a_cspi = +0.171 (SE 0.012) for productivity, but that has not been
// tested in the increment equations.
//
// v8c modifications relative to v8:
//
//   1. Add `vector[N_obs] ln_cspi_shift` to the data block, computed in the
//      R driver as log(max(cspi_v4, 0.1)). cspi_v4 typically ranges 11-40 m,
//      so ln_cspi_shift is roughly 2.3-3.7 with sd ~0.16.
//   2. Add `real a_cspi` parameter with informative prior N(0.15, 0.1)
//      centered on the htdbh result.
//   3. Add `+ a_cspi * ln_cspi_shift` to eta and eta_gq.
//   4. Everything else identical to v8 (3-piece BGI basis, all interactions,
//      species RE, trait gamma, gamma_site, etc.).
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

  vector[N_obs] bgi;
  vector[N_obs] ba_metric;
  vector[N_obs] bal_sw_metric;
  vector[N_obs] rd_additive;
  vector[N_obs] sdi_complexity;
  vector[N_obs] softwood;

  vector[N_obs] ln_cspi_shift;       // v8c: NEW productivity covariate

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;

  real bgi_knot1;
  real bgi_knot2;
}
transformed data {
  vector[N_obs] bgi_b1 = bgi;
  vector[N_obs] bgi_b2;
  vector[N_obs] bgi_b3;
  for (i in 1:N_obs) {
    bgi_b2[i] = fmax(bgi[i] - bgi_knot1, 0.0);
    bgi_b3[i] = fmax(bgi[i] - bgi_knot2, 0.0);
  }
  vector[N_obs] ba_x_rdadd     = ba_metric     .* rd_additive;
  vector[N_obs] balsw_x_rdadd  = bal_sw_metric .* rd_additive;
  vector[N_obs] bgi_x_rdadd    = bgi           .* rd_additive;
  vector[N_obs] bgi_x_lndbh    = bgi           .* ln_dbh;
  vector[N_obs] bgi_x_softwood = bgi           .* softwood;
  vector[N_obs] bgi_x_lncr     = bgi           .* ln_cr_adj;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5;
  real b6;
  real b7; real b8;
  real b9a; real b9b;
  real b11;
  real b12;
  real b13; real b14; real b15;
  real a_cspi;                       // v8c: NEW

  vector[P_trait] gamma;
  vector[P_trait] gamma_site;

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;
  vector[N_L1] z_L1_bgi_raw;

  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;
  real<lower=0> sigma_L1_bgi;

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
  vector[N_sp] z_sp     = sigma_sp     * z_sp_raw;
  vector[N_L1] z_L1     = sigma_L1     * z_L1_raw;
  vector[N_L2] z_L2     = sigma_L2     * z_L2_raw;
  vector[N_L3] z_L3     = sigma_L3     * z_L3_raw;
  vector[N_FT] z_FT     = sigma_FT     * z_FT_raw;
  vector[N_L1] z_L1_bgi = sigma_L1_bgi * z_L1_bgi_raw;
}
model {
  b0  ~ normal(-1.0, 2.0);
  b1  ~ normal( 0.3, 0.5);
  b2  ~ normal(-0.02, 0.05);
  b3  ~ normal( 0.5, 0.5);
  b4  ~ normal(-0.05, 0.1);
  b5  ~ normal(-0.005, 0.02);
  b6  ~ normal( 0.05, 0.10);
  b7  ~ normal( 0.0, 0.02);
  b8  ~ normal( 0.0, 0.02);
  b9a ~ normal( 0.0, 0.10);
  b9b ~ normal( 0.0, 0.10);
  b11 ~ normal( 0.0, 0.3);
  b12 ~ normal( 0.0, 0.1);
  b13 ~ normal( 0.0, 0.05);
  b14 ~ normal( 0.0, 0.10);
  b15 ~ normal( 0.0, 0.10);
  a_cspi ~ normal(0.15, 0.10);     // v8c: informed by htdbh result (0.171)

  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.05);

  z_sp_raw     ~ std_normal();
  z_L1_raw     ~ std_normal();
  z_L2_raw     ~ std_normal();
  z_L3_raw     ~ std_normal();
  z_FT_raw     ~ std_normal();
  z_L1_bgi_raw ~ std_normal();

  sigma_sp     ~ normal(0, 0.15);
  sigma_L1     ~ normal(0, 0.5);
  sigma_L2     ~ normal(0, 0.3);
  sigma_L3     ~ normal(0, 0.3);
  sigma_FT     ~ normal(0, 0.3);
  sigma_L1_bgi ~ normal(0, 0.05);

  sigma ~ normal(0, 0.5);

  vector[N_obs] b_site = b6 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];

  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site .* bgi_b1
    + b9a * bgi_b2 + b9b * bgi_b3
    + a_cspi * ln_cspi_shift             // v8c: NEW
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * bgi_x_rdadd
    + b13 * bgi_x_lndbh
    + b14 * bgi_x_softwood
    + b15 * bgi_x_lncr;

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
  vector[N_obs] b_site_gq = b6 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];
  vector[N_obs] eta_gq =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site_gq .* bgi_b1
    + b9a * bgi_b2 + b9b * bgi_b3
    + a_cspi * ln_cspi_shift             // v8c
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * bgi_x_rdadd
    + b13 * bgi_x_lndbh
    + b14 * bgi_x_softwood
    + b15 * bgi_x_lncr;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
