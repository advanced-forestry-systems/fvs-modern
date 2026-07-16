// ============================================================================
// hg_unified_v8_rd.stan
//
// Height growth, v8_rd architecture = v7 + relative-density-scaled competition.
//
// SAME as v7 (species RE + traits_v3 + 3-piece BGI + ln(HT)^2 asymptote) plus
// two interactions that scale competition by relative density (RD = additive
// SDI / SDImax), mirroring how RD enters DG v8 (BA*RD, BAL*RD) and HT-DBH
// (a_bard, a_blrd). RD is conspicuously absent from HG v7; this tests whether
// density-scaled competition improves height growth. RD is stand-derivable, so
// it is projection-safe (unlike BGI).
//
//   v7:      a5 BAL + a6 BA
//   v8_rd:   a5 BAL + a6 BA + a_bard (BA*RD) + a_blrd (BAL_raw*RD)
//
// Driver supplies rd_additive + bal_raw for any stan whose filename matches v8rd.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector[N_obs] hg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] ln_dbh;
  vector[N_obs] ln_ht;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] bal_log;
  vector[N_obs] bgi;
  vector[N_obs] ba_metric;
  vector[N_obs] slope_pct;
  vector[N_obs] cos_aspect;

  // RD interaction inputs (driver supplies when filename matches v8rd)
  vector[N_obs] rd_additive;
  vector[N_obs] bal_raw;

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
  vector[N_obs] ln_ht_sq;
  for (i in 1:N_obs) {
    bgi_b2[i] = fmax(bgi[i] - bgi_knot1, 0.0);
    bgi_b3[i] = fmax(bgi[i] - bgi_knot2, 0.0);
    ln_ht_sq[i] = ln_ht[i] * ln_ht[i];
  }
  vector[N_obs] bgi_x_bal = bgi .* bal_log;
  vector[N_obs] ba_x_rd   = ba_metric .* rd_additive;
  vector[N_obs] bal_x_rd  = bal_raw   .* rd_additive;
}
parameters {
  real a0; real a1; real a2;
  real a2_quad;
  real a3;
  real a4;
  real a5; real a6; real a7; real a8;
  real a9a; real a9b;
  real a10;
  real a_bard;                       // BA * RD
  real a_blrd;                       // BAL_raw * RD
  vector[P_trait] gamma;
  vector[P_trait] gamma_site;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;
  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_bgi_raw;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;
  real<lower=0> sigma_sp;
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
  vector[N_L1] z_L1     = sigma_L1     * z_L1_raw;
  vector[N_L2] z_L2     = sigma_L2     * z_L2_raw;
  vector[N_L3] z_L3     = sigma_L3     * z_L3_raw;
  vector[N_FT] z_FT     = sigma_FT     * z_FT_raw;
  vector[N_sp] z_sp     = sigma_sp     * z_sp_raw;
  vector[N_L1] z_L1_bgi = sigma_L1_bgi * z_L1_bgi_raw;
}
model {
  a0      ~ normal(-1.0, 2.0);
  a1      ~ normal( 0.0, 0.5);
  a2      ~ normal(-0.3, 0.5);
  a2_quad ~ normal(-0.05, 0.05);
  a3      ~ normal( 0.5, 0.5);
  a4      ~ normal( 0.05, 0.1);
  a5      ~ normal(-0.005, 0.02);
  a6      ~ normal(-0.005, 0.02);
  a7      ~ normal( 0.0, 0.01);
  a8      ~ normal( 0.0, 0.1);
  a9a     ~ normal( 0.0, 0.1);
  a9b     ~ normal( 0.0, 0.1);
  a10     ~ normal( 0.0, 0.05);
  a_bard  ~ normal( 0.0, 0.02);
  a_blrd  ~ normal( 0.0, 0.02);

  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.05);

  z_L1_raw     ~ std_normal();
  z_L2_raw     ~ std_normal();
  z_L3_raw     ~ std_normal();
  z_FT_raw     ~ std_normal();
  z_sp_raw     ~ std_normal();
  z_L1_bgi_raw ~ std_normal();

  sigma_L1     ~ normal(0, 0.5);
  sigma_L2     ~ normal(0, 0.3);
  sigma_L3     ~ normal(0, 0.3);
  sigma_FT     ~ normal(0, 0.3);
  sigma_sp     ~ normal(0, 0.5);
  sigma_L1_bgi ~ normal(0, 0.05);

  sigma ~ normal(0, 0.5);

  vector[N_obs] b_site = a4 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];

  vector[N_obs] eta =
      a0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + a1 * ln_dbh
    + a2 * ln_ht + a2_quad * ln_ht_sq
    + a3 * ln_cr_adj
    + b_site .* bgi_b1
    + a9a * bgi_b2 + a9b * bgi_b3
    + a5 * bal_log + a6 * ba_metric
    + a7 * slope_pct + a8 * cos_aspect
    + a10 * bgi_x_bal
    + a_bard * ba_x_rd + a_blrd * bal_x_rd;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -30.0), 20.0);
  for (i in 1:N_obs) {
    if (hg_obs_a[i] > 0.001) {
      target += lognormal_lpdf(hg_obs_a[i] | eta_safe[i],
                fmin(fmax(sigma / sqrt_years[i], 1e-4), 50.0));
    }
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] b_site_gq = a4 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];
  vector[N_obs] eta_gq =
      a0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + a1 * ln_dbh
    + a2 * ln_ht + a2_quad * ln_ht_sq
    + a3 * ln_cr_adj
    + b_site_gq .* bgi_b1
    + a9a * bgi_b2 + a9b * bgi_b3
    + a5 * bal_log + a6 * ba_metric
    + a7 * slope_pct + a8 * cos_aspect
    + a10 * bgi_x_bal
    + a_bard * ba_x_rd + a_blrd * bal_x_rd;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(hg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
