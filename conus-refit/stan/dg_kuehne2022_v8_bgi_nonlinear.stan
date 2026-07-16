// ============================================================================
// dg_kuehne2022_v8_bgi_nonlinear.stan
//
// Architecture variant 8 (2026-05-15): rich BGI nonlinearity + interactions.
//
// Motivation: in v7 the linear BGI coefficient remained negative even with
// quadratic term, ecoregion-varying slope (z_L1_bgi), species-trait
// modulation (gamma_site), and BGI*RD interaction. Aaron's directive:
// "to me, it suggests a likely interaction with bgi and/or the need for
// multiple expressions of bgi. all of these equations are highly nonlinear."
//
// v8 modifications relative to v7:
//
//   1. BGI basis expansion: instead of bgi + bgi^2, use a 3-piece basis
//      with two interior knots. This lets the response curve be flat at
//      low BGI, accelerate in the middle, and saturate at high BGI without
//      forcing a global parabola.
//
//        bgi_b1 = bgi
//        bgi_b2 = (bgi - k1)_+   (positive part, kink at knot k1)
//        bgi_b3 = (bgi - k2)_+   (positive part, kink at knot k2)
//
//      Knots passed in as data so we can place them at empirical quantiles
//      (e.g. 25th and 75th percentile of bgi).
//
//   2. BGI x ln_DBH interaction. Size-dependent climate response —
//      large trees may have a different relationship with BGI than
//      saplings (deep root systems, light competition, etc.).
//
//   3. BGI x softwood interaction. Functional-group climate response;
//      conifers and broadleaves respond differently to growing season
//      length and biomass-relevant climate.
//
//   4. BGI x ln_cr_adj interaction. Vigorous trees (high crown ratio)
//      can capitalize on better climate; suppressed trees can't.
//
//   5. Keep BGI x rd_additive interaction from v7.
//
//   6. Keep z_L1_bgi (ecoregion-varying linear bgi slope).
//   7. Keep gamma_site (trait-modulated linear bgi slope).
//
// Linear predictor:
//   eta = b0 + trait_effect[sp] + z_sp[sp]
//       + z_L1 + z_L2 + z_L3 + z_FT
//       + b1 ln(DBH) + b2 DBH + b3 ln_cr_adj
//       + b4 ln_bal_sw_adj + b5 bal_hw
//       + b_site * bgi_b1 + b9a bgi_b2 + b9b bgi_b3   <-- 3-piece BGI basis
//       + b7  BA  * rd_additive
//       + b8  BAL_SW * rd_additive
//       + b11 sdi_complexity
//       + b12 bgi * rd_additive
//       + b13 bgi * ln_dbh                            <-- size x climate
//       + b14 bgi * softwood                          <-- functional group x climate
//       + b15 bgi * ln_cr_adj                         <-- vigor x climate
//
// where:
//   b_site_i = b6 + z_L1_bgi[L1_i] + W[sp_i,] * gamma_site
//   softwood comes from per-tree trait join (passed as data)
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
  vector[N_obs] softwood;          // 0/1 per tree (from trait join)

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;

  // Knot locations on bgi scale (pass empirical quantiles from R driver)
  real bgi_knot1;
  real bgi_knot2;
}
transformed data {
  // 3-piece BGI basis: linear + two positive-part hinge terms
  vector[N_obs] bgi_b1 = bgi;
  vector[N_obs] bgi_b2;
  vector[N_obs] bgi_b3;
  for (i in 1:N_obs) {
    bgi_b2[i] = fmax(bgi[i] - bgi_knot1, 0.0);
    bgi_b3[i] = fmax(bgi[i] - bgi_knot2, 0.0);
  }

  // Interaction terms
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
  real b6;                          // bgi linear (b_site mean)
  real b7;                          // BA x rd_additive
  real b8;                          // BAL_SW x rd_additive
  real b9a;                         // hinge above bgi_knot1
  real b9b;                         // hinge above bgi_knot2
  real b11;                         // sdi_complexity
  real b12;                         // bgi x rd_additive
  real b13;                         // bgi x ln_dbh
  real b14;                         // bgi x softwood
  real b15;                         // bgi x ln_cr_adj

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
  b6  ~ normal( 0.05, 0.10);   // bgi linear; allow either sign but weakly positive
  b7  ~ normal( 0.0, 0.02);
  b8  ~ normal( 0.0, 0.02);
  b9a ~ normal( 0.0, 0.10);    // hinges allowed in either direction
  b9b ~ normal( 0.0, 0.10);
  b11 ~ normal( 0.0, 0.3);
  b12 ~ normal( 0.0, 0.1);
  b13 ~ normal( 0.0, 0.05);    // bgi x ln_dbh
  b14 ~ normal( 0.0, 0.10);    // bgi x softwood
  b15 ~ normal( 0.0, 0.10);    // bgi x ln_cr_adj

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
