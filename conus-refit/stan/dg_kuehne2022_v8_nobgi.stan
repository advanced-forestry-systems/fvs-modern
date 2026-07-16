// ============================================================================
// dg_kuehne2022_v8_nobgi.stan
//
// DG v8 with ALL BGI terms removed. Productivity is carried only by the
// ecoregion random effects (z_L1/L2/L3, keyed by mappable location), species,
// and stand structure. This is the "BGI neutralized for projection" model.
//
// elpd(v8 with BGI) - elpd(this) = the full predictive value of BGI, i.e. the
// penalty of neutralizing BGI at projection time (where BGI is unavailable;
// it is unmappable R^2 0.001 and not stand-derivable R^2 0.016).
//
// Data block is identical to dg_kuehne2022_v8_bgi_nonlinear.stan so the driver
// feeds it unchanged; the bgi/knot fields are simply unused here.
//
// Linear predictor:
//   eta = b0 + trait_effect[sp] + z_sp[sp] + z_L1 + z_L2 + z_L3 + z_FT
//       + b1 ln(DBH) + b2 DBH + b3 ln_cr_adj
//       + b4 ln_bal_sw_adj + b5 bal_hw
//       + b7 BA*rd_additive + b8 BAL_SW*rd_additive + b11 sdi_complexity
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

  vector[N_obs] bgi;               // unused (kept for driver compatibility)
  vector[N_obs] ba_metric;
  vector[N_obs] bal_sw_metric;
  vector[N_obs] rd_additive;
  vector[N_obs] sdi_complexity;
  vector[N_obs] softwood;          // unused

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;

  real bgi_knot1;                  // unused
  real bgi_knot2;                  // unused
}
transformed data {
  vector[N_obs] ba_x_rdadd    = ba_metric     .* rd_additive;
  vector[N_obs] balsw_x_rdadd = bal_sw_metric .* rd_additive;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5;
  real b7;                          // BA x rd_additive
  real b8;                          // BAL_SW x rd_additive
  real b11;                         // sdi_complexity

  vector[P_trait] gamma;

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;

  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;

  real<lower=0> sigma;
}
transformed parameters {
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
  vector[N_FT] z_FT = sigma_FT * z_FT_raw;
}
model {
  b0  ~ normal(-1.0, 2.0);
  b1  ~ normal( 0.3, 0.5);
  b2  ~ normal(-0.02, 0.05);
  b3  ~ normal( 0.5, 0.5);
  b4  ~ normal(-0.05, 0.1);
  b5  ~ normal(-0.005, 0.02);
  b7  ~ normal( 0.0, 0.02);
  b8  ~ normal( 0.0, 0.02);
  b11 ~ normal( 0.0, 0.3);

  gamma ~ normal(0, 0.5);

  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();

  sigma_sp ~ normal(0, 0.15);
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);

  sigma ~ normal(0, 0.5);

  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b7 * ba_x_rdadd + b8 * balsw_x_rdadd
    + b11 * sdi_complexity;

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
  vector[N_obs] eta_gq =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b7 * ba_x_rdadd + b8 * balsw_x_rdadd
    + b11 * sdi_complexity;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
