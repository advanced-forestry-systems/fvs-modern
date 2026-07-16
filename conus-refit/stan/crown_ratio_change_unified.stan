// ============================================================================
// crown_ratio_change_unified.stan
//
// Convergence fix for crown_ratio_change_speciesfree.stan. The original 100k
// fits did not converge (max rhat 3.5 to 3.9, all parameters), for two reasons:
//   1. CR was the ONLY component without a species random intercept; all species
//      variation was forced through trait_effect = W*gamma, which trades off
//      against b0 and the ecoregion/FT REs.
//   2. With a near-zero-signal response (sigma ~ 0.023) and four fully additive
//      intercept blocks (b0 + z_L1 + z_L2 + z_L3 + z_FT + trait_effect), only the
//      SUM is identified, not the individual levels -> an intercept ridge.
//
// Fixes (consistent with the unified architecture used by every other component):
//   - add z_sp / sigma_sp species random intercept;
//   - soft sum-to-zero on each RE block so b0 uniquely carries the grand mean;
//   - tighten b0 to the response scale.
// Predictions are unchanged in expectation; this only makes the model identified.
//
// Linear predictor:
//   delta_CR = b0 + trait_effect[sp] + z_sp[sp]
//            + z_L1 + z_L2 + z_L3 + z_FT
//            + b1 DBH + b2 DBH^2 + b3 BA + b4 BAL + b5 CR_init + b6 ln_csi
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector[N_obs] delta_CR_a;
  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;
  vector[N_obs] ba_metric;
  vector[N_obs] bal_metric;
  vector[N_obs] cr_init;
  vector[N_obs] ln_csi;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5; real b6;

  vector[P_trait] gamma;

  vector[N_sp] z_sp_raw;            // NEW species random intercept
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;

  real<lower=0> sigma_sp;           // NEW species RE scale
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
  b0 ~ normal( 0.0, 0.25);    // tightened to the response scale
  b1 ~ normal( 0.0, 0.1);
  b2 ~ normal( 0.0, 0.01);
  b3 ~ normal( 0.0, 0.05);
  b4 ~ normal( 0.0, 0.05);
  b5 ~ normal(-0.1, 0.5);
  b6 ~ normal( 0.0, 0.5);

  gamma ~ normal(0, 0.5);

  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();

  // Soft sum-to-zero: pin each RE block's level so b0 is the unique grand mean.
  sum(z_sp_raw) ~ normal(0, 0.01);
  sum(z_L1_raw) ~ normal(0, 0.01);
  sum(z_L2_raw) ~ normal(0, 0.01);
  sum(z_L3_raw) ~ normal(0, 0.01);
  sum(z_FT_raw) ~ normal(0, 0.01);

  sigma_sp ~ normal(0, 0.3);
  sigma_L1 ~ normal(0, 0.3);
  sigma_L2 ~ normal(0, 0.2);
  sigma_L3 ~ normal(0, 0.2);
  sigma_FT ~ normal(0, 0.2);

  sigma ~ normal(0, 0.5);

  vector[N_obs] mu =
      b0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * ba_metric + b4 * bal_metric
    + b5 * cr_init
    + b6 * ln_csi;

  delta_CR_a ~ normal(mu, sigma);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_pred =
      b0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * ba_metric + b4 * bal_metric
    + b5 * cr_init
    + b6 * ln_csi;
  for (i in 1:N_obs) {
    log_lik[i] = normal_lpdf(delta_CR_a[i] | mu_pred[i], sigma);
  }
}
