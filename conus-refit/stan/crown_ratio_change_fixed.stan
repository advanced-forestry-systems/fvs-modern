// ============================================================================
// crown_ratio_change_fixed.stan
//
// Parsimonious FIXED-EFFECTS crown-recession model. Both hierarchical CR fits
// (rhat 3.5 to 3.96) failed because the near-zero-signal response cannot
// identify five nested RE variances (sigma_sp/L1/L2/L3/FT all collapsed to ~0
// with ESS 4). This drops all random effects; species variation is carried by
// the trait fixed effects (W * gamma). This is identified and will converge.
//
//   delta_CR = b0 + trait_effect[sp] + b1 DBH + b2 DBH^2
//            + b3 BA + b4 BAL + b5 CR_init + b6 ln_csi
//
// Extra data fields the driver passes (RE indices, N_L*) are simply not declared
// here and are ignored.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=0> P_trait;

  vector[N_obs] delta_CR_a;
  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;
  vector[N_obs] ba_metric;
  vector[N_obs] bal_metric;
  vector[N_obs] cr_init;
  vector[N_obs] ln_csi;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5; real b6;
  vector[P_trait] gamma;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N_sp] trait_effect;
  if (P_trait > 0) trait_effect = W * gamma;
  else trait_effect = rep_vector(0.0, N_sp);
}
model {
  b0 ~ normal( 0.0, 0.25);
  b1 ~ normal( 0.0, 0.1);
  b2 ~ normal( 0.0, 0.01);
  b3 ~ normal( 0.0, 0.05);
  b4 ~ normal( 0.0, 0.05);
  b5 ~ normal(-0.1, 0.5);
  b6 ~ normal( 0.0, 0.5);
  gamma ~ normal(0, 0.5);
  sigma ~ normal(0, 0.5);

  vector[N_obs] mu =
      b0 + trait_effect[sp_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * ba_metric + b4 * bal_metric
    + b5 * cr_init + b6 * ln_csi;
  delta_CR_a ~ normal(mu, sigma);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_pred =
      b0 + trait_effect[sp_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * ba_metric + b4 * bal_metric
    + b5 * cr_init + b6 * ln_csi;
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(delta_CR_a[i] | mu_pred[i], sigma);
}
