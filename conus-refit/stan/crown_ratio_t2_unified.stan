// ============================================================================
// crown_ratio_t2_unified.stan
//
// CROWN-RATIO REDESIGN. Modeling the change (delta_CR) failed across three
// architectures (rhat 3.5-3.96) because the differenced, FIA-class-discretized
// response is near-degenerate. This models crown ratio at time 2 (CR2) DIRECTLY
// on the logit scale, with logit(CR1) as the dominant predictor. CR2 has real
// variance, so the unified hierarchical architecture is identifiable here.
//
//   logit(CR2) = b0 + trait_effect[sp] + z_sp[sp]
//              + z_L1 + z_L2 + z_L3 + z_FT
//              + b_cr1 logit(CR1)
//              + b1 DBH + b2 DBH^2 + b3 BA + b4 BAL + b6 ln_csi
//
// Implied recession = CR1 - inv_logit(predicted logit(CR2)).
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector[N_obs] cr2_logit;     // response: logit(CR2)
  vector[N_obs] cr1_logit;     // predictor: logit(CR1)
  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;
  vector[N_obs] ba_metric;
  vector[N_obs] bal_metric;
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
  real b_cr1;
  real b1; real b2; real b3; real b4; real b6;
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
  if (P_trait > 0) trait_effect = W * gamma; else trait_effect = rep_vector(0.0, N_sp);
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
  vector[N_FT] z_FT = sigma_FT * z_FT_raw;
}
model {
  b0    ~ normal(0, 1.0);
  b_cr1 ~ normal(0.8, 0.5);   // CR2 tracks CR1 strongly and positively
  b1 ~ normal(0, 0.1);
  b2 ~ normal(0, 0.01);
  b3 ~ normal(0, 0.05);
  b4 ~ normal(0, 0.05);
  b6 ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);

  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.3);
  sigma_L1 ~ normal(0, 0.3);
  sigma_L2 ~ normal(0, 0.2);
  sigma_L3 ~ normal(0, 0.2);
  sigma_FT ~ normal(0, 0.2);
  sigma ~ normal(0, 1.0);

  vector[N_obs] mu =
      b0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b_cr1 * cr1_logit
    + b1 * dbh + b2 * dbh_sq + b3 * ba_metric + b4 * bal_metric + b6 * ln_csi;
  cr2_logit ~ normal(mu, sigma);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_pred =
      b0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b_cr1 * cr1_logit
    + b1 * dbh + b2 * dbh_sq + b3 * ba_metric + b4 * bal_metric + b6 * ln_csi;
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(cr2_logit[i] | mu_pred[i], sigma);
}
