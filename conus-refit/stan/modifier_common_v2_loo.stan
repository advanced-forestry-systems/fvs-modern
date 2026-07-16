// ============================================================================
// modifier_common_v2_loo.stan
//
// modifier_common.stan + log_lik in generated quantities for LOO comparison.
// Identical model; addition of log_lik unlocks lambda selection via ELPD.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;

  vector[N_obs] residual;
  vector[N_obs] weight;

  vector[N_obs] is_plantation;

  vector[N_obs] d_fire;
  vector[N_obs] d_insect;
  vector[N_obs] d_disease;
  vector[N_obs] d_wind;
  vector[N_obs] d_harvest;
  vector[N_obs] dstrb_decay;

  vector[N_obs] t_cutting;
  vector[N_obs] t_site_prep;
  vector[N_obs] trt_decay;

  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
}
transformed data {
  vector[N_obs] x_fire     = d_fire     .* dstrb_decay;
  vector[N_obs] x_insect   = d_insect   .* dstrb_decay;
  vector[N_obs] x_disease  = d_disease  .* dstrb_decay;
  vector[N_obs] x_wind     = d_wind     .* dstrb_decay;
  vector[N_obs] x_harvest  = d_harvest  .* dstrb_decay;
  vector[N_obs] x_cutting  = t_cutting  .* trt_decay;
  vector[N_obs] x_siteprep = t_site_prep .* trt_decay;

  vector[N_obs] inv_weight_safe;
  for (i in 1:N_obs) inv_weight_safe[i] = 1.0 / fmax(weight[i], 1e-4);
}
parameters {
  real alpha_0;
  real alpha_plant;
  real alpha_fire;
  real alpha_insect;
  real alpha_disease;
  real alpha_wind;
  real alpha_harvest;
  real alpha_cutting;
  real alpha_siteprep;

  vector[N_L1] z_L1_raw;
  real<lower=1e-4> sigma_L1;
  real<lower=1e-3> sigma_resid;
}
transformed parameters {
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
}
model {
  alpha_0        ~ normal(0, 0.1);
  alpha_plant    ~ normal(0, 0.3);
  alpha_fire     ~ normal(0, 0.2);
  alpha_insect   ~ normal(0, 0.2);
  alpha_disease  ~ normal(0, 0.2);
  alpha_wind     ~ normal(0, 0.2);
  alpha_harvest  ~ normal(0, 0.2);
  alpha_cutting  ~ normal(0, 0.2);
  alpha_siteprep ~ normal(0, 0.2);

  z_L1_raw    ~ std_normal();
  sigma_L1    ~ normal(0, 0.1);
  sigma_resid ~ normal(0, 0.3);

  vector[N_obs] delta =
      alpha_0
    + alpha_plant    * is_plantation
    + alpha_fire     * x_fire
    + alpha_insect   * x_insect
    + alpha_disease  * x_disease
    + alpha_wind     * x_wind
    + alpha_harvest  * x_harvest
    + alpha_cutting  * x_cutting
    + alpha_siteprep * x_siteprep
    + z_L1[L1_idx];

  residual ~ normal(delta, sigma_resid * inv_weight_safe);
}
generated quantities {
  vector[N_obs] delta_hat;
  vector[N_obs] log_lik;
  {
    vector[N_obs] delta =
        alpha_0
      + alpha_plant    * is_plantation
      + alpha_fire     * x_fire
      + alpha_insect   * x_insect
      + alpha_disease  * x_disease
      + alpha_wind     * x_wind
      + alpha_harvest  * x_harvest
      + alpha_cutting  * x_cutting
      + alpha_siteprep * x_siteprep
      + z_L1[L1_idx];
    delta_hat = delta;
    for (i in 1:N_obs) {
      log_lik[i] = normal_lpdf(residual[i] | delta[i],
                               sigma_resid * inv_weight_safe[i]);
    }
  }
}
