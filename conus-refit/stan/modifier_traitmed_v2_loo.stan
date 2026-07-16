// ============================================================================
// modifier_traitmed_v2_loo.stan
//
// modifier_traitmed.stan + log_lik in generated quantities for formal LOO
// comparison against modifier_common_v2_loo.stan. Identical model;
// addition of log_lik unlocks ELPD comparison.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;
  int<lower=1> N_sp;
  int<lower=0> P_trait;

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
  array[N_obs] int<lower=1, upper=N_sp> sp_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
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

  vector[P_trait] gamma_alpha_plant;
  vector[P_trait] gamma_alpha_fire;
  vector[P_trait] gamma_alpha_insect;
  vector[P_trait] gamma_alpha_disease;
  vector[P_trait] gamma_alpha_wind;
  vector[P_trait] gamma_alpha_harvest;
  vector[P_trait] gamma_alpha_cutting;
  vector[P_trait] gamma_alpha_siteprep;

  vector[N_L1] z_L1_raw;
  real<lower=1e-4> sigma_L1;
  real<lower=1e-3> sigma_resid;
}
transformed parameters {
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_sp] tm_plant    = W * gamma_alpha_plant;
  vector[N_sp] tm_fire     = W * gamma_alpha_fire;
  vector[N_sp] tm_insect   = W * gamma_alpha_insect;
  vector[N_sp] tm_disease  = W * gamma_alpha_disease;
  vector[N_sp] tm_wind     = W * gamma_alpha_wind;
  vector[N_sp] tm_harvest  = W * gamma_alpha_harvest;
  vector[N_sp] tm_cutting  = W * gamma_alpha_cutting;
  vector[N_sp] tm_siteprep = W * gamma_alpha_siteprep;
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

  gamma_alpha_plant    ~ normal(0, 0.15);
  gamma_alpha_fire     ~ normal(0, 0.10);
  gamma_alpha_insect   ~ normal(0, 0.10);
  gamma_alpha_disease  ~ normal(0, 0.10);
  gamma_alpha_wind     ~ normal(0, 0.10);
  gamma_alpha_harvest  ~ normal(0, 0.10);
  gamma_alpha_cutting  ~ normal(0, 0.10);
  gamma_alpha_siteprep ~ normal(0, 0.10);

  z_L1_raw    ~ std_normal();
  sigma_L1    ~ normal(0, 0.1);
  sigma_resid ~ normal(0, 0.3);

  vector[N_obs] tm_plant_obs;
  vector[N_obs] tm_fire_obs;
  vector[N_obs] tm_insect_obs;
  vector[N_obs] tm_disease_obs;
  vector[N_obs] tm_wind_obs;
  vector[N_obs] tm_harvest_obs;
  vector[N_obs] tm_cutting_obs;
  vector[N_obs] tm_siteprep_obs;
  for (i in 1:N_obs) {
    int s_i = sp_idx[i];
    tm_plant_obs[i]    = tm_plant[s_i];
    tm_fire_obs[i]     = tm_fire[s_i];
    tm_insect_obs[i]   = tm_insect[s_i];
    tm_disease_obs[i]  = tm_disease[s_i];
    tm_wind_obs[i]     = tm_wind[s_i];
    tm_harvest_obs[i]  = tm_harvest[s_i];
    tm_cutting_obs[i]  = tm_cutting[s_i];
    tm_siteprep_obs[i] = tm_siteprep[s_i];
  }
  vector[N_obs] delta =
      alpha_0
    + (alpha_plant    + tm_plant_obs)    .* is_plantation
    + (alpha_fire     + tm_fire_obs)     .* x_fire
    + (alpha_insect   + tm_insect_obs)   .* x_insect
    + (alpha_disease  + tm_disease_obs)  .* x_disease
    + (alpha_wind     + tm_wind_obs)     .* x_wind
    + (alpha_harvest  + tm_harvest_obs)  .* x_harvest
    + (alpha_cutting  + tm_cutting_obs)  .* x_cutting
    + (alpha_siteprep + tm_siteprep_obs) .* x_siteprep
    + z_L1[L1_idx];

  residual ~ normal(delta, sigma_resid * inv_weight_safe);
}
generated quantities {
  vector[N_sp] alpha_plant_sp    = alpha_plant    + tm_plant;
  vector[N_sp] alpha_insect_sp   = alpha_insect   + tm_insect;
  vector[N_sp] alpha_disease_sp  = alpha_disease  + tm_disease;
  vector[N_sp] alpha_cutting_sp  = alpha_cutting  + tm_cutting;
  vector[N_sp] alpha_harvest_sp  = alpha_harvest  + tm_harvest;

  vector[N_obs] log_lik;
  {
    vector[N_obs] delta =
        alpha_0
      + (alpha_plant    + tm_plant[sp_idx])    .* is_plantation
      + (alpha_fire     + tm_fire[sp_idx])     .* x_fire
      + (alpha_insect   + tm_insect[sp_idx])   .* x_insect
      + (alpha_disease  + tm_disease[sp_idx])  .* x_disease
      + (alpha_wind     + tm_wind[sp_idx])     .* x_wind
      + (alpha_harvest  + tm_harvest[sp_idx])  .* x_harvest
      + (alpha_cutting  + tm_cutting[sp_idx])  .* x_cutting
      + (alpha_siteprep + tm_siteprep[sp_idx]) .* x_siteprep
      + z_L1[L1_idx];
    for (i in 1:N_obs) {
      log_lik[i] = normal_lpdf(residual[i] | delta[i],
                               sigma_resid * inv_weight_safe[i]);
    }
  }
}
