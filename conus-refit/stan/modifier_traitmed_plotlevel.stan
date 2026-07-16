// ============================================================================
// modifier_traitmed_plotlevel.stan
//
// Trait-mediated species-specific modifiers for PLOT-LEVEL residuals
// (mortality, where Bernoulli per-tree residuals don't work).
//
// Approach: plot-level residual already aggregates across species.
// We use a plot-level dominant-species trait vector W_plot[plot, ]
// passed as data (the most common species in the plot, or trait-
// weighted by tree count).
//
// Linear predictor for plot-level cloglog residual:
//   delta_plot = alpha_0
//              + (alpha_plant + W_plot * gamma_plant) * is_plantation
//              + (alpha_fire  + W_plot * gamma_fire ) * x_fire
//              + ... etc
//              + z_L1[L1_idx]
//
// W_plot is N_obs x P_trait (already trait-weighted per plot).
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;
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

  // Plot-level trait vector (already aggregated by tree-count weighting)
  matrix[N_obs, P_trait > 0 ? P_trait : 1] W_plot;
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
  vector[N_obs] tm_plant_obs    = W_plot * gamma_alpha_plant;
  vector[N_obs] tm_fire_obs     = W_plot * gamma_alpha_fire;
  vector[N_obs] tm_insect_obs   = W_plot * gamma_alpha_insect;
  vector[N_obs] tm_disease_obs  = W_plot * gamma_alpha_disease;
  vector[N_obs] tm_wind_obs     = W_plot * gamma_alpha_wind;
  vector[N_obs] tm_harvest_obs  = W_plot * gamma_alpha_harvest;
  vector[N_obs] tm_cutting_obs  = W_plot * gamma_alpha_cutting;
  vector[N_obs] tm_siteprep_obs = W_plot * gamma_alpha_siteprep;
}
model {
  alpha_0        ~ normal(0, 0.5);
  alpha_plant    ~ normal(0, 1.0);
  alpha_fire     ~ normal(0, 1.0);
  alpha_insect   ~ normal(0, 1.0);
  alpha_disease  ~ normal(0, 1.0);
  alpha_wind     ~ normal(0, 1.0);
  alpha_harvest  ~ normal(0, 1.0);
  alpha_cutting  ~ normal(0, 1.0);
  alpha_siteprep ~ normal(0, 1.0);

  gamma_alpha_plant    ~ normal(0, 0.3);
  gamma_alpha_fire     ~ normal(0, 0.3);
  gamma_alpha_insect   ~ normal(0, 0.3);
  gamma_alpha_disease  ~ normal(0, 0.3);
  gamma_alpha_wind     ~ normal(0, 0.3);
  gamma_alpha_harvest  ~ normal(0, 0.3);
  gamma_alpha_cutting  ~ normal(0, 0.3);
  gamma_alpha_siteprep ~ normal(0, 0.3);

  z_L1_raw    ~ std_normal();
  sigma_L1    ~ normal(0, 0.5);
  sigma_resid ~ normal(0, 2.0);

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
  vector[N_obs] log_lik;
  {
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
    for (i in 1:N_obs) {
      log_lik[i] = normal_lpdf(residual[i] | delta[i],
                               sigma_resid * inv_weight_safe[i]);
    }
  }
}
