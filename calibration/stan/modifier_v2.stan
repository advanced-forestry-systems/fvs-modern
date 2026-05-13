// ============================================================================
// modifier_v2.stan
//
// Continuous-time successor to modifier_common.stan. Consumes the new
// kernel columns produced by 40c_enrich_cond_modifiers_v2.R:
//
//   acute_kernel_5yr / acute_kernel_10yr  for fire, wind, harvest
//   gamma_kernel_3yr / gamma_kernel_5yr   for insect, disease, treatment
//   plant_age_taper_30, plant_age_taper_50  for plantation main effect
//
// The choice of kernel per event type is encoded in the Stan data:
// each event indicator multiplies its appropriate kernel column. This
// keeps the Stan code clean and lets us swap kernel choice via the
// data side without recompiling.
//
// Linear predictor for residual:
//   delta = alpha_0
//         + alpha_plant * is_plantation * plant_age_taper
//         + alpha_fire     * d_fire     * acute_kernel
//         + alpha_wind     * d_wind     * acute_kernel
//         + alpha_harvest  * d_harvest  * acute_kernel
//         + alpha_insect   * d_insect   * gamma_kernel_dstrb
//         + alpha_disease  * d_disease  * gamma_kernel_dstrb
//         + alpha_cutting  * t_cutting  * gamma_kernel_trt
//         + alpha_siteprep * t_site_prep * acute_kernel_trt
//         + z_L1[L1_idx]
// ============================================================================

data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;

  vector[N_obs] residual;
  vector[N_obs] weight;

  // Plantation
  vector[N_obs] is_plantation;
  vector[N_obs] plant_age_taper;          // sigmoid(STDAGE), passed in by R

  // Disturbance type indicators
  vector[N_obs] d_fire;
  vector[N_obs] d_wind;
  vector[N_obs] d_harvest;
  vector[N_obs] d_insect;
  vector[N_obs] d_disease;
  vector[N_obs] acute_kernel_dstrb;        // exp(-t/tau) for acute events
  vector[N_obs] gamma_kernel_dstrb;        // (t/tau)*exp(1-t/tau) for delayed

  // Treatment type indicators
  vector[N_obs] t_cutting;
  vector[N_obs] t_site_prep;
  vector[N_obs] acute_kernel_trt;
  vector[N_obs] gamma_kernel_trt;

  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
}

transformed data {
  // Pre-multiply indicators by their appropriate kernel
  vector[N_obs] x_plant    = is_plantation .* plant_age_taper;
  vector[N_obs] x_fire     = d_fire        .* acute_kernel_dstrb;
  vector[N_obs] x_wind     = d_wind        .* acute_kernel_dstrb;
  vector[N_obs] x_harvest  = d_harvest     .* acute_kernel_dstrb;
  vector[N_obs] x_insect   = d_insect      .* gamma_kernel_dstrb;
  vector[N_obs] x_disease  = d_disease     .* gamma_kernel_dstrb;
  vector[N_obs] x_cutting  = t_cutting     .* gamma_kernel_trt;
  vector[N_obs] x_siteprep = t_site_prep   .* acute_kernel_trt;

  vector[N_obs] inv_weight_safe;
  for (i in 1:N_obs) inv_weight_safe[i] = 1.0 / fmax(weight[i], 1e-4);
}

parameters {
  real alpha_0;
  real alpha_plant;
  real alpha_fire;
  real alpha_wind;
  real alpha_harvest;
  real alpha_insect;
  real alpha_disease;
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
  alpha_wind     ~ normal(0, 0.2);
  alpha_harvest  ~ normal(0, 0.2);
  alpha_insect   ~ normal(0, 0.2);
  alpha_disease  ~ normal(0, 0.2);
  alpha_cutting  ~ normal(0, 0.2);
  alpha_siteprep ~ normal(0, 0.2);

  z_L1_raw    ~ std_normal();
  sigma_L1    ~ normal(0, 0.1);
  sigma_resid ~ normal(0, 0.3);

  vector[N_obs] delta =
      alpha_0
    + alpha_plant    * x_plant
    + alpha_fire     * x_fire
    + alpha_wind     * x_wind
    + alpha_harvest  * x_harvest
    + alpha_insect   * x_insect
    + alpha_disease  * x_disease
    + alpha_cutting  * x_cutting
    + alpha_siteprep * x_siteprep
    + z_L1[L1_idx];

  residual ~ normal(delta, sigma_resid * inv_weight_safe);
}
