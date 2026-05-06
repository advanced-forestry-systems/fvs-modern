// ============================================================================
// modifier_binary.stan
//
// FVS-CONUS Phase 3 residual-fit modifier model for BINARY base
// outcomes (mortality / survival). Unlike modifier_common.stan which
// regresses a precomputed continuous residual on modifier indicators
// and L1 RE, this model carries the binary outcome through a Bernoulli
// logit likelihood with the base-fit eta as an offset. That preserves
// the correct binomial uncertainty rather than approximating
// logit(observed_indicator) as a continuous residual (logit(0) and
// logit(1) are undefined).
//
// Linear predictor on the logit scale:
//   logit(p_i) = eta_base_i + log_years_i + delta_i
//   delta_i    = alpha_plant * is_plantation_i
//              + sum_k alpha_dstrb_k * d_k_i * dstrb_decay_i
//              + sum_j alpha_trt_j   * t_j_i * trt_decay_i
//              + z_L1[l1_i]
//
// alive_i ~ Bernoulli(p_i)
//
// eta_base is the precomputed posterior-mean linear predictor from the
// base mortality fit (without modifiers and without the log_years
// offset). log_years is added separately so this model can also be
// applied to base fits that absorbed log_years into eta_base; in that
// case pass log_years = 0.
// ============================================================================

data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;

  array[N_obs] int<lower=0, upper=1> alive;
  vector[N_obs] eta_base;
  vector[N_obs] log_years;

  // Stand origin
  vector[N_obs] is_plantation;

  // Disturbance indicators x decay envelope
  vector[N_obs] d_fire;
  vector[N_obs] d_insect;
  vector[N_obs] d_disease;
  vector[N_obs] d_wind;
  vector[N_obs] d_harvest;
  vector[N_obs] dstrb_decay;

  // Treatment indicators x decay envelope
  vector[N_obs] t_cutting;
  vector[N_obs] t_site_prep;
  vector[N_obs] trt_decay;

  // EPA L1 random effect key
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
}

transformed parameters {
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
}

model {
  // Priors mirror modifier_common.stan: tight around zero on the
  // logit scale so effects only emerge when the data demand them.
  alpha_0        ~ normal(0, 0.1);
  alpha_plant    ~ normal(0, 0.3);
  alpha_fire     ~ normal(0, 0.2);
  alpha_insect   ~ normal(0, 0.2);
  alpha_disease  ~ normal(0, 0.2);
  alpha_wind     ~ normal(0, 0.2);
  alpha_harvest  ~ normal(0, 0.2);
  alpha_cutting  ~ normal(0, 0.2);
  alpha_siteprep ~ normal(0, 0.2);

  z_L1_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.1);

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

  vector[N_obs] eta_total = eta_base + log_years + delta;

  // Clamp eta_total to [-30, +30] to avoid logit-scale overflow during
  // warmup; mirrors the HG/DG eta-clamp pattern.
  vector[N_obs] eta_clamped;
  for (i in 1:N_obs) eta_clamped[i] = fmin(fmax(eta_total[i], -30.0), 30.0);

  alive ~ bernoulli_logit(eta_clamped);
}

generated quantities {
  vector[N_obs] delta_hat =
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
}
