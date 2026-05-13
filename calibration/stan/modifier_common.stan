// ============================================================================
// modifier_common.stan
//
// FVS-CONUS Phase 3 residual-fit modifier model. Takes a precomputed
// per-observation residual (base-model response minus base-model eta,
// on whatever scale the base model uses) and regresses it on stand-
// origin, disturbance, and treatment indicators with exponential
// time-since decay envelopes and an EPA Level 1 random effect.
//
// This model is intentionally agnostic about the base model that
// produced the residuals. For a log-response base model (DG, HG),
// residual = log(obs) - eta_base. For a logit-response base model
// (HCB, mortality, CR), residual = logit(obs) - eta_base. The
// modifier coefficients are interpreted on the same scale as the base
// eta, so plugging them back into the projection engine is a simple
// additive correction: eta_total = eta_base + delta.
//
// Linear predictor for the residual:
//   delta_i = alpha_plant * is_plantation_i
//           + sum_k  alpha_dstrb_k * d_k_i * dstrb_decay_i
//           + sum_j  alpha_trt_j   * t_j_i * trt_decay_i
//           + z_L1[l1_i]
//
// where d_k_i is the indicator for disturbance type k (fire, insect,
// disease, wind, harvest) and t_j_i for treatment type j (cutting,
// site prep). Decay envelopes are passed in as data; fitting the
// lambdas themselves is deferred to the integrated Phase 4 re-fit.
// ============================================================================

data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;

  vector[N_obs] residual;
  vector[N_obs] weight;                 // weight_i (e.g., 1 / sd_base_i) or all 1

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
  // Precompute the interactions once
  vector[N_obs] x_fire     = d_fire     .* dstrb_decay;
  vector[N_obs] x_insect   = d_insect   .* dstrb_decay;
  vector[N_obs] x_disease  = d_disease  .* dstrb_decay;
  vector[N_obs] x_wind     = d_wind     .* dstrb_decay;
  vector[N_obs] x_harvest  = d_harvest  .* dstrb_decay;
  vector[N_obs] x_cutting  = t_cutting  .* trt_decay;
  vector[N_obs] x_siteprep = t_site_prep .* trt_decay;

  // Precompute clamped inverse weight. weight should be a small
  // positive number; fmax guards against zero / negative input.
  vector[N_obs] inv_weight_safe;
  for (i in 1:N_obs) inv_weight_safe[i] = 1.0 / fmax(weight[i], 1e-4);
}

parameters {
  real alpha_0;                 // intercept (residual bias term)
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
  real<lower=1e-3> sigma_resid;   // strict lower bound keeps the normal scale finite
}

transformed parameters {
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
}

model {
  // Priors: tight around zero so the modifier pulls effect only when
  // the data demand it. These are on the same scale as eta_base.
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

  // Vectorized scale using precomputed inverse weight. sigma_resid is
  // already bounded below in the parameters block so the product is
  // strictly positive.
  residual ~ normal(delta, sigma_resid * inv_weight_safe);
}

generated quantities {
  // Posterior-mean modifier per observation for diagnostic purposes
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
