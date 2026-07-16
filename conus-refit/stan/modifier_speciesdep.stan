// ============================================================================
// modifier_speciesdep.stan
//
// SPECIES-DEPENDENT (b2-analog) extension of modifier_traitmed.stan. Adds a
// per-species random intercept on each disturbance/treatment alpha with a
// TRAIT-INFORMED prior mean, so the modifier framework gets the same
// species-free (b1) vs species-specific (b2) treatment as the base equations:
//
//   alpha_k_i = alpha_k_global + W[sp_i,]*gamma_k + z_sp_k[sp_i]
//   z_sp_k[s] ~ Normal(0, sigma_sp_k)        (free per-species deviation)
//
// The trait part W*gamma_k is the informed mean (carries rare/unsampled
// species); z_sp_k is the per-species deviation that well-sampled species can
// take. Setting sigma_sp_k -> 0 recovers modifier_traitmed (species-free), so
// the two are directly comparable by held-out ELPD. Deploy via a per-species
// shrinkage weight w = n/(n+kappa) blending this against the species-free leg.
//
// Same residual data as modifier_traitmed.stan (no new inputs).
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

  // NEW: per-species deviations on each alpha (non-centered)
  vector[N_sp] zsp_plant_raw;
  vector[N_sp] zsp_fire_raw;
  vector[N_sp] zsp_insect_raw;
  vector[N_sp] zsp_disease_raw;
  vector[N_sp] zsp_wind_raw;
  vector[N_sp] zsp_harvest_raw;
  vector[N_sp] zsp_cutting_raw;
  vector[N_sp] zsp_siteprep_raw;
  real<lower=1e-4> sigma_sp_plant;
  real<lower=1e-4> sigma_sp_fire;
  real<lower=1e-4> sigma_sp_insect;
  real<lower=1e-4> sigma_sp_disease;
  real<lower=1e-4> sigma_sp_wind;
  real<lower=1e-4> sigma_sp_harvest;
  real<lower=1e-4> sigma_sp_cutting;
  real<lower=1e-4> sigma_sp_siteprep;

  vector[N_L1] z_L1_raw;
  real<lower=1e-4> sigma_L1;
  real<lower=1e-3> sigma_resid;
}

transformed parameters {
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;

  // species effective shift = trait mean (W*gamma) + free deviation (sigma*zsp)
  vector[N_sp] sp_plant    = W * gamma_alpha_plant    + sigma_sp_plant    * zsp_plant_raw;
  vector[N_sp] sp_fire     = W * gamma_alpha_fire     + sigma_sp_fire     * zsp_fire_raw;
  vector[N_sp] sp_insect   = W * gamma_alpha_insect   + sigma_sp_insect   * zsp_insect_raw;
  vector[N_sp] sp_disease  = W * gamma_alpha_disease  + sigma_sp_disease  * zsp_disease_raw;
  vector[N_sp] sp_wind     = W * gamma_alpha_wind     + sigma_sp_wind     * zsp_wind_raw;
  vector[N_sp] sp_harvest  = W * gamma_alpha_harvest  + sigma_sp_harvest  * zsp_harvest_raw;
  vector[N_sp] sp_cutting  = W * gamma_alpha_cutting  + sigma_sp_cutting  * zsp_cutting_raw;
  vector[N_sp] sp_siteprep = W * gamma_alpha_siteprep + sigma_sp_siteprep * zsp_siteprep_raw;
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

  // partial-pooling priors: tight sigma so rare-disturbance species shrink
  // toward the trait mean; well-sampled species pull free.
  zsp_plant_raw ~ std_normal();    sigma_sp_plant    ~ normal(0, 0.10);
  zsp_fire_raw ~ std_normal();     sigma_sp_fire     ~ normal(0, 0.10);
  zsp_insect_raw ~ std_normal();   sigma_sp_insect   ~ normal(0, 0.10);
  zsp_disease_raw ~ std_normal();  sigma_sp_disease  ~ normal(0, 0.10);
  zsp_wind_raw ~ std_normal();     sigma_sp_wind     ~ normal(0, 0.10);
  zsp_harvest_raw ~ std_normal();  sigma_sp_harvest  ~ normal(0, 0.10);
  zsp_cutting_raw ~ std_normal();  sigma_sp_cutting  ~ normal(0, 0.10);
  zsp_siteprep_raw ~ std_normal(); sigma_sp_siteprep ~ normal(0, 0.10);

  z_L1_raw    ~ std_normal();
  sigma_L1    ~ normal(0, 0.1);
  sigma_resid ~ normal(0, 0.3);

  vector[N_obs] sp_plant_obs;
  vector[N_obs] sp_fire_obs;
  vector[N_obs] sp_insect_obs;
  vector[N_obs] sp_disease_obs;
  vector[N_obs] sp_wind_obs;
  vector[N_obs] sp_harvest_obs;
  vector[N_obs] sp_cutting_obs;
  vector[N_obs] sp_siteprep_obs;
  for (i in 1:N_obs) {
    int s_i = sp_idx[i];
    sp_plant_obs[i]    = sp_plant[s_i];
    sp_fire_obs[i]     = sp_fire[s_i];
    sp_insect_obs[i]   = sp_insect[s_i];
    sp_disease_obs[i]  = sp_disease[s_i];
    sp_wind_obs[i]     = sp_wind[s_i];
    sp_harvest_obs[i]  = sp_harvest[s_i];
    sp_cutting_obs[i]  = sp_cutting[s_i];
    sp_siteprep_obs[i] = sp_siteprep[s_i];
  }
  vector[N_obs] delta =
      alpha_0
    + (alpha_plant    + sp_plant_obs)    .* is_plantation
    + (alpha_fire     + sp_fire_obs)     .* x_fire
    + (alpha_insect   + sp_insect_obs)   .* x_insect
    + (alpha_disease  + sp_disease_obs)  .* x_disease
    + (alpha_wind     + sp_wind_obs)     .* x_wind
    + (alpha_harvest  + sp_harvest_obs)  .* x_harvest
    + (alpha_cutting  + sp_cutting_obs)  .* x_cutting
    + (alpha_siteprep + sp_siteprep_obs) .* x_siteprep
    + z_L1[L1_idx];

  residual ~ normal(delta, sigma_resid * inv_weight_safe);
}

generated quantities {
  // log_lik for held-out ELPD vs the species-free (traitmed/common) legs
  vector[N_obs] log_lik;
  {
    vector[N_obs] delta;
    for (i in 1:N_obs) {
      int s = sp_idx[i];
      delta[i] =
          alpha_0
        + (alpha_plant    + sp_plant[s])    * is_plantation[i]
        + (alpha_fire     + sp_fire[s])     * x_fire[i]
        + (alpha_insect   + sp_insect[s])   * x_insect[i]
        + (alpha_disease  + sp_disease[s])  * x_disease[i]
        + (alpha_wind     + sp_wind[s])     * x_wind[i]
        + (alpha_harvest  + sp_harvest[s])  * x_harvest[i]
        + (alpha_cutting  + sp_cutting[s])  * x_cutting[i]
        + (alpha_siteprep + sp_siteprep[s]) * x_siteprep[i]
        + z_L1[L1_idx[i]];
      log_lik[i] = normal_lpdf(residual[i] | delta[i], sigma_resid * inv_weight_safe[i]);
    }
  }
  // per-species effective alphas for the deployment table
  vector[N_sp] alpha_fire_sp    = alpha_fire    + sp_fire;
  vector[N_sp] alpha_insect_sp  = alpha_insect  + sp_insect;
  vector[N_sp] alpha_harvest_sp = alpha_harvest + sp_harvest;
  vector[N_sp] alpha_cutting_sp = alpha_cutting + sp_cutting;
}
