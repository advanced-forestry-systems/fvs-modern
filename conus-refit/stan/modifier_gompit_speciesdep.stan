// ============================================================================
// modifier_gompit_speciesdep.stan
//
// SPECIES-DEPENDENT disturbance modifier for the GOMPIT survival base
// (survival_unified_v2_crz, the deployed CONUS mortality form), matching its
// link rather than the logit link of modifier_binary. The base model is
//   log P_surv = -exp(-eta_base) * T_years
// so the modifier enters on the eta (log-hazard) scale:
//   log P_surv_i = -exp( -(eta_base_i + delta_i) ) * T_years_i
//   delta_i = alpha_0
//           + (alpha_k + W[sp_i,]*gamma_k + z_sp_k[sp_i]) * x_k_i   (per term)
//           + z_L1[l1_i]
//   z_sp_k[s] ~ Normal(0, sigma_sp_k)    (trait-informed per-species deviation)
//   alive_i ~ Bernoulli(P_surv_i)
//
// Same b2 species-dependent structure as modifier_binary_speciesdep, only the
// likelihood changes from bernoulli_logit to the gompit survival, so the modifier
// alphas are on the gompit hazard scale and compose correctly with the deployed
// surv_crz eta_base. eta_base is the precomputed posterior-mean linear predictor
// from the surv_crz fit (param_draws + meta), supplied per observation; no
// log_years offset because T_years multiplies the hazard directly.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;
  int<lower=1> N_sp;
  int<lower=0> P_trait;

  array[N_obs] int<lower=0, upper=1> alive;
  vector[N_obs] eta_base;             // surv_crz linear predictor (log-hazard scale)
  vector<lower=0>[N_obs] T_years;     // remeasurement interval

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
}

parameters {
  real alpha_0;
  real alpha_plant; real alpha_fire; real alpha_insect; real alpha_disease;
  real alpha_wind; real alpha_harvest; real alpha_cutting; real alpha_siteprep;

  vector[P_trait] gamma_alpha_plant;   vector[P_trait] gamma_alpha_fire;
  vector[P_trait] gamma_alpha_insect;  vector[P_trait] gamma_alpha_disease;
  vector[P_trait] gamma_alpha_wind;    vector[P_trait] gamma_alpha_harvest;
  vector[P_trait] gamma_alpha_cutting; vector[P_trait] gamma_alpha_siteprep;

  vector[N_sp] zsp_plant_raw;   vector[N_sp] zsp_fire_raw;
  vector[N_sp] zsp_insect_raw;  vector[N_sp] zsp_disease_raw;
  vector[N_sp] zsp_wind_raw;    vector[N_sp] zsp_harvest_raw;
  vector[N_sp] zsp_cutting_raw; vector[N_sp] zsp_siteprep_raw;
  real<lower=1e-4> sigma_sp_plant;  real<lower=1e-4> sigma_sp_fire;
  real<lower=1e-4> sigma_sp_insect; real<lower=1e-4> sigma_sp_disease;
  real<lower=1e-4> sigma_sp_wind;   real<lower=1e-4> sigma_sp_harvest;
  real<lower=1e-4> sigma_sp_cutting;real<lower=1e-4> sigma_sp_siteprep;

  vector[N_L1] z_L1_raw;
  real<lower=1e-4> sigma_L1;
}

transformed parameters {
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
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
  alpha_0 ~ normal(0, 0.1);
  alpha_plant ~ normal(0,0.3); alpha_fire ~ normal(0,0.2); alpha_insect ~ normal(0,0.2);
  alpha_disease ~ normal(0,0.2); alpha_wind ~ normal(0,0.2); alpha_harvest ~ normal(0,0.2);
  alpha_cutting ~ normal(0,0.2); alpha_siteprep ~ normal(0,0.2);

  gamma_alpha_plant ~ normal(0,0.15); gamma_alpha_fire ~ normal(0,0.10);
  gamma_alpha_insect ~ normal(0,0.10); gamma_alpha_disease ~ normal(0,0.10);
  gamma_alpha_wind ~ normal(0,0.10); gamma_alpha_harvest ~ normal(0,0.10);
  gamma_alpha_cutting ~ normal(0,0.10); gamma_alpha_siteprep ~ normal(0,0.10);

  zsp_plant_raw ~ std_normal(); sigma_sp_plant ~ normal(0,0.10);
  zsp_fire_raw ~ std_normal(); sigma_sp_fire ~ normal(0,0.10);
  zsp_insect_raw ~ std_normal(); sigma_sp_insect ~ normal(0,0.10);
  zsp_disease_raw ~ std_normal(); sigma_sp_disease ~ normal(0,0.10);
  zsp_wind_raw ~ std_normal(); sigma_sp_wind ~ normal(0,0.10);
  zsp_harvest_raw ~ std_normal(); sigma_sp_harvest ~ normal(0,0.10);
  zsp_cutting_raw ~ std_normal(); sigma_sp_cutting ~ normal(0,0.10);
  zsp_siteprep_raw ~ std_normal(); sigma_sp_siteprep ~ normal(0,0.10);

  z_L1_raw ~ std_normal(); sigma_L1 ~ normal(0,0.1);

  for (i in 1:N_obs) {
    int s = sp_idx[i];
    real delta =
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
    real eta = fmin(fmax(eta_base[i] + delta, -5.0), 20.0);
    real log_p_surv = -exp(-eta) * T_years[i];
    target += alive[i] == 1 ? log_p_surv : log1m_exp(log_p_surv);
  }
}

generated quantities {
  vector[N_obs] log_lik;
  for (i in 1:N_obs) {
    int s = sp_idx[i];
    real delta =
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
    real eta = fmin(fmax(eta_base[i] + delta, -5.0), 20.0);
    real log_p_surv = -exp(-eta) * T_years[i];
    log_lik[i] = alive[i] == 1 ? log_p_surv : log1m_exp(log_p_surv);
  }
}
