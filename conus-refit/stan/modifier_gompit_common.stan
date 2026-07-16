// modifier_gompit_common.stan
// COMMON (species-free, trait-free) gompit disturbance modifier: each
// disturbance effect is a single global alpha shared across all species, plus an
// ecoregion (L1) random effect. No trait projection, no per-species z_sp. The
// most parsimonious leg of the 3-way comparison (common < traitmed < speciesdep).
// W is accepted in data for an identical interface but is unused here.
data {
  int<lower=1> N_obs;
  int<lower=1> N_L1;
  int<lower=1> N_sp;
  int<lower=0> P_trait;

  array[N_obs] int<lower=0, upper=1> alive;
  vector[N_obs] eta_base;
  vector<lower=0>[N_obs] T_years;

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

  vector[N_L1] z_L1_raw;
  real<lower=1e-4> sigma_L1;
}

transformed parameters {
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
}

model {
  alpha_0 ~ normal(0, 0.1);
  alpha_plant ~ normal(0,0.3); alpha_fire ~ normal(0,0.2); alpha_insect ~ normal(0,0.2);
  alpha_disease ~ normal(0,0.2); alpha_wind ~ normal(0,0.2); alpha_harvest ~ normal(0,0.2);
  alpha_cutting ~ normal(0,0.2); alpha_siteprep ~ normal(0,0.2);

  z_L1_raw ~ std_normal(); sigma_L1 ~ normal(0,0.1);

  for (i in 1:N_obs) {
    real delta =
        alpha_0
      + alpha_plant    * is_plantation[i]
      + alpha_fire     * x_fire[i]
      + alpha_insect   * x_insect[i]
      + alpha_disease  * x_disease[i]
      + alpha_wind     * x_wind[i]
      + alpha_harvest  * x_harvest[i]
      + alpha_cutting  * x_cutting[i]
      + alpha_siteprep * x_siteprep[i]
      + z_L1[L1_idx[i]];
    real eta = fmin(fmax(eta_base[i] + delta, -5.0), 20.0);
    real log_p_surv = -exp(-eta) * T_years[i];
    target += alive[i] == 1 ? log_p_surv : log1m_exp(log_p_surv);
  }
}

generated quantities {
  vector[N_obs] log_lik;
  for (i in 1:N_obs) {
    real delta =
        alpha_0
      + alpha_plant    * is_plantation[i]
      + alpha_fire     * x_fire[i]
      + alpha_insect   * x_insect[i]
      + alpha_disease  * x_disease[i]
      + alpha_wind     * x_wind[i]
      + alpha_harvest  * x_harvest[i]
      + alpha_cutting  * x_cutting[i]
      + alpha_siteprep * x_siteprep[i]
      + z_L1[L1_idx[i]];
    real eta = fmin(fmax(eta_base[i] + delta, -5.0), 20.0);
    real log_p_surv = -exp(-eta) * T_years[i];
    log_lik[i] = alive[i] == 1 ? log_p_surv : log1m_exp(log_p_surv);
  }
}
