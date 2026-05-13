// ============================================================================
// hg_organon.stan
//
// FVS-CONUS height growth, ORGANON-analogue form with estimated K1h, K2h
// and the same nested EPA L1 / L2 / L3 + species random intercept
// hierarchy used by DG. This replaces the legacy scripts/organon_hg_conus.stan
// (crossed species + ecodiv only) as production HG.
//
// Linear predictor (annual height growth, m/yr):
//   ln(dHT_a) = a0 + z_sp[sp] + z_L1[l1] + z_L2[l2] + z_L3[l3]
//             + a1 ln(HT + K1h)
//             + a2 HT^K2h
//             + a3 ln((CR + 0.2)/1.2)
//             + a4 ln(CSPI_shift)
//             + a5 CCFL / ln(HT + K4h)
//             + a6 sqrt(BA)
//             + a7 BA  x RD
//             + a8 CCFL x RD
//
// Annualization: preflight uses Cao variance scaling (sigma / sqrt(years));
// production wraps the prediction in the iterative Cao/Weiskittel loop on
// the R side or via a brms custom family.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;               // species-level trait columns (0 = disable prior)

  vector[N_obs] hg_obs_a;             // annual HT increment, m/yr
  vector[N_obs] sqrt_years;

  vector[N_obs] ht;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] ln_cspi_shift;
  vector[N_obs] ccfl;
  vector[N_obs] sqrt_ba;
  vector[N_obs] ba_x_rd;
  vector[N_obs] ccfl_x_rd;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;  // species x trait design (centered, scaled)
}
parameters {
  real a0;
  real a1;
  real a2;
  real a3;
  real a4;
  real a5;
  real a6;
  real a7;
  real a8;

  real<lower=0.1>  K1h;
  real<lower=0.05> K2h;
  real<lower=0.5>  K4h;

  vector[P_trait] gamma;              // trait loadings mapping W -> species RE mean

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

  real<lower=0> sigma;
}
transformed parameters {
  vector[N_sp] mu_sp;
  if (P_trait > 0) {
    mu_sp = W * gamma;
  } else {
    mu_sp = rep_vector(0.0, N_sp);
  }
  vector[N_sp] z_sp = mu_sp + sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
}
model {
  a0 ~ normal(-2.5, 1.5);
  a1 ~ normal( 0.5, 0.5);
  a2 ~ normal(-0.03, 0.03);
  a3 ~ normal( 0.8, 0.5);
  a4 ~ normal( 0.4, 0.5);
  a5 ~ normal(-0.005, 0.02);
  a6 ~ normal(-0.03, 0.05);
  a7 ~ normal( 0.0, 0.02);
  a8 ~ normal( 0.0, 0.02);

  K1h ~ normal(1.5, 0.75);
  K2h ~ normal(0.8, 0.3);
  K4h ~ normal(2.7, 1.0);

  gamma ~ normal(0, 0.5);
  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.5);
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);

  sigma ~ normal(0, 0.5);

  vector[N_obs] ln_ht_k1 = log(ht + K1h);
  vector[N_obs] ht_k2;
  for (i in 1:N_obs) ht_k2[i] = pow(ht[i], K2h);
  vector[N_obs] comp = ccfl ./ log(ht + K4h);

  vector[N_obs] eta =
      a0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a1 * ln_ht_k1
    + a2 * ht_k2
    + a3 * ln_cr_adj
    + a4 * ln_cspi_shift
    + a5 * comp
    + a6 * sqrt_ba
    + a7 * ba_x_rd
    + a8 * ccfl_x_rd;

  // Cap eta to prevent exp() overflow during warmup (exp(20) ~ 5e8)
  vector[N_obs] mu_safe;
  for (i in 1:N_obs) mu_safe[i] = exp(fmin(eta[i], 20.0));

  hg_obs_a ~ normal(mu_safe, sigma ./ sqrt_years);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_a;
  {
    vector[N_obs] ln_ht_k1 = log(ht + K1h);
    vector[N_obs] ht_k2;
    for (i in 1:N_obs) ht_k2[i] = pow(ht[i], K2h);
    vector[N_obs] comp = ccfl ./ log(ht + K4h);
    {
      vector[N_obs] eta_gq =
          a0
        + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
        + a1 * ln_ht_k1 + a2 * ht_k2 + a3 * ln_cr_adj
        + a4 * ln_cspi_shift + a5 * comp + a6 * sqrt_ba
        + a7 * ba_x_rd + a8 * ccfl_x_rd;
      for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
    }
  }
  for (i in 1:N_obs) {
    log_lik[i] = normal_lpdf(hg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
  }
}
