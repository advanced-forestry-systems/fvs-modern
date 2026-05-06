// ============================================================================
// hg_organon_fixedK_plantcspi.stan
//
// Variant of hg_organon_fixedK_quadcspi.stan that adds:
//   1. Plantation main effect (a_plant * is_plantation)
//   2. Plantation x centered-CSPI interaction
//      (a_plant_cspi * is_plantation * ln_cspi_c)
//
// The quadratic preflight returned a U-shaped CSPI response (a4 negative
// at the mean, a4b positive). That curvature is well identified but
// ecologically suspicious. Most likely interpretation: the fit is
// being pulled into a U because the plantation subset (concentrated
// at high CSPI in the SE pine belt) has a different slope on log-CSPI
// than the natural-stand subset.
//
// This preflight tests that hypothesis by letting the plantation
// subset have its own intercept AND its own linear-CSPI slope.
//
// Linear predictor (log annual height growth, log m/yr):
//   ln(dHT_a) = a0 + trait_effect[sp] + z_sp[sp]
//             + z_L1[l1] + z_L2[l2] + z_L3[l3]
//             + a1  ln(HT + K1h) + a2  HT^K2h + a3  ln((CR+0.2)/1.2)
//             + a4  (ln_cspi_shift - mean)
//             + a4b (ln_cspi_shift - mean)^2
//             + a_plant      is_plantation
//             + a_plant_cspi is_plantation * (ln_cspi_shift - mean)
//             + a5  CCFL/ln(HT+K4h) + a6 sqrt(BA)
//             + a7  BA*RD + a8 CCFL*RD
//
// Decision rule for the next round:
//   * If a_plant_cspi is large and positive, plantations have a
//     stronger positive CSPI slope than natural stands. Production HG
//     should split slope by STDORGCD or fit a stratified model.
//   * If a_plant_cspi is near zero, the plantation effect is on the
//     intercept only and the U-shape has another driver (likely L1).
//   * If a4b is no longer needed once the interaction is in (posterior
//     for a4b shrinks toward zero), the curvature was a confound from
//     the plantation slope difference.
// ============================================================================

data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;

  vector[N_obs] log_hg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] ht;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] ln_cspi_shift;
  vector[N_obs] ccfl;
  vector[N_obs] sqrt_ba;
  vector[N_obs] ba_x_rd;
  vector[N_obs] ccfl_x_rd;

  vector<lower=0, upper=1>[N_obs] is_plantation;

  real<lower=0.01> K1h;
  real<lower=0.01> K2h;
  real<lower=0.01> K4h;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}

transformed data {
  vector[N_obs] ln_ht_k1 = log(ht + K1h);
  vector[N_obs] ht_k2;
  vector[N_obs] comp;
  for (i in 1:N_obs) ht_k2[i] = pow(ht[i], K2h);
  comp = ccfl ./ log(ht + K4h);

  real ln_cspi_mean = mean(ln_cspi_shift);
  vector[N_obs] ln_cspi_c = ln_cspi_shift - ln_cspi_mean;
  vector[N_obs] ln_cspi_c2 = ln_cspi_c .* ln_cspi_c;
  vector[N_obs] plant_x_cspi = is_plantation .* ln_cspi_c;
}

parameters {
  real a0;
  real a1;
  real a2;
  real a3;
  real a4;
  real a4b;
  real a_plant;
  real a_plant_cspi;
  real a5;
  real a6;
  real a7;
  real a8;

  vector[P_trait] gamma;
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
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
}

model {
  a0           ~ normal(-2.5, 1.5);
  a1           ~ normal( 0.5, 0.5);
  a2           ~ normal(-0.03, 0.03);
  a3           ~ normal( 0.8, 0.5);
  a4           ~ normal( 0.4, 0.5);
  a4b          ~ normal( 0.0, 0.3);
  a_plant      ~ normal( 0.0, 0.3);
  a_plant_cspi ~ normal( 0.0, 0.3);
  a5           ~ normal(-0.005, 0.02);
  a6           ~ normal(-0.03, 0.05);
  a7           ~ normal( 0.0, 0.02);
  a8           ~ normal( 0.0, 0.02);

  gamma     ~ normal(0, 0.5);
  z_sp_raw  ~ std_normal();
  z_L1_raw  ~ std_normal();
  z_L2_raw  ~ std_normal();
  z_L3_raw  ~ std_normal();

  sigma_sp ~ normal(0, 0.15);
  sigma_L1 ~ normal(0, 0.3);
  sigma_L2 ~ normal(0, 0.2);
  sigma_L3 ~ normal(0, 0.2);

  sigma    ~ normal(0, 0.5);

  vector[N_obs] eta =
      a0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a1           * ln_ht_k1
    + a2           * ht_k2
    + a3           * ln_cr_adj
    + a4           * ln_cspi_c
    + a4b          * ln_cspi_c2
    + a_plant      * is_plantation
    + a_plant_cspi * plant_x_cspi
    + a5           * comp
    + a6           * sqrt_ba
    + a7           * ba_x_rd
    + a8           * ccfl_x_rd;

  vector[N_obs] eta_clamped;
  for (i in 1:N_obs) eta_clamped[i] = fmax(fmin(eta[i], 30.0), -30.0);
  log_hg_obs_a ~ normal(eta_clamped, sigma ./ sqrt_years);
}

generated quantities {
  real ln_cspi_mean_out = mean(ln_cspi_shift);
  real plantation_share = mean(is_plantation);
}
