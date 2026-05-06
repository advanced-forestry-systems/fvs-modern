// ============================================================================
// hg_organon_fixedK.stan  (v4 — traits as fixed effects)
//
// FVS-CONUS height growth, ORGANON-analogue form with K1h, K2h, K4h FIXED
// as data inputs rather than estimated.
//
// v2: Log-normal likelihood, tighter RE priors.
// v3: eta cap in generated quantities.
// v4: Trait architecture refactored. W * gamma now enters the linear
//     predictor as a FIXED EFFECT (trait_effect), separate from the species
//     random intercept. z_sp is purely residual (zero-mean) with a tighter
//     sigma_sp prior to encourage shrinkage. This gives interpretable
//     gamma coefficients and a model that generalizes across species via
//     functional traits, with the species RE limited to capturing residual
//     variation not explained by traits.
//
// Linear predictor (log annual height growth, log m/yr):
//   ln(dHT_a) = a0 + trait_effect[sp] + z_sp[sp] + z_L1[l1] + z_L2[l2] + z_L3[l3]
//             + a1 ln(HT + K1h) + a2 HT^K2h + a3 ln((CR + 0.2)/1.2)
//             + a4 ln(CSPI_shift) + a5 CCFL / ln(HT + K4h)
//             + a6 sqrt(BA) + a7 BA x RD + a8 CCFL x RD
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
  // Traits as fixed effects (separate from species RE)
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
  // Species RE: purely residual, zero-mean
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
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

  gamma ~ normal(0, 0.5);
  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();

  sigma_sp ~ normal(0, 0.15);        // tight: traits should explain most species variation
  sigma_L1 ~ normal(0, 0.3);
  sigma_L2 ~ normal(0, 0.2);
  sigma_L3 ~ normal(0, 0.2);

  sigma ~ normal(0, 0.5);

  vector[N_obs] eta =
      a0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a1 * ln_ht_k1
    + a2 * ht_k2
    + a3 * ln_cr_adj
    + a4 * ln_cspi_shift
    + a5 * comp
    + a6 * sqrt_ba
    + a7 * ba_x_rd
    + a8 * ccfl_x_rd;

  // Clamp eta to keep the likelihood location finite during
  // warmup. Growth on a log(m/yr) scale should live in roughly
  // [-8, 3], so bounding at [-30, 30] is a no-op in the well-
  // identified region but prevents init rejection cascades.
  vector[N_obs] eta_clamped;
  for (i in 1:N_obs) eta_clamped[i] = fmax(fmin(eta[i], 30.0), -30.0);
  log_hg_obs_a ~ normal(eta_clamped, sigma ./ sqrt_years);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_a;
  {
    vector[N_obs] eta =
        a0
      + trait_effect[sp_idx] + z_sp[sp_idx]
      + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + a1 * ln_ht_k1 + a2 * ht_k2 + a3 * ln_cr_adj
      + a4 * ln_cspi_shift + a5 * comp + a6 * sqrt_ba
      + a7 * ba_x_rd + a8 * ccfl_x_rd;
    for (i in 1:N_obs) mu_a[i] = exp(fmin(eta[i], 20.0));
    for (i in 1:N_obs) {
      log_lik[i] = normal_lpdf(log_hg_obs_a[i] | eta[i], sigma / sqrt_years[i]);
    }
  }
}
