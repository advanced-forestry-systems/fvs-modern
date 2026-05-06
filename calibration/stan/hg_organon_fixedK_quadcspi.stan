// ============================================================================
// hg_organon_fixedK_quadcspi.stan
//
// Variant of hg_organon_fixedK.stan that adds a centered quadratic
// site-productivity term to test the curvilinear hypothesis. The
// production HG fit returned a4 = -0.23 on linear ln(CSPI_shift),
// opposite the expected sign. Two diagnoses:
//   1. The relationship between log growth and log CSPI is curvilinear
//      (peaks at intermediate productivity, declines at extremes)
//   2. Hidden interactions confound the linear term
//
// This preflight tests #1 by adding `a4b * (ln_cspi_shift - mean)^2`.
// Centering before squaring decorrelates the linear and quadratic
// columns so a4 (linear) and a4b (curvature) are individually
// identifiable.
//
// Linear predictor (log annual height growth, log m/yr):
//   ln(dHT_a) = a0 + trait_effect[sp] + z_sp[sp]
//             + z_L1[l1] + z_L2[l2] + z_L3[l3]
//             + a1 ln(HT + K1h) + a2 HT^K2h + a3 ln((CR + 0.2)/1.2)
//             + a4 (ln_cspi_shift - mean)
//             + a4b (ln_cspi_shift - mean)^2
//             + a5 CCFL / ln(HT + K4h)
//             + a6 sqrt(BA) + a7 BA x RD + a8 CCFL x RD
//
// Note: a4 here is the slope at the CSPI mean, NOT the same as a4 in
// the production model where ln_cspi_shift was uncentered. Compare
// curvature (a4b) and signs to interpret.
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

  // Centered ln_CSPI to decorrelate linear and quadratic terms
  real ln_cspi_mean = mean(ln_cspi_shift);
  vector[N_obs] ln_cspi_c = ln_cspi_shift - ln_cspi_mean;
  vector[N_obs] ln_cspi_c2 = ln_cspi_c .* ln_cspi_c;
}

parameters {
  real a0;
  real a1;
  real a2;
  real a3;
  real a4;       // slope at the CSPI mean
  real a4b;      // curvature (negative => concave-down, peak at center)
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
  a0  ~ normal(-2.5, 1.5);
  a1  ~ normal( 0.5, 0.5);
  a2  ~ normal(-0.03, 0.03);
  a3  ~ normal( 0.8, 0.5);
  a4  ~ normal( 0.4, 0.5);    // slope at mean CSPI; expect positive
  a4b ~ normal( 0.0, 0.3);    // curvature; symmetric prior
  a5  ~ normal(-0.005, 0.02);
  a6  ~ normal(-0.03, 0.05);
  a7  ~ normal( 0.0, 0.02);
  a8  ~ normal( 0.0, 0.02);

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
    + a1  * ln_ht_k1
    + a2  * ht_k2
    + a3  * ln_cr_adj
    + a4  * ln_cspi_c
    + a4b * ln_cspi_c2
    + a5  * comp
    + a6  * sqrt_ba
    + a7  * ba_x_rd
    + a8  * ccfl_x_rd;

  vector[N_obs] eta_clamped;
  for (i in 1:N_obs) eta_clamped[i] = fmax(fmin(eta[i], 30.0), -30.0);
  log_hg_obs_a ~ normal(eta_clamped, sigma ./ sqrt_years);
}

generated quantities {
  real ln_cspi_mean_out = mean(ln_cspi_shift);
}
