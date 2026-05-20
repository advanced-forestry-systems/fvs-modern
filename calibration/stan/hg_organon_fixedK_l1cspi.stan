// ============================================================================
// hg_organon_fixedK_l1cspi.stan
//
// HG preflight that adds an EPA L1 ecoregion-varying slope on
// centered ln(CSPI). Tests whether the U-shape in the quadratic
// preflight (a4 negative, a4b positive) is driven by different CSPI
// responses across ecoregions.
//
// The plant x CSPI preflight ruled out plantation-driven slope
// variation. The next plausible driver is regional heterogeneity:
// PNW high-CSPI Douglas-fir likely has a different log-HG response
// than SE high-CSPI loblolly pine, and the L1 random intercept
// already in production absorbs the means but not the slopes.
//
// Linear predictor (log annual height growth, log m/yr):
//   eta = a0 + trait_effect[sp] + z_sp[sp]
//       + z_L1[l1] + z_L2[l2] + z_L3[l3]
//       + a1  ln(HT + K1h) + a2  HT^K2h + a3  ln((CR+0.2)/1.2)
//       + a4   * ln_cspi_c                          # global slope
//       + a4b  * ln_cspi_c^2                        # global curvature
//       + a_l1_cspi[l1] * ln_cspi_c                 # L1-varying slope
//       + a5 CCFL/ln(HT+K4h) + a6 sqrt(BA)
//       + a7 BA*RD + a8 CCFL*RD
//
// Decision rule:
//   * If sigma_a_l1_cspi > 0.10 (substantial L1 variation in slope),
//     the U-shape is regional. Production HG should have an
//     L1-varying CSPI slope.
//   * If sigma_a_l1_cspi < 0.05 and individual a_l1_cspi entries
//     are tight around zero, regional slope variation is small;
//     U-shape has yet another driver (stand age, species clusters,
//     or a CSPI definition artifact).
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

  real ln_cspi_mean = mean(ln_cspi_shift);
  vector[N_obs] ln_cspi_c = ln_cspi_shift - ln_cspi_mean;
  vector[N_obs] ln_cspi_c2 = ln_cspi_c .* ln_cspi_c;
}

parameters {
  real a0;
  real a1;
  real a2;
  real a3;
  real a4;
  real a4b;
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

  // L1-varying CSPI slope (non-centered)
  vector[N_L1] a_l1_cspi_raw;
  real<lower=0> sigma_a_l1_cspi;

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
  vector[N_L1] a_l1_cspi = sigma_a_l1_cspi * a_l1_cspi_raw;
}

model {
  a0  ~ normal(-2.5, 1.5);
  a1  ~ normal( 0.5, 0.5);
  a2  ~ normal(-0.03, 0.03);
  a3  ~ normal( 0.8, 0.5);
  a4  ~ normal( 0.4, 0.5);
  a4b ~ normal( 0.0, 0.3);
  a5  ~ normal(-0.005, 0.02);
  a6  ~ normal(-0.03, 0.05);
  a7  ~ normal( 0.0, 0.02);
  a8  ~ normal( 0.0, 0.02);

  gamma     ~ normal(0, 0.5);
  z_sp_raw  ~ std_normal();
  z_L1_raw  ~ std_normal();
  z_L2_raw  ~ std_normal();
  z_L3_raw  ~ std_normal();

  a_l1_cspi_raw    ~ std_normal();
  sigma_a_l1_cspi  ~ normal(0, 0.2);

  sigma_sp ~ normal(0, 0.15);
  sigma_L1 ~ normal(0, 0.3);
  sigma_L2 ~ normal(0, 0.2);
  sigma_L3 ~ normal(0, 0.2);

  sigma    ~ normal(0, 0.5);

  // Pre-extract the L1-varying CSPI slope to a per-obs vector
  // to avoid a CmdStan compiler optimization trigger when an
  // indexed vector is multiplied elementwise with another vector.
  vector[N_obs] a_l1_slope_obs;
  for (i in 1:N_obs)
    a_l1_slope_obs[i] = a_l1_cspi[L1_idx[i]] * ln_cspi_c[i];

  vector[N_obs] eta =
      a0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a1  * ln_ht_k1
    + a2  * ht_k2
    + a3  * ln_cr_adj
    + a4  * ln_cspi_c
    + a4b * ln_cspi_c2
    + a_l1_slope_obs
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
