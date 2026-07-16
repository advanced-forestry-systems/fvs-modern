// ============================================================================
// dg_kuehne2022_speciesfree_v6_data_driven.stan
//
// Architecture variant 6: data-driven covariate set selected from the
// 2026-05-15 covariate exploration analysis. Builds on v5 by:
//
//   - Site variable: SICOND (FIA Site Index, ft @ base age 50). Univariate
//     R^2 = 0.044 for DG with positive sign; strongest among tested site
//     metrics. Replaces ln_csi (climate_si) which had R^2 = 0.007.
//
//   - Competition: CCFL1 (Crown Competition Factor for larger trees). Best
//     univariate predictor of DG (R^2 = 0.107). Used as primary competition
//     measure alongside the existing ln_bal_sw_adj and bal_hw.
//
//   - Plantation main effect: is_plantation (0/1 indicator from STDORGCD).
//     Conditional R^2 added = 0.107 for DG. Currently absorbed into modifier
//     alphas; promoted to main effect.
//
//   - Elevation: ELEV (m). Conditional R^2 added = 0.070 for DG. Topographic
//     position matters beyond climate.
//
//   - Stand structural complexity: sdi_complexity (sdi_additive1 / SDI1).
//     Continuous, captures uneven-aged structure. Modest signal (R^2 added
//     ~0.01 for DG) but interpretable.
//
//   - Continuous RD interactions (Aaron's biometric guidance):
//       b * ln_sicond * rd_additive  (continuous site-density modifier)
//
// All RD-based modifiers are continuous variable multipliers, not bounded
// indicators.
//
// Linear predictor:
//   eta = b0 + trait_effect[sp] + z_L1 + z_L2 + z_L3
//       + b1 ln(DBH) + b2 DBH + b3 ln_cr_adj
//       + b4 ln_bal_sw_adj + b5 bal_hw
//       + b_site * ln_sicond + b9 ln_sicond^2
//       + b7 ccfl1
//       + b8 is_plantation
//       + b10 ln_elev
//       + b11 sdi_complexity
//       + b12 ln_sicond * rd_additive   <-- continuous site x density modifier
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=0> P_trait;

  vector[N_obs] dg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] ln_dbh;
  vector[N_obs] dbh;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] ln_bal_sw_adj;
  vector[N_obs] bal_hw;

  vector[N_obs] ln_sicond;              // log of FIA Site Index (SICOND)
  vector[N_obs] ccfl1;                  // crown competition factor (larger trees)
  vector[N_obs] is_plantation;          // 0 / 1
  vector[N_obs] ln_elev;                // log of elevation (m)
  vector[N_obs] sdi_complexity;         // sdi_additive1 / SDI_Reineke
  vector[N_obs] rd_additive;            // sdi_additive1 / SDImax_brms

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  vector[N_obs] ln_sicond_sq = ln_sicond .* ln_sicond;
  vector[N_obs] sicond_x_rdadd = ln_sicond .* rd_additive;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5;
  real b6;                              // mean linear site coef on ln_sicond
  real b7;                              // ccfl1
  real b8;                              // is_plantation
  real b9;                              // quadratic site
  real b10;                             // ln_elev
  real b11;                             // sdi_complexity
  real b12;                             // ln_sicond x rd_additive

  vector[P_trait] gamma;
  vector[P_trait] gamma_site;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_L1] z_L1_csi_raw;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_L1_csi;

  real<lower=0> sigma;
}
transformed parameters {
  vector[N_sp] trait_effect;
  vector[N_sp] species_site_slope;
  if (P_trait > 0) {
    trait_effect       = W * gamma;
    species_site_slope = W * gamma_site;
  } else {
    trait_effect       = rep_vector(0.0, N_sp);
    species_site_slope = rep_vector(0.0, N_sp);
  }
  vector[N_L1] z_L1     = sigma_L1     * z_L1_raw;
  vector[N_L2] z_L2     = sigma_L2     * z_L2_raw;
  vector[N_L3] z_L3     = sigma_L3     * z_L3_raw;
  vector[N_L1] z_L1_csi = sigma_L1_csi * z_L1_csi_raw;
}
model {
  b0  ~ normal(-1.0, 2.0);
  b1  ~ normal( 0.3, 0.5);
  b2  ~ normal(-0.02, 0.05);
  b3  ~ normal( 0.5, 0.5);
  b4  ~ normal(-0.05, 0.1);
  b5  ~ normal(-0.005, 0.02);
  b6  ~ normal( 0.4, 0.3);        // positive prior on ln(SICOND) coef
  b7  ~ normal(-0.001, 0.005);    // CCFL1 negative effect (competition)
  b8  ~ normal( 0.3, 0.3);        // is_plantation: positive prior (plantations grow faster)
  b9  ~ normal(-0.1, 0.3);        // quadratic site
  b10 ~ normal(-0.1, 0.3);        // ELEV: weakly negative (cooler, less productive)
  b11 ~ normal( 0.0, 0.3);        // sdi_complexity
  b12 ~ normal( 0.0, 0.3);        // site x rd interaction

  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.3);
  z_L1_raw     ~ std_normal();
  z_L2_raw     ~ std_normal();
  z_L3_raw     ~ std_normal();
  z_L1_csi_raw ~ std_normal();
  sigma_L1     ~ normal(0, 0.5);
  sigma_L2     ~ normal(0, 0.3);
  sigma_L3     ~ normal(0, 0.3);
  sigma_L1_csi ~ normal(0, 0.3);
  sigma        ~ normal(0, 0.5);

  vector[N_obs] b_site = b6 + z_L1_csi[L1_idx] + species_site_slope[sp_idx];

  vector[N_obs] eta =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site .* ln_sicond + b9 * ln_sicond_sq
    + b7  * ccfl1
    + b8  * is_plantation
    + b10 * ln_elev
    + b11 * sdi_complexity
    + b12 * sicond_x_rdadd;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -30.0), 20.0);
  for (i in 1:N_obs) {
    if (dg_obs_a[i] > 0.001) {
      target += lognormal_lpdf(dg_obs_a[i] | eta_safe[i],
                fmin(fmax(sigma / sqrt_years[i], 1e-4), 50.0));
    }
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] b_site_gq = b6 + z_L1_csi[L1_idx] + species_site_slope[sp_idx];
  vector[N_obs] eta_gq =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site_gq .* ln_sicond + b9 * ln_sicond_sq
    + b7  * ccfl1
    + b8  * is_plantation
    + b10 * ln_elev
    + b11 * sdi_complexity
    + b12 * sicond_x_rdadd;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
