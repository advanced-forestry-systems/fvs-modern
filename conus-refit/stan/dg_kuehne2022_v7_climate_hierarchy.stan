// ============================================================================
// dg_kuehne2022_v7_climate_hierarchy.stan
//
// Architecture variant 7: Aaron's revised strategy (2026-05-15 evening):
//   - Drop SICOND (circular: site index derived from height/age measurements)
//   - Drop ELEV (not trusted)
//   - Use BGI as primary climate-derived site variable (non-circular)
//   - Restore species random effect (z_sp) alongside trait fixed effect
//   - ADD forest type random effect (z_FORTYPCD), 156 levels
//   - Keep ecoregion random effects (z_L1/L2/L3)
//   - Plantation back to modifier (not main effect)
//   - Keep continuous RD modifiers per earlier biometric guidance
//   - Keep additive SDI / SDImax ratio as RD measure but acknowledge it
//     adds little; sdi_complexity stays for structure
//
// Linear predictor:
//   eta = b0 + trait_effect[sp] + z_sp[sp]
//       + z_L1 + z_L2 + z_L3 + z_FORTYPCD
//       + b1 ln(DBH) + b2 DBH + b3 ln_cr_adj
//       + b4 ln_bal_sw_adj + b5 bal_hw
//       + b_site * bgi + b9 bgi^2
//       + b7 BA * rd_additive       (continuous RD interaction)
//       + b8 BAL_SW * rd_additive   (continuous RD interaction)
//       + b11 sdi_complexity
//       + b12 bgi * rd_additive     (continuous site x density modifier)
//
// where:
//   b_site_i = b6 + z_L1_csi[L1_i] + W[sp_i,] * gamma_site
//   bgi: biomass growth index (climate-derived, non-circular)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;             // forest type levels
  int<lower=0> P_trait;

  vector[N_obs] dg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] ln_dbh;
  vector[N_obs] dbh;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] ln_bal_sw_adj;
  vector[N_obs] bal_hw;

  vector[N_obs] bgi;             // biomass growth index (climate-derived)
  vector[N_obs] ba_metric;
  vector[N_obs] bal_sw_metric;
  vector[N_obs] rd_additive;     // sdi_additive1 / SDImax_brms
  vector[N_obs] sdi_complexity;  // sdi_additive1 / SDI_Reineke

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;  // forest type index

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  vector[N_obs] bgi_sq = bgi .* bgi;
  vector[N_obs] ba_x_rdadd     = ba_metric     .* rd_additive;
  vector[N_obs] balsw_x_rdadd  = bal_sw_metric .* rd_additive;
  vector[N_obs] bgi_x_rdadd    = bgi           .* rd_additive;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5;
  real b6;                          // mean linear bgi coef
  real b7;                          // BA x rd_additive
  real b8;                          // BAL_SW x rd_additive
  real b9;                          // quadratic bgi
  real b11;                         // sdi_complexity
  real b12;                         // bgi x rd_additive

  vector[P_trait] gamma;
  vector[P_trait] gamma_site;       // trait-modulated bgi slope

  // Random effects
  vector[N_sp] z_sp_raw;            // species RE (restored)
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;            // forest type RE (NEW)
  vector[N_L1] z_L1_bgi_raw;        // L1-varying bgi slope

  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;           // forest type scale (NEW)
  real<lower=0> sigma_L1_bgi;

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
  vector[N_sp] z_sp     = sigma_sp     * z_sp_raw;
  vector[N_L1] z_L1     = sigma_L1     * z_L1_raw;
  vector[N_L2] z_L2     = sigma_L2     * z_L2_raw;
  vector[N_L3] z_L3     = sigma_L3     * z_L3_raw;
  vector[N_FT] z_FT     = sigma_FT     * z_FT_raw;
  vector[N_L1] z_L1_bgi = sigma_L1_bgi * z_L1_bgi_raw;
}
model {
  b0  ~ normal(-1.0, 2.0);
  b1  ~ normal( 0.3, 0.5);
  b2  ~ normal(-0.02, 0.05);
  b3  ~ normal( 0.5, 0.5);
  b4  ~ normal(-0.05, 0.1);
  b5  ~ normal(-0.005, 0.02);
  b6  ~ normal( 0.05, 0.05);  // bgi coef on its natural scale, positive prior
  b7  ~ normal( 0.0, 0.02);
  b8  ~ normal( 0.0, 0.02);
  b9  ~ normal(-0.01, 0.05);  // quadratic bgi
  b11 ~ normal( 0.0, 0.3);
  b12 ~ normal( 0.0, 0.1);

  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.05);   // tighter prior on trait-bgi interaction

  z_sp_raw     ~ std_normal();
  z_L1_raw     ~ std_normal();
  z_L2_raw     ~ std_normal();
  z_L3_raw     ~ std_normal();
  z_FT_raw     ~ std_normal();
  z_L1_bgi_raw ~ std_normal();

  sigma_sp     ~ normal(0, 0.15);  // tight: traits explain most species variance
  sigma_L1     ~ normal(0, 0.5);
  sigma_L2     ~ normal(0, 0.3);
  sigma_L3     ~ normal(0, 0.3);
  sigma_FT     ~ normal(0, 0.3);   // forest type
  sigma_L1_bgi ~ normal(0, 0.05);

  sigma ~ normal(0, 0.5);

  vector[N_obs] b_site = b6 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];

  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]                                   // forest type RE
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site .* bgi + b9 * bgi_sq
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * bgi_x_rdadd;

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
  vector[N_obs] b_site_gq = b6 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];
  vector[N_obs] eta_gq =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + z_FT[FT_idx]
    + b1 * ln_dbh + b2 * dbh + b3 * ln_cr_adj
    + b4 * ln_bal_sw_adj + b5 * bal_hw
    + b_site_gq .* bgi + b9 * bgi_sq
    + b7  * ba_x_rdadd
    + b8  * balsw_x_rdadd
    + b11 * sdi_complexity
    + b12 * bgi_x_rdadd;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
