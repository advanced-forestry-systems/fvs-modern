// ============================================================================
// dg_weiskittel2016.stan
//
// FVS-CONUS diameter growth (DG) candidate form: Weiskittel 2016 / Kuehne
// 2020 SS_REAL "single-species realized" unified equation, extended with
// nested EPA L1 / L2 / L3 ecoregion random intercepts and BAPA x RD,
// BAL x RD interactions.
//
// Units (enforced upstream by 30_build_conus_dataset.R Step 5b):
//   DBH in cm, HT in m, BA / BAL in m^2/ha, CSPI unitless z composite,
//   RD additive SDI / SDImax_brms (plot-level).
//
// Reparameterization (2026-04-16):
//   The original form used g3 * CR (raw, 0-1) which allowed the linear CR
//   coefficient to push exp(eta) to extreme values when g3 is large,
//   causing divergences in the preflight. The fix is to use the same
//   log-adjusted CR transform as the organon and kuehne forms:
//     ln_cr_adj = log((CR + 0.2) / 1.2)
//   This maps (0, 1) -> (-inf, 0) with a small shift, stabilizing sampling.
//
// Linear predictor (annual diameter growth, cm/yr):
//   ln(dDBH_a) = g0 + z_sp[sp] + z_L1[l1] + z_L2[l2] + z_L3[l3]
//              + g1 ln(DBH)
//              + g2 DBH
//              + g3 ln((CR + 0.2) / 1.2)
//              + g4 BAL_SW
//              + g5 BAL_HW
//              + g6 RD
//              + g7 ln(CSPI_shift)
//              + g8 BA  x RD
//              + g9 BAL x RD
//
// Observed annual growth dDBH_obs_a = (DBH2 - DBH1) / YEARS.
// Observation variance scales as sigma / sqrt(YEARS) per Cao (2000).
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;         // species-level trait columns (0 = disable prior)

  vector[N_obs] dg_obs_a;       // observed annual DBH growth (cm/yr)
  vector[N_obs] sqrt_years;     // sqrt(YEARS) for variance scaling

  vector[N_obs] ln_dbh;
  vector[N_obs] dbh;
  vector[N_obs] ln_cr_adj;          // log((CR + 0.2) / 1.2), computed R-side
  vector[N_obs] bal_sw;
  vector[N_obs] bal_hw;
  vector[N_obs] rd;
  vector[N_obs] ln_cspi_shift;  // log(CSPI + shift) computed R-side
  vector[N_obs] ba_x_rd;
  vector[N_obs] bal_x_rd;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;  // species x trait design (centered, scaled)
}
parameters {
  real g0;
  real g1;
  real g2;
  real g3;
  real g4;
  real g5;
  real g6;
  real g7;
  real g8;
  real g9;

  vector[P_trait] gamma;        // trait loadings mapping W -> species RE mean

  // Non-centered RE
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
  // Priors: weak but informative, metric scale
  g0 ~ normal(-1.5, 2.0);
  g1 ~ normal( 0.4, 0.5);
  g2 ~ normal(-0.02, 0.05);
  g3 ~ normal( 0.5, 0.5);   // tighter after log-CR reparameterization
  g4 ~ normal(-0.02, 0.05);
  g5 ~ normal(-0.02, 0.05);
  g6 ~ normal(-1.0, 1.0);
  g7 ~ normal( 0.3, 0.5);
  g8 ~ normal( 0.0, 0.02);
  g9 ~ normal( 0.0, 0.02);

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

  vector[N_obs] eta =
      g0
    + z_sp[sp_idx]
    + z_L1[L1_idx]
    + z_L2[L2_idx]
    + z_L3[L3_idx]
    + g1 * ln_dbh
    + g2 * dbh
    + g3 * ln_cr_adj
    + g4 * bal_sw
    + g5 * bal_hw
    + g6 * rd
    + g7 * ln_cspi_shift
    + g8 * ba_x_rd
    + g9 * bal_x_rd;

  // Normal likelihood on annual growth with Cao (2000) variance scaling.
  // sigma_eff = sigma / sqrt(years); vectorized via elementwise division.
  dg_obs_a ~ normal(exp(eta), sigma ./ sqrt_years);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_a = exp(
      g0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + g1 * ln_dbh + g2 * dbh + g3 * ln_cr_adj
    + g4 * bal_sw + g5 * bal_hw + g6 * rd
    + g7 * ln_cspi_shift + g8 * ba_x_rd + g9 * bal_x_rd
  );
  for (i in 1:N_obs) {
    log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
  }
}
