// ============================================================================
// dg_organon_gj.stan
//
// FVS-CONUS diameter growth candidate form: Greg Johnson simplified
// ORGANON. Adapted from the draft per-species DG fits in the
// fvs_remodeling repository (github.com/gregjohnsonbiometrics/fvs_remodeling).
//
// Key simplifications relative to the full ORGANON / SWO form:
//   * No site productivity term (CSPI / BGI / ClimateSI / BYI all omitted).
//   * No BA*RD or BAL*RD interactions.
//   * K1 and K4 fixed at 1 and 1 rather than estimated.
//   * K2 fixed at 1 (no DBH^K2 curvature; linear DBH effect only).
//
// Serves as a null / baseline candidate in the five-form comparison:
// if the four CSPI-bearing forms fail to beat this one on hold-out MAB,
// the site productivity contribution does not survive.
//
// Linear predictor (annual diameter growth, cm/yr):
//   ln(dDBH_a) = s0 + z_sp[sp] + z_L1[l1] + z_L2[l2] + z_L3[l3]
//              + s1 ln(DBH + 1)
//              + s2 DBH
//              + s3 ln((CR + 0.2)/1.2)
//              + s4 BAL / ln(DBH + 1)
//              + s5 sqrt(BA)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;            // species-level trait columns (0 = disable prior)

  vector[N_obs] dg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] dbh;
  vector[N_obs] ln_cr_adj;
  vector[N_obs] bal;
  vector[N_obs] sqrt_ba;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;  // species x trait design (centered, scaled)
}
parameters {
  real s0;
  real s1;
  real s2;
  real s3;
  real s4;
  real s5;

  vector[P_trait] gamma;           // trait loadings mapping W -> species RE mean

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
  s0 ~ normal(-1.5, 2.0);
  s1 ~ normal( 0.4, 0.5);
  s2 ~ normal(-0.02, 0.05);
  s3 ~ normal( 0.8, 0.5);
  s4 ~ normal(-0.005, 0.02);
  s5 ~ normal(-0.03, 0.05);

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

  vector[N_obs] ln_dbh_p1 = log(dbh + 1.0);
  vector[N_obs] comp      = bal ./ ln_dbh_p1;

  vector[N_obs] eta =
      s0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + s1 * ln_dbh_p1
    + s2 * dbh
    + s3 * ln_cr_adj
    + s4 * comp
    + s5 * sqrt_ba;

  dg_obs_a ~ normal(exp(eta), sigma ./ sqrt_years);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_a;
  {
    vector[N_obs] ln_dbh_p1 = log(dbh + 1.0);
    vector[N_obs] comp      = bal ./ ln_dbh_p1;
    mu_a = exp(
        s0
      + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + s1 * ln_dbh_p1 + s2 * dbh + s3 * ln_cr_adj
      + s4 * comp + s5 * sqrt_ba
    );
  }
  for (i in 1:N_obs) {
    log_lik[i] = normal_lpdf(dg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
  }
}
