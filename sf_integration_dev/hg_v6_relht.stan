// ============================================================================
// hg_organon_speciesfree_v5_bgi.stan
//
// Height growth (ORGANON form) species-free B1, v5 architecture.
// Modernizes v4_full by:
//   - Replacing ln_csi with BGI piecewise basis (matching DG_Kuehne v8):
//     3-piece basis (linear + 2 hinges at empirical knots)
//   - Adding forest type RE (z_FT)
//   - Keeping L1/L2/L3 ecoregion + trait-modulated site slope
//
// Linear predictor:
//   ln(HTG_a) = a0 + trait_effect[sp]
//             + z_L1 + z_L2 + z_L3 + z_FT
//             + a1 ln(DBH) + a2 ln(HT) + a3 ln((CR+0.2)/1.2)
//             + b_site * bgi + a9a bgi_b2 + a9b bgi_b3
//             + a5 BAL + a6 BA + a7 SLOPE + a8 cos(ASP)
//             + a10 bgi * bal_log    <-- continuous site x competition
// where b_site_i = a4 + z_L1_bgi[L1_i] + W[sp_i,] * gamma_site
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector[N_obs] hg_obs_a;
  vector[N_obs] sqrt_years;

  vector[N_obs] ln_dbh;
  vector<lower=0>[N_obs] ht_obs_m;          // raw height in meters
  vector<lower=0>[N_sp] max_ht_sp_m;         // species max height (raw, meters)
  vector[N_obs] ln_cr_adj;
  vector[N_obs] bal_log;
  vector[N_obs] bgi;
  vector[N_obs] ba_metric;
  vector[N_obs] slope_pct;
  vector[N_obs] cos_aspect;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;

  real bgi_knot1;
  real bgi_knot2;
}
transformed data {
  vector[N_obs] bgi_b1 = bgi;
  vector[N_obs] bgi_b2;
  vector[N_obs] bgi_b3;
  for (i in 1:N_obs) {
    bgi_b2[i] = fmax(bgi[i] - bgi_knot1, 0.0);
    bgi_b3[i] = fmax(bgi[i] - bgi_knot2, 0.0);
  }
  vector[N_obs] bgi_x_bal = bgi .* bal_log;
  vector[N_obs] relht;
  vector[N_obs] log_one_minus_relht;
  for (i in 1:N_obs) {
    real r = ht_obs_m[i] / fmax(max_ht_sp_m[sp_idx[i]], 0.1);
    relht[i] = r;
    log_one_minus_relht[i] = log(fmax(1.0 - r, 0.01));
  }
}
parameters {
  real a0; real a1; real a2; real a3;
  real a4;
  real a5; real a6; real a7; real a8;
  real a9a; real a9b;
  real a10;
  vector[P_trait] gamma;
  vector[P_trait] gamma_site;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;
  vector[N_sp] z_sp_raw;            // NEW species RE
  vector[N_L1] z_L1_bgi_raw;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;
  real<lower=0> sigma_sp;           // NEW species RE scale
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
  vector[N_L1] z_L1     = sigma_L1     * z_L1_raw;
  vector[N_L2] z_L2     = sigma_L2     * z_L2_raw;
  vector[N_L3] z_L3     = sigma_L3     * z_L3_raw;
  vector[N_FT] z_FT     = sigma_FT     * z_FT_raw;
  vector[N_sp] z_sp     = sigma_sp     * z_sp_raw;
  vector[N_L1] z_L1_bgi = sigma_L1_bgi * z_L1_bgi_raw;
}
model {
  a0  ~ normal(-1.0, 2.0);
  a1  ~ normal( 0.0, 0.5);
  a2  ~ normal( 0.5, 1.0);     // asymptotic shutdown strength
  a3  ~ normal( 0.5, 0.5);
  a4  ~ normal( 0.05, 0.1);     // bgi mean slope, weak positive
  a5  ~ normal(-0.005, 0.02);   // BAL
  a6  ~ normal(-0.005, 0.02);   // BA
  a7  ~ normal( 0.0, 0.01);     // slope
  a8  ~ normal( 0.0, 0.1);      // aspect
  a9a ~ normal( 0.0, 0.1);
  a9b ~ normal( 0.0, 0.1);
  a10 ~ normal( 0.0, 0.05);

  gamma      ~ normal(0, 0.5);
  gamma_site ~ normal(0, 0.05);

  z_L1_raw     ~ std_normal();
  z_L2_raw     ~ std_normal();
  z_L3_raw     ~ std_normal();
  z_FT_raw     ~ std_normal();
  z_sp_raw     ~ std_normal();
  z_L1_bgi_raw ~ std_normal();

  sigma_L1     ~ normal(0, 0.5);
  sigma_L2     ~ normal(0, 0.3);
  sigma_L3     ~ normal(0, 0.3);
  sigma_FT     ~ normal(0, 0.3);
  sigma_sp     ~ normal(0, 0.5);
  sigma_L1_bgi ~ normal(0, 0.05);

  sigma ~ normal(0, 0.5);

  vector[N_obs] b_site = a4 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];

  vector[N_obs] eta =
      a0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + a1 * ln_dbh + a2 * log_one_minus_relht + a3 * ln_cr_adj
    + b_site .* bgi_b1
    + a9a * bgi_b2 + a9b * bgi_b3
    + a5 * bal_log + a6 * ba_metric
    + a7 * slope_pct + a8 * cos_aspect
    + a10 * bgi_x_bal;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -30.0), 20.0);
  for (i in 1:N_obs) {
    if (hg_obs_a[i] > 0.001) {
      target += lognormal_lpdf(hg_obs_a[i] | eta_safe[i],
                fmin(fmax(sigma / sqrt_years[i], 1e-4), 50.0));
    }
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] b_site_gq = a4 + z_L1_bgi[L1_idx] + species_site_slope[sp_idx];
  vector[N_obs] eta_gq =
      a0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + a1 * ln_dbh + a2 * log_one_minus_relht + a3 * ln_cr_adj
    + b_site_gq .* bgi_b1
    + a9a * bgi_b2 + a9b * bgi_b3
    + a5 * bal_log + a6 * ba_metric
    + a7 * slope_pct + a8 * cos_aspect
    + a10 * bgi_x_bal;
  vector[N_obs] mu_a;
  for (i in 1:N_obs) mu_a[i] = exp(fmin(eta_gq[i], 20.0));
  for (i in 1:N_obs) log_lik[i] = normal_lpdf(hg_obs_a[i] | mu_a[i], sigma / sqrt_years[i]);
}
