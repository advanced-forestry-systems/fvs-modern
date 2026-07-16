// ORGANON-form Diameter Growth Model for CONUS-wide FVS -- THREADED v2
// Identical model to organon_dg_conus.stan (Hann SWO form, species x ecodiv
// hierarchy, K1/K2 estimated, direct-PAI likelihood) but the per-observation
// likelihood is wrapped in reduce_sum() so each chain runs multi-threaded.
// Parameter names are byte-for-byte identical to the single-thread model so
// the downstream extractor (61_extract_conus_summaries.R) and validate
// functions in 31_fit_dg_organon_v2.R work unchanged.

functions {
  // Partial sum over a slice of observations [start:end] (indices into the
  // sliced arrays). reduce_sum passes the first sliced arg (dg_obs_s) plus
  // start/end; every other tree-level array is sliced in lock-step.
  real dg_partial_lpdf(array[] real dg_obs_s, int start, int end,
                       real mu_b0,
                       vector b0_sp, vector b0_eco,
                       real b1, real b2, real b3, real b4, real b5,
                       real b6, real b7, real b8, real b9, real b10,
                       real K1, real K2, real sigma,
                       array[] real dbh_s, array[] real ln_cr_adj_s,
                       array[] real ln_site_prod_s, array[] real bal_ratio_s,
                       array[] real ln_bal_s, array[] real sqrt_ba_s,
                       array[] real clim1_s, array[] real clim2_s,
                       array[] real rd_s, array[] real years_s,
                       array[] int species_id_s, array[] int ecodiv_id_s) {
    real lp = 0;
    int M = end - start + 1;
    for (i in 1:M) {
      real ln_dg = mu_b0 + b0_sp[species_id_s[i]] + b0_eco[ecodiv_id_s[i]]
                   + b1 * log(dbh_s[i] + K1)
                   + b2 * pow(dbh_s[i], K2)
                   + b3 * ln_cr_adj_s[i]
                   + b4 * ln_site_prod_s[i]
                   + b5 * bal_ratio_s[i]
                   + b6 * sqrt_ba_s[i]
                   + b7 * clim1_s[i]
                   + b8 * clim2_s[i]
                   + b9 * rd_s[i]
                   + b10 * rd_s[i] * ln_bal_s[i];
      real ln_dg_safe = fmin(fmax(ln_dg, -30.0), 20.0);
      real dg_pred = exp(ln_dg_safe) * years_s[i];
      if (dg_pred > 0.001) {
        lp += normal_lpdf(dg_obs_s[i] | dg_pred,
                          sigma * sqrt(fmax(dg_pred, 0.01)));
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> N_species;
  int<lower=1> N_ecodiv;
  int<lower=1> grainsize;            // reduce_sum chunk size (tuning knob)

  vector[N] dg_obs;
  vector[N] dbh;
  vector[N] ln_cr_adj;
  vector[N] ln_site_prod;
  vector[N] bal_ratio;
  vector[N] ln_bal;
  vector[N] sqrt_ba;
  vector[N] clim1;
  vector[N] clim2;
  vector[N] rd;
  vector<lower=0>[N] years;
  array[N] int<lower=1, upper=N_species> species_id;
  array[N] int<lower=1, upper=N_ecodiv> ecodiv_id;
}

transformed data {
  // reduce_sum slices the *first* argument, which must be an array. Convert
  // the tree-level vectors to real arrays once here so they can be sliced.
  array[N] real dg_obs_a       = to_array_1d(dg_obs);
  array[N] real dbh_a          = to_array_1d(dbh);
  array[N] real ln_cr_adj_a    = to_array_1d(ln_cr_adj);
  array[N] real ln_site_prod_a = to_array_1d(ln_site_prod);
  array[N] real bal_ratio_a    = to_array_1d(bal_ratio);
  array[N] real ln_bal_a       = to_array_1d(ln_bal);
  array[N] real sqrt_ba_a      = to_array_1d(sqrt_ba);
  array[N] real clim1_a        = to_array_1d(clim1);
  array[N] real clim2_a        = to_array_1d(clim2);
  array[N] real rd_a           = to_array_1d(rd);
  array[N] real years_a        = to_array_1d(years);
}

parameters {
  real mu_b0;
  real<lower=0.01> sigma_sp;
  real<lower=0.01> sigma_eco;
  vector[N_species] z_sp;
  vector[N_ecodiv] z_eco;
  real<lower=0.1> K1;
  real<lower=0.01> K2;
  real b1; real b2; real b3; real b4; real b5;
  real b6; real b7; real b8; real b9; real b10;
  real<lower=0.01> sigma;
}

transformed parameters {
  vector[N_species] b0_sp = sigma_sp * z_sp;
  vector[N_ecodiv]  b0_eco = sigma_eco * z_eco;
}

model {
  // ---- Priors (identical to single-thread model) ----
  mu_b0 ~ normal(-5.0, 2.0);
  sigma_sp ~ exponential(1);
  sigma_eco ~ exponential(2);
  z_sp ~ std_normal();
  z_eco ~ std_normal();
  K1 ~ normal(6.0, 2.0);
  K2 ~ normal(1.0, 0.5);
  b1 ~ normal(0.8, 0.5);
  b2 ~ normal(-0.04, 0.03);
  b3 ~ normal(1.0, 0.5);
  b4 ~ normal(0.8, 0.5);
  b5 ~ normal(-0.008, 0.005);
  b6 ~ normal(-0.04, 0.03);
  b7 ~ normal(0, 1);
  b8 ~ normal(0, 1);
  b9 ~ normal(-0.5, 0.5);
  b10 ~ normal(0, 0.3);
  sigma ~ exponential(1);

  // ---- Threaded likelihood ----
  target += reduce_sum(dg_partial_lpdf, dg_obs_a, grainsize,
                       mu_b0, b0_sp, b0_eco,
                       b1, b2, b3, b4, b5, b6, b7, b8, b9, b10,
                       K1, K2, sigma,
                       dbh_a, ln_cr_adj_a, ln_site_prod_a, bal_ratio_a,
                       ln_bal_a, sqrt_ba_a, clim1_a, clim2_a, rd_a, years_a,
                       species_id, ecodiv_id);
}

generated quantities {
}
