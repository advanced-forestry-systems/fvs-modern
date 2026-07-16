// crown_ratio_change_fixed_z.stan
// Standardized variant of crown_ratio_change_fixed.stan. The unstandardized
// model failed to mix (rhat 2.4-3.8, ESS ~4 on ALL params incl sigma/intercept)
// because raw predictors span wildly different scales (dbh~10s, dbh^2~1000s,
// ba/bal large) and dbh/dbh^2 are collinear. This z-scores every predictor in
// transformed data and builds the quadratic from the CENTERED dbh so linear
// and quadratic terms decorrelate. Identical likelihood, well-conditioned.
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=0> P_trait;
  vector[N_obs] delta_CR_a;
  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;      // passed by driver; ignored (recomputed from z_dbh)
  vector[N_obs] ba_metric;
  vector[N_obs] bal_metric;
  vector[N_obs] cr_init;
  vector[N_obs] ln_csi;
  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  vector[N_obs] z_dbh; vector[N_obs] z_dbh2; vector[N_obs] z_ba;
  vector[N_obs] z_bal; vector[N_obs] z_cr;   vector[N_obs] z_csi;
  {
    real s;
    z_dbh = (dbh - mean(dbh)) / sd(dbh);
    vector[N_obs] dq = z_dbh .* z_dbh;
    z_dbh2 = (dq - mean(dq)) / sd(dq);
    z_ba  = (ba_metric  - mean(ba_metric))  / sd(ba_metric);
    z_bal = (bal_metric - mean(bal_metric)) / sd(bal_metric);
    z_cr  = (cr_init    - mean(cr_init))    / sd(cr_init);
    s = sd(ln_csi); if (s <= 0) s = 1;
    z_csi = (ln_csi - mean(ln_csi)) / s;
  }
}
parameters {
  real b0; real b1; real b2; real b3; real b4; real b5; real b6;
  vector[P_trait] gamma;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N_sp] trait_effect;
  if (P_trait > 0) trait_effect = W * gamma; else trait_effect = rep_vector(0.0, N_sp);
}
model {
  b0 ~ normal(0.0, 0.25);
  b1 ~ normal(0.0, 0.1); b2 ~ normal(0.0, 0.1);
  b3 ~ normal(0.0, 0.1); b4 ~ normal(0.0, 0.1);
  b5 ~ normal(0.0, 0.1); b6 ~ normal(0.0, 0.1);
  gamma ~ normal(0, 0.5);
  sigma ~ normal(0, 0.5);
  vector[N_obs] mu = b0 + trait_effect[sp_idx]
    + b1*z_dbh + b2*z_dbh2 + b3*z_ba + b4*z_bal + b5*z_cr + b6*z_csi;
  delta_CR_a ~ normal(mu, sigma);
}
generated quantities {
  vector[N_obs] log_lik;
  {
    vector[N_obs] mu_pred = b0 + trait_effect[sp_idx]
      + b1*z_dbh + b2*z_dbh2 + b3*z_ba + b4*z_bal + b5*z_cr + b6*z_csi;
    for (i in 1:N_obs) log_lik[i] = normal_lpdf(delta_CR_a[i] | mu_pred[i], sigma);
  }
}
