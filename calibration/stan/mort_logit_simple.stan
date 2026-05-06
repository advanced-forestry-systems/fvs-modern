// ============================================================================
// mort_logit_simple.stan
//
// FVS-CONUS individual tree mortality with logit link and simplified
// hierarchy: species + EPA L1 random intercepts only (L2/L3 dropped).
//
// Motivation: mort_gompit.stan (cloglog, full L1/L2/L3 hierarchy) failed
// to converge even with tightened priors (25% divergences, rhat = Inf).
// The deep nesting (L1->L2->L3 + species) created near-singularity in
// the random effects structure given the low per-cell mortality rates.
// Reducing to L1 + species dramatically cuts the number of latent
// parameters and eliminates the nested confounding.
//
// The logit link simplifies the likelihood (Stan's bernoulli_logit is
// vectorized and numerically stable) and is standard for binary outcomes.
// Annualization uses a log(YEARS) offset on the logit scale, which is an
// approximation but performs well empirically when intervals are short
// (1-20 years) and mortality rates are low (<10% annual).
//
// Probability of dying during the interval:
//   eta[i]    = x[i] . beta + z_sp[sp_i] + z_L1[L1_i] + log(YEARS[i])
//   p_die[i]  = inv_logit(eta[i])
//   alive[i] ~ bernoulli(1 - p_die[i])
//
// Predictors:
//   dbh, dbh^2, bal_over_ba, cr, sqrt_ba_rd, ln_cspi_shift
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;

  int<lower=0> P_trait;                      // species-level trait columns (0 = disable)

  array[N_obs] int<lower=0, upper=1> alive;  // 1 = survived, 0 = died
  vector[N_obs] log_years;                   // exposure offset

  vector[N_obs] dbh;
  vector[N_obs] dbh2;
  vector[N_obs] bal_over_ba;
  vector[N_obs] cr;
  vector[N_obs] sqrt_ba_rd;
  vector[N_obs] ln_cspi_shift;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;  // species x trait design
}
parameters {
  real m0;
  real m1;  // DBH
  real m2;  // DBH^2
  real m3;  // BAL / BA
  real m4;  // CR
  real m5;  // sqrt(BA * RD)
  real m6;  // ln(CSPI_shift)

  vector[P_trait] gamma;                     // trait loadings

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
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
}
model {
  // Fixed effect priors (same sign expectations as gompit version,
  // but on the logit scale rather than log-hazard scale)
  m0 ~ normal(-4.0, 2.0);        // baseline logit(p_die) — most trees survive
  m1 ~ normal(-0.05, 0.05);      // smaller trees die more
  m2 ~ normal( 0.0, 0.001);      // DBH^2 concavity
  m3 ~ normal( 1.0, 1.0);        // high competition boosts mortality
  m4 ~ normal(-2.0, 1.0);        // higher CR lowers mortality
  m5 ~ normal( 0.1, 0.1);        // density effect
  m6 ~ normal(-0.2, 0.3);        // better site lowers mortality

  gamma ~ normal(0, 0.5);

  // Non-centered random effects
  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();

  // Moderate RE priors — with only 2 RE groups (not 4), the sampler
  // has much less risk of confounding. Keep these moderately tight.
  sigma_sp ~ normal(0, 0.5);
  sigma_L1 ~ normal(0, 0.3);

  // Logit link: bernoulli_logit is vectorized and numerically stable.
  // Note: alive=1 means survived, so we model p(survive) = 1 - p(die).
  // bernoulli_logit(alive | -eta) is equivalent to bernoulli(alive | 1-invlogit(eta)).
  vector[N_obs] eta =
      m0
    + z_sp[sp_idx] + z_L1[L1_idx]
    + m1 * dbh
    + m2 * dbh2
    + m3 * bal_over_ba
    + m4 * cr
    + m5 * sqrt_ba_rd
    + m6 * ln_cspi_shift
    + log_years;                                  // exposure offset

  // Model survival: alive ~ bernoulli_logit(-eta)
  // This is equivalent to: p_die = inv_logit(eta), alive ~ bernoulli(1-p_die)
  alive ~ bernoulli_logit(-eta);
}
generated quantities {
  vector[N_obs] p_die;
  vector[N_obs] log_lik;
  {
    vector[N_obs] eta =
        m0
      + z_sp[sp_idx] + z_L1[L1_idx]
      + m1 * dbh + m2 * dbh2 + m3 * bal_over_ba + m4 * cr
      + m5 * sqrt_ba_rd + m6 * ln_cspi_shift
      + log_years;
    for (i in 1:N_obs) {
      p_die[i] = inv_logit(eta[i]);
      // log_lik for LOO: log P(alive[i] | params)
      if (alive[i] == 1) {
        log_lik[i] = log1m(p_die[i]);    // log(1 - p_die)
      } else {
        log_lik[i] = log(p_die[i]);      // log(p_die)
      }
    }
  }
}
