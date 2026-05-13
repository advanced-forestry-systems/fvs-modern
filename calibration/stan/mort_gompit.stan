// ============================================================================
// mort_gompit.stan
//
// FVS-CONUS individual tree mortality with complementary log-log (gompit)
// link and a log(YEARS) exposure offset, harmonized to the nested EPA
// L1 / L2 / L3 + species random intercept hierarchy shared with DG and HG.
//
// This is the Gompit.e.m1 formulation that won the mortality example
// comparison against logit alternatives: the cloglog link derives from a
// continuous-time proportional hazards model, so exp(x * beta) is
// interpreted as an annual hazard and any measurement interval converts
// correctly through the exposure offset log(YEARS).
//
// Probability of dying during the interval:
//   eta[i]    = x[i] . beta + REs[i] + log(YEARS[i])
//   p_die[i]  = 1 - exp(-exp(eta[i]))
//   alive[i] ~ bernoulli(1 - p_die[i])
//
// Predictors:
//   dbh, dbh^2, bal_over_ba, cr, sqrt_ba_rd, ln_cspi_shift
//   (softwood/hardwood BAL split not used for mortality)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;                      // species-level trait columns (0 = disable prior)

  array[N_obs] int<lower=0, upper=1> alive;  // 1 = survived, 0 = died (AGENTCD < 80)
  vector[N_obs] log_years;                   // exposure offset

  vector[N_obs] dbh;
  vector[N_obs] dbh2;
  vector[N_obs] bal_over_ba;
  vector[N_obs] cr;
  vector[N_obs] sqrt_ba_rd;
  vector[N_obs] ln_cspi_shift;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;  // species x trait design (centered, scaled)
}
parameters {
  real m0;
  real m1;  // DBH
  real m2;  // DBH^2
  real m3;  // BAL / BA
  real m4;  // CR
  real m5;  // sqrt(BA * RD)
  real m6;  // ln(CSPI_shift)

  vector[P_trait] gamma;                     // trait loadings mapping W -> species RE mean

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
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
  m0 ~ normal(-5.0, 2.0);        // baseline annual log-hazard
  m1 ~ normal(-0.05, 0.05);      // smaller trees die more
  m2 ~ normal( 0.0, 0.001);      // DBH^2 adds concavity
  m3 ~ normal( 1.0, 1.0);        // BAL / BA ratio boosts hazard
  m4 ~ normal(-1.5, 1.0);        // higher CR lowers hazard
  m5 ~ normal( 0.1, 0.1);
  m6 ~ normal(-0.2, 0.3);        // better site lowers hazard

  gamma ~ normal(0, 0.5);
  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  // Tighter RE priors (v2): prevent random effects from absorbing all signal.
  // Original values (1.0 / 0.5 / 0.3 / 0.3) led to sigma_sp > 3, rhat = Inf.
  sigma_sp ~ normal(0, 0.5);         // was 1.0
  sigma_L1 ~ normal(0, 0.3);         // was 0.5
  sigma_L2 ~ normal(0, 0.15);        // was 0.3
  sigma_L3 ~ normal(0, 0.15);        // was 0.3

  vector[N_obs] eta =
      m0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + m1 * dbh
    + m2 * dbh2
    + m3 * bal_over_ba
    + m4 * cr
    + m5 * sqrt_ba_rd
    + m6 * ln_cspi_shift
    + log_years;                                  // exposure offset

  // cloglog likelihood: bernoulli(p_die) with p_die = 1 - exp(-exp(eta))
  // Stan's bernoulli_logit works on logit; we use the direct cloglog form.
  for (i in 1:N_obs) {
    real log_surv = -exp(eta[i]);                 // log(1 - p_die)
    if (alive[i] == 1) {
      target += log_surv;
    } else {
      target += log1m_exp(log_surv);              // log(p_die)
    }
  }
}
generated quantities {
  vector[N_obs] p_die;
  vector[N_obs] log_lik;
  {
    vector[N_obs] eta =
        m0
      + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + m1 * dbh + m2 * dbh2 + m3 * bal_over_ba + m4 * cr
      + m5 * sqrt_ba_rd + m6 * ln_cspi_shift
      + log_years;
    for (i in 1:N_obs) {
      real log_surv = -exp(eta[i]);
      p_die[i] = 1 - exp(log_surv);
      if (alive[i] == 1) {
        log_lik[i] = log_surv;
      } else {
        log_lik[i] = log1m_exp(log_surv);
      }
    }
  }
}
