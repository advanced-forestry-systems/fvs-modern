// ORGANON form Height to Crown Base Model for CONUS wide FVS
// Based on Hanus, Hann, and Marshall (2000) HCB_SWO equation
// Extended with species x ecodivision hierarchy and relative density
//
// ALL UNITS METRIC: HT in m, DBH in cm, BA in m2/ha, site productivity in m
//
// Model form (logistic proportion of total height):
//   HCB = HT / (1 + exp(b0[s,d] + b1*HT + b2*CCFL + b3*ln(BA)
//                       + b4*(DBH/HT) + b5*SITE_PROD + b6*RD))
//
// This predicts HCB as a static function of current tree and stand attributes.
// Crown ratio is derived as: CR = 1 - HCB/HT
// For projection, crown recession: HCBG = HCB(t2) - HCB(t1), constrained >= 0
//
// Site productivity: ClimateSI, BGI (biomass growth index), or climate PCA
// (NO traditional FIA SI)
//
// Hierarchy:
//   b0[s,d] = mu_b0 + sigma_sp * z_sp[s] + sigma_eco * z_eco[d]

data {
  int<lower=1> N;                      // Number of tree observations
  int<lower=1> N_species;
  int<lower=1> N_ecodiv;

  // Response: observed height to crown base (m)
  vector<lower=0>[N] hcb_obs;

  // Tree and stand attributes (all metric)
  vector<lower=0>[N] ht;              // Total height (m)
  vector[N] ccfl;                      // Crown competition factor in larger trees
  vector[N] ln_ba;                     // ln(basal area in m2/ha)
  vector[N] dbh_ht_ratio;             // DBH(cm) / HT(m) (form ratio, mixed units OK)
  vector[N] site_prod;                 // ClimateSI (m), BGI, or climate PC score
  vector[N] rd;                        // Relative density (Curtis: BA/sqrt(QMD))

  // Group indices
  array[N] int<lower=1, upper=N_species> species_id;
  array[N] int<lower=1, upper=N_ecodiv> ecodiv_id;
}

transformed data {
  // Cap posterior-predictive and log_lik arrays at 2000 obs so they fit in
  // memory and the Stan compiler accepts a data-time size expression.
  int<lower=0> M = min(N, 2000);
}

parameters {
  // Hierarchical intercept
  real mu_b0;
  real<lower=0.01> sigma_sp;
  real<lower=0.01> sigma_eco;
  vector[N_species] z_sp;
  vector[N_ecodiv] z_eco;

  // Fixed slopes
  real b1;                             // HT effect (negative: taller = higher HCB)
  real b2;                             // CCFL effect (negative: more competition = higher HCB)
  real b3;                             // ln(BA) effect (negative: denser = higher HCB)
  real b4;                             // DBH/HT effect (positive: stocky trees = lower HCB)
  real b5;                             // Site productivity effect
  real b6;                             // RD effect (negative: denser = higher HCB)

  // Observation error (m)
  real<lower=0.01> sigma;
}

transformed parameters {
  vector[N_species] b0_sp = sigma_sp * z_sp;
  vector[N_ecodiv] b0_eco = sigma_eco * z_eco;
}

model {
  // ---- Priors ----
  // Informed by ORGANON HCB_SWO coefficient ranges, adjusted for metric units
  mu_b0 ~ normal(2.0, 2.0);           // ~1.8 to 3.8 in SWO
  sigma_sp ~ exponential(1);
  sigma_eco ~ exponential(2);
  z_sp ~ std_normal();
  z_eco ~ std_normal();

  // HT prior adjusted for meters (original SWO was per foot, ~0.008/ft)
  // In meters: b1 ~ -0.008 * 3.281 ~ -0.026
  b1 ~ normal(-0.026, 0.015);         // Negative HT effect (per meter)
  b2 ~ normal(-0.003, 0.002);         // Negative CCFL effect (unitless)
  b3 ~ normal(-0.5, 0.3);             // Negative ln(BA) effect
  b4 ~ normal(3.5, 2.0);              // Positive form ratio effect
  b5 ~ normal(0.005, 0.005);          // Site productivity effect
  b6 ~ normal(-0.3, 0.3);             // Negative RD effect (denser = higher HCB)

  sigma ~ exponential(1);

  // ---- Likelihood ----
  {
    vector[N] hcb_pred;
    for (n in 1:N) {
      real logit_cr = mu_b0 + b0_sp[species_id[n]] + b0_eco[ecodiv_id[n]]
                    + b1 * ht[n]
                    + b2 * ccfl[n]
                    + b3 * ln_ba[n]
                    + b4 * dbh_ht_ratio[n]
                    + b5 * site_prod[n]
                    + b6 * rd[n];

      // HCB = HT / (1 + exp(logit_cr))
      hcb_pred[n] = ht[n] / (1.0 + exp(logit_cr));
    }

    hcb_obs ~ normal(hcb_pred, sigma);
  }
}

generated quantities {
  vector[M] hcb_rep;
  vector[M] log_lik;

  for (n in 1:M) {
    real logit_cr_n = mu_b0 + b0_sp[species_id[n]] + b0_eco[ecodiv_id[n]]
                    + b1 * ht[n]
                    + b2 * ccfl[n]
                    + b3 * ln_ba[n]
                    + b4 * dbh_ht_ratio[n]
                    + b5 * site_prod[n]
                    + b6 * rd[n];
    real hcb_pred_n = ht[n] / (1.0 + exp(logit_cr_n));

    hcb_rep[n] = normal_rng(hcb_pred_n, sigma);
    log_lik[n] = normal_lpdf(hcb_obs[n] | hcb_pred_n, sigma);
  }
}
