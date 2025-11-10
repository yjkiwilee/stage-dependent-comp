// STAN script for testing the use of a Bayesian model to determine stage-dependence
// This model is to be used for continuous growth.
// Written by Young Jun Lee
// Feb 2024

data {
  // Data attributes
  int n; // Total number of focal observations
  int n_n; // Number of neighbours to calculate crowding effect from
  int n_plots; // Number of plots
  int n_years; // Number of years
  
  // Plant lengths
  vector[n] l_f;
  vector[n_n] l_n;
  
  // Overlap matrix
  matrix[n, n_n] overlap;
  
  // Overlap from non-focal neighbours
  vector[n] overlap_other;
  
  // Response vital rate
  vector[n] vr;
  
  // Random effects; must be 1 ~ n_plots/n_years
  int plot[n];
  int year[n];
  
  // Parameters to control priors
  
  // Intercept without random effect of year & plot
  real rand_intercept;
  // Random effects of year & plot
  vector[n_plots] plot_effs;
  vector[n_years] year_effs;
  
  // Vital rate sd
  real vr_sd;
}

parameters {
  // Intercepts
  real a_o;
  real a;
  
  // Slopes
  real b_f;
  real b_n;
  real b_i;
  real b_o;
  real b_l;
  
  // Global sigma
  real<lower=0> s;
  
  // Random effect intercepts
  vector[n_plots] a_plot;
  vector[n_years] a_year;
}

model {
  // Weakly informative priors
  // Priors are calculated to match the expected statistical properties of the
  // vital rate.
  
  a_o ~ normal(0, 1);
  
  // Use intercept after random effect
  a ~ normal(rand_intercept, 1);
  
  b_f ~ normal(0, 1);
  b_n ~ normal(0, 1);
  b_i ~ normal(0, 1);
  b_l ~ normal(0, 1);
  b_o ~ normal(0, 1);
  
  s ~ normal(vr_sd, 1);
  
  // Random effect priors
  a_plot ~ normal(plot_effs, rep_vector(1, n_plots));
  a_year ~ normal(year_effs, rep_vector(1, n_years));
  
  // Calculate contribution from crowding
  
  row_vector[n_n] a_o_vec;
  row_vector[n_n] b_f_vec;
  row_vector[n_n] b_n_vec;
  row_vector[n_n] b_i_vec;
  
  a_o_vec = rep_row_vector(a_o, n_n);
  b_f_vec = rep_row_vector(b_f, n_n);
  b_n_vec = rep_row_vector(b_n, n_n);
  b_i_vec = rep_row_vector(b_i, n_n);
  
  vector[n] ones_n;
  ones_n = rep_vector(1, n);
  vector[n_n] ones_n_n;
  ones_n_n = rep_vector(1, n_n);
  
  vector[n] crwd;
  
  crwd = (overlap .* (
    (ones_n * a_o_vec) +
    (l_f * b_f_vec) +
    (ones_n * (l_n' .* b_n_vec)) +
    (l_f * (l_n' .* b_i_vec))
  )) * ones_n_n;
  
  // Calculate vectors with random effects
  
  vector[n] r_eff;
  for(i in 1:n) {
    r_eff[i] = a_plot[plot[i]] + a_year[year[i]];
  }
  
  // Calculate final vital rates
  
  vector[n] mu;
  mu = crwd + b_o * overlap_other + b_l * l_f + r_eff + a;
  vr ~ normal(mu, s);
}



