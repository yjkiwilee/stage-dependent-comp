// Script for testing the use of a Bayesian model to determine stage-dependence
// Written by Young Jun Lee
// Feb 2024

data {
  // Data attributes
  int n_f; // Total number of focal plants
  int n_c; // Number of conspecifics to calculate crowding effect from
  int n_h; // Number of specific competitors to calculate crowding effect from
  int n_o; // Number of other competitors to calculate crowding effect from
  int n_plots; // Number of plots
  int n_years; // Number of years
  
  // Plant lengths
  vector[n_f] l_f;
  vector[n_c] l_c;
  vector[n_h] l_h;
  vector[n_o] l_o;
  
  // Distance matrices
  matrix[n_f, n_c] dist_c;
  matrix[n_f, n_h] dist_h;
  matrix[n_f, n_o] dist_o;
  
  // Response vital rate
  vector[n_f] f_vr;
  
  // Random effects; must be 1 ~ n_plots/n_years
  int plot[n_f];
  int year[n_f];
  
  // Parameters for prior specification
  real e_vr; // Mean vital rate
  
}

parameters {
  // See daily log 08032025 for notation
  
  // Shape parameter for decay kernel
  real d_t;
  
  // Slopes
  real b_c;
  real b_co;
  real b_ci;
  real b_h;
  real b_ho;
  real b_hi;
  real b_o;
  real b_oo;
  real b_oi;
  real b_l;
  
  // Intercepts
  real a_c;
  real a_h;
  real a_o;
  real err;
  
  // Random effect intercepts
  vector[n_plots] a_plot;
  vector[n_years] a_year;
}

model {
  // Priors
  
  real crwd_sd;
  crwd_sd = 0.5;
  
  d_t ~ normal(30, 5);
  
  a_c ~ normal(0, 0.5);
  a_h ~ normal(0, 0.5);
  a_o ~ normal(0, 0.5);
  
  b_c ~ normal(0, )
  
  // Transform distance matrices with Gaussian decay kernel
  
  matrix[n_f, n_c] dist_c_t;
  matrix[n_f, n_h] dist_h_t;
  matrix[n_f, n_o] dist_o_t;
  
  dist_c_st = exp((log(0.5) / d_t^2) .* (dist_c_st .* dist_c_st));
  dist_h_st = exp((log(0.5) / d_t^2) .* (dist_h_st .* dist_h_st));
  dist_o_st = exp((log(0.5) / d_t^2) .* (dist_o_st .* dist_o_st));
  
  // Calculate contribution from crowding
  
  vector[n_f] crwd_c;
  vector[n_f] crwd_h;
  vector[n_f] crwd_o;
  
  crwd_c = (dist_c_st .* (
    a_c +
    (b_c * (l_f * rep_row_vector(1, n_c))) +
    (b_co * (rep_vector(1, n_f) * l_c)) +
    (b_ci * (l_f * l_c'))
  )) * rep_vector(1, n_c);
  
  crwd_h = (dist_h_st .* (
    a_h +
    (b_h * (l_f * rep_row_vector(1, n_h))) +
    (b_ho * (rep_vector(1, n_f) * l_h)) +
    (b_hi * (l_f * l_h'))
  )) * rep_vector(1, n_h);
  
  crwd_o = (dist_o_st .* (
    a_o +
    (b_o * (l_f * rep_row_vector(1, n_o))) +
    (b_oo * (rep_vector(1, n_f) * l_o)) +
    (b_oi * (l_f * l_o'))
  )) * rep_vector(1, n_o);
  
  // Calculate vectors with random effects
  
  vector[n_f] r_eff;
  
  for(i in 1:n_f) {
    r_eff[i] = a_plot[plot[i]] + a_year[year[i]];
  }
  
  // Calculate final vital rates
  
  f_vr = crwd_c + crwd_h + crwd_o + b_l * l_f + r_eff + err;
}



