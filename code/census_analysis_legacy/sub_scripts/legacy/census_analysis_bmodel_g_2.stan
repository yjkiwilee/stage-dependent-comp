// STAN script for testing the use of a Bayesian model to determine stage-dependence
// This model is to be used for continuous growth.
// Written by Young Jun Lee
// Feb 2024

data {
  // Data attributes
  int n; // Total number of observations
  int n_c; // Number of conspecifics to calculate crowding effect from
  int n_h; // Number of specific competitors to calculate crowding effect from
  int n_o; // Number of other competitors to calculate crowding effect from
  int n_tot; // n_c + n_h + n_o
  int n_plots; // Number of plots
  int n_years; // Number of years
  
  // Plant lengths
  vector[n] l_f;
  vector[n_tot] l_all;
  // vector[n_c] l_c;
  // vector[n_h] l_h;
  // vector[n_o] l_o;
  
  // Distance matrix
  matrix[n, n_tot] dist;
  // matrix[n, n_c] dist_c;
  // matrix[n, n_h] dist_h;
  // matrix[n, n_o] dist_o;
  
  // Response vital rate
  vector[n] vr;
  
  // Random effects; must be 1 ~ n_plots/n_years
  int plot[n];
  int year[n];
  
  // Parameters to control priors
  
  // Biologically reasonable mean value for shape parameter of kernel
  real exp_mean_d_t;
  
  // Intercept without random effect of year & plot
  real rand_intercept;
  // Random effects of year & plot
  vector[n_plots] plot_effs;
  vector[n_years] year_effs;
  
  // Vital rate sd
  real vr_sd;
}

parameters {
  // See daily log 08032025 for notation
  
  // Shape parameter for decay kernel
  real<lower=0, upper=200> d_t;
  
  // Intercepts
  real a_c;
  real a_h;
  real a_o;
  real a;
  
  // Slopes
  real b_c;
  real b_h;
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
  
  d_t ~ normal(exp_mean_d_t, 10);
  
  a_c ~ normal(0, 1);
  a_h ~ normal(0, 1);
  a_o ~ normal(0, 1);
  
  // Use intercept after random effect
  a ~ normal(rand_intercept, 1);
  
  b_c ~ normal(0, 1);
  b_h ~ normal(0, 1);
  b_o ~ normal(0, 1);
  b_l ~ normal(0, 1);
  
  s ~ normal(vr_sd, 1);
  
  // Random effect priors
  a_plot ~ normal(plot_effs, rep_vector(1, n_plots));
  a_year ~ normal(year_effs, rep_vector(1, n_years));
  
  // Transform distance matrices with Gaussian decay kernel
  
  matrix[n, n_tot] dist_t;
  dist_t = exp((-0.30103 / d_t^2) .* (dist .* dist));
  
  // Calculate contribution from crowding
  
  row_vector[n_tot] a_crwd;
  row_vector[n_tot] b_crwd;
  row_vector[n_tot] bo_crwd;
  row_vector[n_tot] bi_crwd;
  
  a_crwd = append_col(
    rep_row_vector(a_c, n_c),
    append_col(
      rep_row_vector(a_h, n_h),
      rep_row_vector(a_o, n_o)
    )
  );
  
  b_crwd = append_col(
    rep_row_vector(b_c, n_c),
    append_col(
      rep_row_vector(b_h, n_h),
      rep_row_vector(b_o, n_o)
    )
  );
  
  vector[n] ones_n;
  ones_n = rep_vector(1, n);
  row_vector[n_tot] ones_n_tot;
  ones_n_tot = rep_row_vector(1, n_tot);
  
  vector[n] crwd;
  
  crwd = (dist_t .* (
    (ones_n * a_crwd) +
    (ones_n * (l_all' .* b_crwd))
  )) * ones_n_tot';
  
  // vector[n] crwd_c;
  // vector[n] crwd_h;
  // vector[n] crwd_o;
  // 
  // crwd_c = (dist_c_t .* (
  //   a_c +
  //   (b_c * (l_f * rep_row_vector(1, n_c))) +
  //   (b_co * (rep_vector(1, n) * l_c')) +
  //   (b_ci * (l_f * l_c'))
  // )) * rep_vector(1, n_c);
  // 
  // crwd_h = (dist_h_t .* (
  //   a_h +
  //   (b_h * (l_f * rep_row_vector(1, n_h))) +
  //   (b_ho * (rep_vector(1, n) * l_h')) +
  //   (b_hi * (l_f * l_h'))
  // )) * rep_vector(1, n_h);
  // 
  // crwd_o = (dist_o_t .* (
  //   a_o +
  //   (b_o * (l_f * rep_row_vector(1, n_o))) +
  //   (b_oo * (rep_vector(1, n) * l_o')) +
  //   (b_oi * (l_f * l_o'))
  // )) * rep_vector(1, n_o);
  
  // Calculate vectors with random effects
  
  vector[n] r_eff;
  
  for(i in 1:n) {
    r_eff[i] = a_plot[plot[i]] + a_year[year[i]];
  }
  
  // Calculate final vital rates
  vector[n] mu;
  mu = crwd + (b_l * l_f) + r_eff + a;
  vr ~ normal(mu, s);
}



