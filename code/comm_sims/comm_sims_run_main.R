################################################################################
#' Run script for two-species community simulations
#'
#' Written by Young Jun Lee
#' Mar 2025

# Set run state to true
RUN_STATE = TRUE

# Install here library
if(system.file(package="here") == "") {
  install.packages("here")
}
if(system.file(package="pacman") == "") {
  install.packages("pacman")
}
require("here")

# Load packages
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_init_setup.R"))

# Specify models to use
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_def_models.R"))

# Sample time series for each level of stage-dependence
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_sample_tseries.R"))

# Fit regression models
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_ddvr_est.R"))

# Generate predictions & calculate model errors
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_proj_from_est_long.R"))

# Calculate statistical summaries
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_calc_stats.R"))

# Plot visualisations
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_vis.R"))

# Carry out meta-statistical analyses
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_meta_stats.R"))
# 
# # Q2.1: Portion of variance explained by population size vs structure
# source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_var_decomp.R"))
# 
# # H2.3: S-D vs. probability of stochastic extinction
# source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_sd_prob_ext.R"))
# 
# # H2.4: Interaction with life history
# source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_sd_life_history.R"))



