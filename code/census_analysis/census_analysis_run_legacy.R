################################################################################
# Run script for census analysis
#
# Written by Young Jun Lee
# Oct 2024

# Set run state to true
RUN_STATE = TRUE

# Install here library
if(system.file(package="here") == "") {
  install.packages("here")
}
library("here")

# Load packages
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_packages.R"))

# Load data
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_data.R"))

# Process census data & save output
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_process_census.R"))

# Set focal species
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_set_sp.R"))

# Setup ggplot
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_setup_ggplot.R"))

# Plot vital rates of focal species
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_plot_vrs.R"))

# Calculate species stage & vital rates
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_calc_stage.R"))

# Calculate recruitment factor for species
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_calc_rec_fact.R"))

# Calculate matrices for species
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_calc_mat.R"))

# Process data for stats
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_process_stat.R"))

# Fit models for a range of crowding factors
source(here("StageDependentComp","code","census_analysis","census_analysis_run_crowding_coeff.R"))

# # Calculate crowding factors
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_calc_crowding.R"))
# 
# # Fit frequentist models
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm.R"))
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm_nontrans.R"))



