################################################################################
#' Run script for census analysis
#'
#' Written by Young Jun Lee
#' Apr 2025

# Set run state to true
RUN_STATE = TRUE

# Install here library
if(system.file(package="here") == "") {
  install.packages("here")
}
library("here")

# Load packages
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_packages.R"))
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_stat_packages.R"))

# Load data
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_data.R"))

# Process census data & save output
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_process_census.R"))

# Setup ggplot
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_setup_ggplot.R"))

# Run statistical analyses
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_overlap.R"))
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_merge_bootstrap.R"))

# Visualise analysis
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_vis.R"))