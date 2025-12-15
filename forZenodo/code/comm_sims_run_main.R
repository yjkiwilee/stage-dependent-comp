################################################################################
#' Run script for two-species community simulations
#'
#' Written by Young Jun Lee
#' Mar 2025

# Set run state to true
RUN_STATE = TRUE

# Check if pacman package is installed. If not, install it.
if(system.file(package="pacman") == "") {
  install.packages("pacman")
}

# Set the working directory to the location where you have unzipped the code and data archive
setwd("/Users/chrissy/Oxford_Postdoc/MBiol_JunLee/stage-dependent-comp/forZenodo/")

# If the results directories do not exist, then make them:
if (!dir.exists(file.path("result_data"))){
  dir.create("result_data")
}
if (!dir.exists(file.path("figures"))){
  dir.create("figures")
}

# Load packages
source(file.path("code","sub_scripts","comm_sims_init_setup.R"))

# Specify models to use
source(file.path("code","sub_scripts","comm_sims_def_models.R"))

# Sample time series for each level of stage-dependence
source(file.path("code","sub_scripts","comm_sims_sample_tseries.R"))

# Fit regression models
source(file.path("code","sub_scripts","comm_sims_ddvr_est.R"))

# Generate predictions & calculate model errors
source(file.path("code","sub_scripts","comm_sims_proj_from_est_long.R"))

# Calculate statistical summaries
source(file.path("code","sub_scripts","comm_sims_calc_stats.R"))

# Process statistical summaries
source(file.path("code","sub_scripts","comm_sims_proc_stats.R"))

# Plot visualisations
source(file.path("code","sub_scripts","comm_sims_vis.R"))

# Carry out meta-statistical analyses
source(file.path("code","sub_scripts","comm_sims_meta_stats.R"))

# Make the table of species vital rates
source(file.path("code", "sub_scripts","comm_sims_tables.R" ))

