################################################################################
# Run script for testing crowding coefficients
#
# Written by Young Jun Lee
# Oct 2024

# Set run state to true
RUN_STATE = TRUE

# ===== Calculate crowding factors =====

pacman::p_load("here", "tidyverse")

# Unload stats packages that cause collision
pacman::p_unload("arm", "lme4", "MASS")

# Load required packages
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_packages.R"))

# Load census data
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))

# Load species vital rates
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_sp_vr.R"))

# Range of crowding function coefficients
# crowding_func_coeffs <- seq(0, 2, 0.025)
crowding_func_coeffs <- seq(5, 50, 2.5)

# For each crowding coefficient:
for(cfunc_coeff in crowding_func_coeffs) {
  # Set ALPHA
  ALPHA <- cfunc_coeff
  
  cat(paste0("Calculating crowding coefficients for ALPHA = ", ALPHA, "\n"))
  
  # Print current time
  cat(paste0("CURRENT TIME: ", Sys.time(), "\n"))
  
  # Calculate crowding factors
  # source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_calc_crowding_coeff_test.R"))
  source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_calc_crowding_trunc_coeff_test.R"))
}

# ====== Process data for stats =======

# For each crowding coefficient:
for(cfunc_coeff in crowding_func_coeffs) {
  # Set ALPHA
  ALPHA <- cfunc_coeff
  
  cat(paste0("Stat processing for ALPHA = ", ALPHA, "\n"))
  
  # Print current time
  cat(paste0("CURRENT TIME: ", Sys.time(), "\n"))
  
  # Calculate stats data
  source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_process_stat_coeff.R"))
}

# ===== Fit GLMMs =====

# Load stats packages
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
            "census_analysis_load_stat_packages.R"))

# Load stats data as dummy
source(here("StageDependentComp","code","census_analysis","sub_scripts",
            "census_analysis_load_stat_data.R"))

# Load stats functions
source(here("StageDependentComp","code","census_analysis","func",
            "stat_fitting.R"))

# For each crowding coefficient:
for(cfunc_coeff in crowding_func_coeffs) {
  # Set ALPHA
  ALPHA <- cfunc_coeff
  
  cat(paste0("Fitting GLMMs for ALPHA = ", ALPHA, "\n"))
  
  # Print current time
  cat(paste0("CURRENT TIME: ", Sys.time(), "\n"))
  
  # Load stat data
  source(here("StageDependentComp","code","census_analysis","sub_scripts",
              "census_analysis_load_stat_data_coeff.R"))
  
  # Fit frequentist models
  MODEL_TYPE <- "lmm"
  source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_coeff_test.R"))
  MODEL_TYPE <- "glmm"
  source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_coeff_test.R"))
  # source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm_nontrans_coeff_test.R"))
}

# Set ALPHA
ALPHA <- 5

cat(paste0("Fitting GLMMs for ALPHA = ", ALPHA, "\n"))

# Print current time
cat(paste0("CURRENT TIME: ", Sys.time(), "\n"))

# Load stat data
source(here("StageDependentComp","code","census_analysis","sub_scripts",
            "census_analysis_load_stat_data_coeff.R"))

# Fit frequentist models
MODEL_TYPE <- "lmm"
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_coeff_test.R"))
MODEL_TYPE <- "glmm"
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_coeff_test.R"))
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm_nontrans_coeff_test.R"


# Create summary comparing alpha values
MODEL_TYPE <- "lmm"

source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_summarise_coeff_test.R"))

MODEL_TYPE <- "glmm"

source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_summarise_coeff_test.R"))

# 
# # Fit frequentist models
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm_coeff_test.R"))
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm_nontrans_coeff_test.R"))
# 
# # Calculate crowding factors
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_calc_crowding.R"))
# 
# # Fit frequentist models
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm.R"))
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_stat_glmm_nontrans.R"))