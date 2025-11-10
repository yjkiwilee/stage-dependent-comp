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

# Load config parameters
source(here("StageDependentComp","code","census_analysis","census_analysis_config.R"))

# Load packages
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_packages.R"))

# Setup ggplot
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_setup_ggplot.R"))

# Process raw environmental data
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_process_env.R"))

# Load raw environmental data
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_ann_env_data.R"))

# Load raw census data
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_data.R"))

# Process census data & save output
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_process_census.R"))

# Run statistical analyses

# FOCAL_VR <- "X..Capitulescences"
FOCAL_VRS <- c(
  "Survival.future",
  # "Growth.future",
  "Is.flowering"
  # "X..Capitulescences"
)

FOCAL_VR <- NULL

SD_SDS <- c(
  0.01,
  0.01,
  0.01,
  0.01
)

SD_SD <- NULL

SD_SPPS <- c(
  0.0002,
  0.0006,
  0.0005,
  0.0006
)

SD_SPP <- NULL

CHAIN_IDS <- c(
  1, 2, 3
)

CHAIN_ID <- NULL

for(i in 1:length(FOCAL_VRS)) {
  for(cid in CHAIN_IDS) {
    FOCAL_VR <- FOCAL_VRS[i]
    SD_SD <- SD_SDS[i]
    SD_SPP <- SD_SPPS[i]
    CHAIN_ID <- cid
    
    jobRunScript(
      here("StageDependentComp","code","census_analysis","sub_scripts",
           "census_analysis_overlap_mcmc_mult.R"),
      name = paste(FOCAL_VR, CHAIN_ID),
      importEnv = TRUE
    )
  }
}

# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_overlap_mcmc_mult.R"))

# Visualise analysis
# source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_vis.R"))


