################################################################################
# Script for loading the packages needed for census data analysis as part of WP2.
#
# Written by Young Jun Lee
# Oct 2024

# Check if pacman is installed, if not, install pacman
if(system.file(package="pacman") == "") {
  install.packages("pacman")
}

# Install & load the rest of the packages
pacman::p_load(
  "here", # For generating file paths consistent across OS
  "popbio", # For demographic calculations
  "tidyverse", # For plotting & data manipulation
  "ggforce", # Additional features for ggplot
  "patchwork", # For merging plots together
  "RColorBrewer", # For Brewer palette in ggplot
  "parallel", # For parallel processing
  "boot", # For bootstrapping analysis
  "lme4", # For GLMM fitting
  "lmtest", # For testing linear hypotheses
  "car", # Extension for lme4
  "rstudioapi" # For running script as background job
)
