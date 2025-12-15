################################################################################
#' Load package for community simulations & do initial setup
#'
#' Written by Young Jun Lee
#' Mar 2025
#' 

# Load necessary packages
pacman::p_load("tidyverse", "Rage", "here", "ggpubr", "popbio", "tseries",
               "sarima", "lme4", "car", "patchwork", "cowplot", "flextable",
               "brm")

# Load ggplot themes
source(file.path("code","custom_themes.R"))

# Load sdcompr functions
source(file.path("code","func_simple","sdcompr.R"))

# Load comm-sims run script functions
source(file.path("code", "sub_scripts", "comm_sims_funcs.R"))

# Set default figure size
FIG_W <- 6
FIG_H <- 4

# Number of life stages
N_STAGES <- 2

