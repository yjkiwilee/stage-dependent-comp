#'
#' This version is 'simplified' in the sense that the 'reproduction factor' has been removed.
#' This is because fecundity and reproduction factor can be collapsed into one.
#'
#' TODO
#' - Implement multidimensional k_mat altering function
#' - Implement function for adjusting vital rates while keeping lambda constant

# Number of timesteps used for long-term behaviour
SIMEQ_TIMESTEP <- 1000
# Number of timesteps used for robust k_mat estimation
SIMEQ_TIMESTEP_ROBUST <- 7000
# Number of timesteps to be averaged for getting asymptotic behaviour of stochastic / cyclic / chaotic models
SIMEQ_AVGTF <- 5000
# Minimum stage abundance
MIN_SIZE <- 1
# Default population structure
DEFAULT_POP_STRUCT <- c(100, 100)

# Import sdcompr_def_model.R with functions that define a model: DONE SIMPLIFYING
source(here("StageDependentComp", "code", "comm_sims", "func_simple", "sdcompr_def_model.R"))

# Import sdcompr_mpm.R with functions related to generating Rage mpms: DONE SIMPLIFYING
source(here("StageDependentComp", "code", "comm_sims", "func_simple", "sdcompr_mpm.R"))

# Import sdcompr_projection.R with functions that project populations using a sdcomp model: DONE SIMPLIFYING
source(here("StageDependentComp", "code", "comm_sims", "func_simple", "sdcompr_projection.R"))

# Import sdcompr_long.R with functions for long-term dynamics analysis: DONE SIMPLIFYING
source(here("StageDependentComp", "code", "comm_sims", "func_simple", "sdcompr_long.R"))

# Import sdcompr_comp_manip.R with functions for manipulating competitive coefficients: DONE SIMPLIFYING
source(here("StageDependentComp", "code", "comm_sims", "func_simple", "sdcompr_comp_manip.R"))

# Import sdcompr_vr_manip.R with functions for manipulating vital rates
source(here("StageDependentComp", "code", "comm_sims", "func_simple", "sdcompr_vr_manip.R"))





