#' TODO
#' - Implement multidimensional k_mat altering function
#' - Implement function for adjusting vital rates while keeping lambda constant

# Number of timesteps used for long-term behaviour
SIMEQ_TIMESTEP <- 400
# Number of timesteps used for robust k_mat estimation
SIMEQ_TIMESTEP_ROBUST <- 1000
# Number of timesteps to be averaged for getting asymptotic behaviour of stochastic / cyclic / chaotic models
SIMEQ_AVGTF <- 700
# Minimum stage density, below which density is set to 0
MIN_SIZE <- 0.1
# Default population structure
DEFAULT_POP_STRUCT <- c(100, 100)

# Import sdcompr_def_model.R with functions that define a model
source(here("StageDependentComp", "code", "comm_sims", "func", "sdcompr_def_model.R"))

# Import sdcompr_mpm.R with functions related to generating Rage mpms
source(here("StageDependentComp", "code", "comm_sims", "func", "sdcompr_mpm.R"))

# Import sdcompr_projection.R with functions that project populations using a sdcomp model
source(here("StageDependentComp", "code", "comm_sims", "func", "sdcompr_projection.R"))

# Import sdcompr_long.R with functions for long-term dynamics analysis
source(here("StageDependentComp", "code", "comm_sims", "func", "sdcompr_long.R"))

# Import sdcompr_comp_manip.R with functions for manipulating competitive coefficients
source(here("StageDependentComp", "code", "comm_sims", "func", "sdcompr_comp_manip.R"))

# Import sdcompr_vr_manip.R with functions for manipulating vital rates
source(here("StageDependentComp", "code", "comm_sims", "func", "sdcompr_vr_manip.R"))





