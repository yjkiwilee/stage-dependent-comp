# TODO
# - The current alter_k_mat function is ill equipped for situations where simulations 'blow up'.
# This can be manually circumvented if the scal_fact_l / scal_fact_r range is made narrower, but could I implement automatic determination of those limit?
# - The alter_k_mat function is also limited where the asymptotic behaviour is a cycle or chaotic.

# Number of timesteps used for long-term behaviour
SIMEQ_TIMESTEP <- 400
# Number of timesteps to be averaged for getting asymptotic behaviour of stochastic models
SIMEQ_AVGTF <- 100
# Minimum stage size, below which sizes are set to 0
MIN_SIZE <- 0.1
# Default population structure
DEFAULT_POP_STRUCT <- c(100, 100)

#' Define a pre-breeding stage-dependent competition model
#'
#' This function defines a stage-dependent competition model with a set of density-independent
#' vital rates and competitive parameters (provided as a matrix). This model includes
#' seed survival. This model currently includes sexual reproduction only and does not include shrinkage.
#' The model also assumes that individuals of the smallest stage are non-reproductive.
#' With n stages, the vital rates of each species should include: s (survival probabilities including seed survival, n+1),
#' g (maturation probabilities, n-1), and f (seed production, n-1).
#' With n stages, the competitive coefficient matrices is a list that contains s, g, and f.
#' Each contain the competitive coefficients controlling survival, maturation or fecundity, respectively.
#' Each are given as a 2m x 2n matrix, where m is the number of corresponding vital rates per species.
#' The final column in the survival and maturation matrices are the intercepts.
#' The vital rates are ordered by species and then by stage. Positive value denotes negative effect.
#'
#' @param n_stages Number of stages (n) in the model
#' @param sp1_vrs List containing the density-independent vital rates of species 1
#' @param sp2_vrs List containing the density-independent vital rates of species 2
#' @param k_mat Matrix containing competitive coefficients
#' @param k_mat_sd Matrix containing the standard deviations for the competitive coefficients (assumed to be normally distributed)
#' 
#' @return The pre-breeding stage-dependent competition model, given as a list
#' 
#' @export
def_sdcomp_model <- function(n_stages, sp1_vrs, sp2_vrs, k_mat, k_mat_sd = NULL) {
  # Calculate intercepts from vrs for each vital rate type
  for(vr_type in names(k_mat)) {
    is_logistic <- (vr_type == "s" | vr_type == "g")
    
    comm_vrs <- matrix(c(sp1_vrs[[vr_type]], sp2_vrs[[vr_type]]), ncol = 1)
    
    intrcpt <- matrix()
    
    if(is_logistic) {
      intrcpt <- -qlogis(comm_vrs) # Inverse function of negative logistic
    } else {
      intrcpt <- -log(comm_vrs) # Inverse function of negative exponential
    }
    
    # Bind intercepts to k_mat
    k_mat[[vr_type]] <- cbind(k_mat[[vr_type]], intrcpt)
  }
  
  # Store model (NB: k_mat_sd will be NULL if not specified)
  sdcomp_model <- list(
    n_stages = n_stages,
    sp1_vrs = sp1_vrs,
    sp2_vrs = sp2_vrs,
    k_mat = k_mat,
    k_mat_sd = k_mat_sd
  )
  
  # Return model
  sdcomp_model
}

#' Define a pre-breeding stage-dependent competition model with pre-defined intercepts
#'
#' This function is the same as def_sdcomp_model, except the k_mat is fully specified
#' with intercepts. This means that sp1_vrs and sp2_vrs are not necessary.
#' 
#' @param n_stages Number of stages (n) in the model
#' @param k_mat Matrix containing competitive coefficients
#' 
#' @return The pre-breeding stage-dependent competition model, given as a list
#' 
#' @export
def_sdcomp_model_intrcpt <- function(n_stages, k_mat) {
  # Lists to store species vital rates
  sp1_vrs <- list()
  sp2_vrs <- list()
  
  # Calculate vrs from intercepts
  for(vr_type in names(k_mat)) {
    is_logistic <- (vr_type == "s" | vr_type == "g")
    
    intrcpt <- k_mat[[vr_type]][, 2 * n_stages + 1]
    comm_vrs <- c()
    
    if(is_logistic) {
      comm_vrs <- plogis(-intrcpt) # Negative logistic function of intercept
    } else {
      comm_vrs <- exp(-intrcpt) # Negative exponential function of intercept
    }
    
    # Store into sp1 & sp2_vrs
    sp1_vrs[[vr_type]] <- as.vector(comm_vrs[1:(length(comm_vrs) / 2)])
    sp2_vrs[[vr_type]] <- as.vector(comm_vrs[(length(comm_vrs) / 2 + 1):length(comm_vrs)])
  }
  
  # List to store the model in (NB: k_mat_sd will be NULL if not specified)
  sdcomp_model <- list(
    n_stages = n_stages,
    sp1_vrs = sp1_vrs,
    sp2_vrs = sp2_vrs,
    k_mat = k_mat,
    k_mat_sd = k_mat_sd
  )
  
  # Return model
  sdcomp_model
}

#' Calculate a matrix population model from vital rates
#'
#' This function calculates an n-stage matrix population model from vital rates.
#' The vital rates are given as a list, which should include:
#' s (survival probabilities including seed survival, n+1),
#' g (maturation probabilities, n-1), and f (seed production, n-1).
#'
#' @param n_stages The number of stages
#' @param vrs Vital rates given in a list
#' 
#' @return An Rage matrix, given as a list
#'
#' @export
vrs_to_mpm <- function(n_stages, vrs) {
  # Extract vital rates
  surv <- vrs$s
  matur <- vrs$g
  fec <- vrs$f
  
  # Survival and fecundity matrices
  matU <- matrix(rep(0, n_stages^2), nrow = n_stages, ncol = n_stages)
  matF <- matrix(rep(0, n_stages^2), nrow = n_stages, ncol = n_stages)
  
  # ===== Calculate matU =====
  # For each stage
  for(i in 1:n_stages) {
    if(i == n_stages) { # If final stage
      # Stasis only
      matU[i, i] <- surv[i+1]
    } else { # If not final stage
      # Stasis
      matU[i, i] <- surv[i+1] * (1 - matur[i])
      # Maturation
      matU[i+1, i] <- surv[i+1] * matur[i]
    }
  }
  
  # ===== Calculate matF =====
  # For each reproductive stage
  for(i in 2:n_stages) {
    # Reproduction (with consideration of seed survival)
    matF[1, i] <- surv[1] * fec[i-1]
  }
  
  # Create & return mpm
  mpm <- list(
    matU = matU,
    matF = matF,
    matA = matU + matF
  )
  
  mpm
}

#' Calculate a two-species projection matrix from population structure
#'
#' This function calculates a two-species projection matrix from a stage-dependent competition model,
#' given some population structure. Probabilities use a logistic density-dependent function,
#' while fecundity uses a negative exponential function.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_vec Vector with the population structure of species 1
#' @param sp2_vec Vector with the population structure of species 2
#' 
#' @return An Rage two-species matrix, given as a list
#'
#' @export
sdcomp_to_mpm <- function(sdcomp_model, sp1_pop, sp2_pop) {
  # Number of stages
  n_stages <- sdcomp_model$n_stages
  
  # Is stochastic model?
  is_stoch <- !is.null(sdcomp_model$k_mat_sd)
  
  # Lists to store vital rates after density-dependent scaling
  sp1_dd_vrs <- list()
  sp2_dd_vrs <- list()
  
  # Competitive coefficients + intercepts
  k_mat <- sdcomp_model$k_mat
  
  # For each vital rate type, calculate the version after density-dependent scaling
  for(vr_type in names(k_mat)) {
    # Alter k_mat if stochastic
    if(is_stoch) {
      # Set k_mat according to a normal distribution
      k_mat[[vr_type]] <- matrix(rnorm(length(k_mat[[vr_type]]), mean = k_mat[[vr_type]], sd = k_mat_sd[[vr_type]]),
                                 ncol = 2 * n_stages + 1)
    }
    
    is_logistic <- (vr_type == "s" | vr_type == "g")
    
    # Convert sp1/2_pop to column vector
    # Add extra element if logistic for intercept
    comm_cv <- matrix(c(sp1_pop, sp2_pop, 1), ncol = 1)
    
    # Variable to store scaling factor for vital rates
    dd_vrs <- c()
    
    # Calculate density-dependent vital rates
    if(is_logistic) { # If survival or maturation prob.
      dd_vrs <- as.vector(plogis(-k_mat[[vr_type]] %*% comm_cv)) # Logistic scaling
    } else { # If fecundity
      dd_vrs <- as.vector(exp(-k_mat[[vr_type]] %*% comm_cv)) # Negative exponential scaling
    }
    
    # Store density-dependent vital rates
    sp1_dd_vrs[[vr_type]] <- dd_vrs[1:(length(dd_vrs) / 2)]
    sp2_dd_vrs[[vr_type]] <- dd_vrs[(length(dd_vrs) / 2 + 1):length(dd_vrs)]
  }
  
  # Construct matrix from vital rates
  sp1_mpm <- vrs_to_mpm(n_stages, sp1_dd_vrs)
  sp2_mpm <- vrs_to_mpm(n_stages, sp2_dd_vrs)
  
  # Combine matrices
  comm_mpm <- list()
  comm_mpm$matU <- cbind(
    rbind(sp1_mpm$matU, matrix(rep(0, n_stages^2), nrow = n_stages)),
    rbind(matrix(rep(0, n_stages^2), nrow = n_stages), sp2_mpm$matU)
  )
  comm_mpm$matF <- cbind(
    rbind(sp1_mpm$matF, matrix(rep(0, n_stages^2), nrow = n_stages)),
    rbind(matrix(rep(0, n_stages^2), nrow = n_stages), sp2_mpm$matF)
  )
  
  # Row/col names
  rc_names <- c(paste0("sp1_s", 1:n_stages), paste0("sp2_s", 1:n_stages))
  rownames(comm_mpm$matU) <- rc_names
  colnames(comm_mpm$matU) <- rc_names
  rownames(comm_mpm$matF) <- rc_names
  colnames(comm_mpm$matF) <- rc_names
  
  comm_mpm$matA <- comm_mpm$matU + comm_mpm$matF
  
  # Return MPM
  comm_mpm
}

#' Project the community one time step forward in time
#'
#' This function projects the community forward in time by one time step.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_vec Vector with the population structure of species 1
#' @param sp2_vec Vector with the population structure of species 2
#' @param min_size Minimum size of stage
#' 
#' @return An Rage two-species matrix, given as a list
#'
#' @export
sdcomp_project_step <- function(sdcomp_model, sp1_pop, sp2_pop, min_size = MIN_SIZE) {
  n_stages <- sdcomp_model$n_stages
  
  # Get MPM for timestep
  comm_mpm <- sdcomp_to_mpm(sdcomp_model, sp1_pop, sp2_pop)
  
  # Merge population vectors into column vector
  comm_cv <- matrix(c(sp1_pop, sp2_pop), ncol = 1)
  
  # Project by one time step
  comm_cv_f <- comm_mpm$matA %*% comm_cv
  comm_vec <- as.vector(comm_cv_f)
  
  # Set size as zero if less than one
  comm_vec[comm_vec < min_size] <- 0
  
  # Split and return as list
  sp1_pop_f <- comm_vec[1:n_stages]
  sp2_pop_f <- comm_vec[(n_stages + 1):(2 * n_stages)]
  list(
    sp1_pop = sp1_pop_f,
    sp2_pop = sp2_pop_f
  )
}

#' Project the community through time to get a time-series
#'
#' This function projects the community forward in time and returns the trajectory.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_0 Vector with the initial population structure of species 1
#' @param sp2_0 Vector with the initial population structure of species 2
#' @param timestep Number of time steps to simulate
#' 
#' @return A wide dataframe with columns for time and population structure of the two species
#'
#' @export
sdcomp_project <- function(sdcomp_model, sp1_0, sp2_0, timestep) {
  n_stages <- sdcomp_model$n_stages
  
  # Vector for storing time
  time_vec <- 0:timestep
  
  # Vector for storing current population structure
  sp1_curr <- sp1_0
  sp2_curr <- sp2_0
  
  # Matrix for storing the population structures through time
  sp1_mat <- matrix(sp1_0, ncol = n_stages)
  sp2_mat <- matrix(sp2_0, ncol = n_stages)
  
  # Simulate
  for(t in 1:timestep) {
    # Project by one time step
    comm_curr <- sdcomp_project_step(sdcomp_model, sp1_curr, sp2_curr)
    sp1_curr <- comm_curr$sp1_pop
    sp2_curr <- comm_curr$sp2_pop
    
    # Store population structure
    sp1_mat <- rbind(sp1_mat, matrix(sp1_curr, ncol = n_stages))
    sp2_mat <- rbind(sp2_mat, matrix(sp2_curr, ncol = n_stages))
  }
  
  # Bind population trajectories
  comm_mat <- cbind(sp1_mat, sp2_mat)
  
  # Rename columns and coerce into dataframe
  colnames(comm_mat) <- c(paste0("sp1_s", 1:n_stages), paste0("sp2_s", 1:n_stages))
  res_df <- data.frame(comm_mat)
  
  # Put time into dataframe
  res_df$time <- time_vec
  
  # Return
  res_df
}

#' Project the community through time to get the end state
#' 
#' This function projects the community forward in time and returns the end-state.
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_0 Vector with the initial population structure of species 1
#' @param sp2_0 Vector with the initial population structure of species 2
#' @param timestep Number of time steps to simulate
#' @param stochmean Boolean determining whether we want to retrieve the averaged 'stochastic mean' or the true end state at timestep
#' @param avgtf Number of time steps before 'timestep' over which the population sizes will be averaged (set as half of timestep if not specified)
#' 
#' @return A list containing the end-state of each of the population
#' 
#' @export
sdcomp_project_end <- function(sdcomp_model, sp1_0, sp2_0, timestep) {
  n_stages <- sdcomp_model$n_stages
  
  # Get simulation trajectory using sdcomp_project
  sim_traj <- sdcomp_project(sdcomp_model, sp1_0, sp2_0, timestep)
  sp1_traj <- sim_traj %>%
    select(starts_with("sp1"))
  sp2_traj <- sim_traj %>%
    select(starts_with("sp2"))
  
  return( # Return last timestep
    list(
      sp1_pop = as.vector(t(sp1_traj[nrow(sp1_traj),])),
      sp2_pop = as.vector(t(sp2_traj[nrow(sp2_traj),]))
    )
  )
}

#' Project the community through time to get the average across the last n timesteps
#' 
#' Identical to sdcomp_project_end, but more suited for stochastic simulations
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_0 Vector with the initial population structure of species 1
#' @param sp2_0 Vector with the initial population structure of species 2
#' @param timestep Number of time steps to simulate
#' @param avgtf Number of time steps before 'timestep' over which the population sizes will be averaged (set as half of timestep if not specified)
#' 
#' @return A list containing the averaged end-state of each of the population
#' 
#' @export
sdcomp_project_end_avg <- function(sdcomp_model, sp1_0, sp2_0, timestep, avgtf = NULL) {
  n_stages <- sdcomp_model$n_stages
  
  if(is.null(avgtf)) { avgtf <- timestep %/% 2 } # Set avgtf as half of timestep if not provided
  
  # Get simulation trajectory using sdcomp_project
  sim_traj <- sdcomp_project(sdcomp_model, sp1_0, sp2_0, timestep)
  sim_traj_end <- sim_traj[(nrow(sim_traj) - avgtf + 1):nrow(sim_traj),] # Subset last 'avgtf' rows
  sp1_traj <- sim_traj_end %>%
    select(starts_with("sp1"))
  sp2_traj <- sim_traj_end %>%
    select(starts_with("sp2"))
  
  # Get average
  sp1_pop_avg <- as.vector(t(colMeans(sp1_traj)))
  sp2_pop_avg <- as.vector(t(colMeans(sp2_traj)))
  
  # Return as list
  list(
    sp1_pop = sp1_pop_avg,
    sp2_pop = sp2_pop_avg
  )
}

#' Determine whether stable coexistence or competitive exclusion occurs
#' 
#' This function either simulates or takes in a time-series and determine the long-term behaviour of the system.
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_0 Vector with the initial population structure of species 1 (Default = DEFAULT_POP_STRUCT)
#' @param sp2_0 Vector with the initial population structure of species 2 (Default = DEFAULT_POP_STRUCT)
#' @param timestep Number of time steps to simulate
#' @param proj_end End-state data; if this is not NULL, the other arguments are ignored
#' 
#' @return A boolean value, TRUE for stable coexistence and FALSE for competitive exclusion
#' 
#' @export
does_coexist <- function(sdcomp_model, sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT, timestep = SIMEQ_TIMESTEP, proj_end = NULL) {
  if(is.null(proj_data)) { # If the projection data is not already given
    proj_end <- sdcomp_project_end(sdcomp_model, sp1_0, sp2_0, timestep) # Populate projection data
  }
  
  if(sum(proj_end$sp1_pop) == 0 | sum(proj_end$sp2_pop) == 0) { # Competitive exclusion
    return(FALSE)
  } else { # Stable coexistence
    return(TRUE)
  }
}

#' Calculate average per-capita interaction strengths of a species given a sdcomp model
#' 
#' This function calculates the average per-capita competitive strength of a species on a specific vital rate at equilibrium.
#' NOT COMPATIBLE WITH STOCHASTIC MODEL
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param vr_type Type of vital rate being affected
#' @param mat_row Row of matrix corresponding to vital rate being affected
#' @param species Species exerting the interactive effect. 1 = sp1, 2 = sp2
#' @param timestep Number of time steps to simulate to reach equilibrium (Default = SIMEQ_TIMESTEP)
#' @param sp1_0 Initial population structure for species 1 (Default = DEFAULT_POP_STRUCT)
#' @param sp2_0 Initial population structure for species 2 (Default = DEFAULT_POP_STRUCT)
#' 
#' @return The average per-capita competitive strength or NA if competitive exclusion occurs
#' 
#' @export
calc_avg_int_str <- function(sdcomp_model, vr_type, mat_row, species, timestep = SIMEQ_TIMESTEP, sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT) {
  n_stages <- sdcomp_model$n_stages
  
  # Remove stochastic component
  k_mat_sd <- sdcomp_model$k_mat_sd
  sdcomp_model$k_mat_sd <- NULL
  
  # Find 'equilibrium' of species with current parameters
  end_state <- sdcomp_project_end(sdcomp_model, sp1_0, sp2_0, timestep)
  
  end_state_sp <- end_state[[ifelse(species == 1, "sp1_pop", "sp2_pop")]] # End state for the specific species
  # Escape if the focal species goes extinct
  if(sum(end_state_sp) <= 0) { return(NA) }
  
  # Normalise end state
  end_state_sp_norm <- end_state_sp / sum(end_state_sp)
  
  # Determine the current per-capia interaction strength of the species on the specific vital rate
  mat_cols <- c((1:n_stages) + ((species - 1) * n_stages))
  k_mat_elems <- sdcomp_model$k_mat[[vr_type]][mat_row, mat_cols] # Extract k_mat elements for the species & vital rate
  avg_int_str <- sum(end_state_sp_norm * k_mat_elems) # Average per-capita interaction strength
  
  return(avg_int_str)
}

#' Adjust a competitive coefficient while keeping species-level per-capita interaction strength equal
#' 
#' This function takes in a sdcomp model, and a competitive coefficient to be adjusted.
#' It then returns a new sdcomp model where the per-capita interaction strength at the stable population structure is equal.
#' This function uses the bisection method to approximate the appropriate competitive coefficient values.
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param vr_type Type of vital rate being affected
#' @param mat_row Row of matrix element to alter
#' @param mat_col Column of matrix element to alter
#' @param alt_val Alternate value of the competitive coefficient
#' @param timestep Number of time steps to simulate to reach equilibrium (Default = SIMEQ_TIMESTEP)
#' @param iter Number of iterations used to approsimate the solution
#' @param scal_fact_l Initial minimum scale factor to try for bisection
#' @param scal_fact_r Initial minimum scale factor to try for bisection
#' @param sp1_0 Initial population structure of species 1
#' @param sp2_0 Initial population structure of species 2
#' 
#' @return A sdcomp model with the adjustment controlling for average per-capita interaction strength
#' 
#' @export
alter_k_mat <- function(sdcomp_model, vr_type, mat_row, mat_col, alt_val, timestep = SIMEQ_TIMESTEP, avgtf = SIMEQ_AVGTF, iter = 40, scal_fact_l = 0, scal_fact_r = 5, sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT) {
  n_stages <- sdcomp_model$n_stages
  species <- (mat_col - 1) %/% n_stages + 1 # 1 if species 1, 2 if species 2
  
  # Determine the target per-capita interaction strength
  targ_int_str <- calc_avg_int_str(sdcomp_model, vr_type, mat_row, species, timestep = timestep, sp1_0 = sp1_0, sp2_0 = sp2_0)
  
  # Escape if the focal species goes extinct
  if(is.na(targ_int_str)) {
    warning("Focal species extinction")
    return(NA)
  }
  
  # Identify the columns of elements to be adjusted to accommodate the altered value
  adj_cols <- (1:n_stages) + ((species - 1) * n_stages)
  adj_cols <- adj_cols[adj_cols != mat_col]
  
  # Variable to store the appropriate scaling factor for the elements to be adjusted
  scal_fact_m <- mean(c(scal_fact_l, scal_fact_r))
  
  # Alter k_mat element
  sdcomp_model$k_mat[[vr_type]][mat_row, mat_col] <- alt_val
  
  # Variables to store sdcomp for left and right bounds
  l_model <- NA
  m_model <- NA
  r_model <- NA
  
  # Bisection method to find root
  for(i in 1:iter) {
    l_model <- sdcomp_model
    m_model <- sdcomp_model
    r_model <- sdcomp_model
    
    l_model$k_mat[[vr_type]][mat_row, adj_cols] <- l_model$k_mat[[vr_type]][mat_row, adj_cols] * scal_fact_l
    r_model$k_mat[[vr_type]][mat_row, adj_cols] <- r_model$k_mat[[vr_type]][mat_row, adj_cols] * scal_fact_r
    m_model$k_mat[[vr_type]][mat_row, adj_cols] <- m_model$k_mat[[vr_type]][mat_row, adj_cols] * scal_fact_m
    
    l_int_str <- calc_avg_int_str(l_model, vr_type, mat_row, species, timestep, sp1_0, sp2_0)
    r_int_str <- calc_avg_int_str(r_model, vr_type, mat_row, species, timestep, sp1_0, sp2_0)
    m_int_str <- calc_avg_int_str(m_model, vr_type, mat_row, species, timestep, sp1_0, sp2_0)
    
    if(l_int_str > targ_int_str | r_int_str < targ_int_str) {
      warning("Failed to converge")
      return(NA)
    } # Failing to converge, panic break...
    else if(m_int_str < targ_int_str) { scal_fact_l <- scal_fact_m }
    else { scal_fact_r <- scal_fact_m }
    
    scal_fact_m <- mean(c(scal_fact_l, scal_fact_r))
  }
  
  # Return resulting model
  return(m_model)
}








