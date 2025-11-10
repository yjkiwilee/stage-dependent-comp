#' Project the community through time to get the end state
#' 
#' This function projects the community forward in time and returns the end-state.
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_0 Vector with the initial population structure of species 1
#' @param sp2_0 Vector with the initial population structure of species 2
#' @param timestep Number of time steps to simulate
#' @param stochmean Boolean determining whether we want to retrieve the averaged 'stochastic mean' or the true end state at timestep
#' 
#' @return A list containing the end-state of each of the population
#' 
#' @export
sdcomp_project_end <- function(sdcomp_model, sp1_0, sp2_0, timestep, dstb_func = NULL, vrstoch_func = NULL) {
  n_stages <- sdcomp_model$n_stages
  
  # Get simulation trajectory using sdcomp_project
  sim_traj <- NULL
  sim_traj <- sdcomp_project_dstb(sdcomp_model, sp1_0, sp2_0, timestep, dstb_func, vrstoch_func)
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
#' @param dstb_func Disturbance function
#' @param vrstoch_func Vital rate stochasticity function
#' 
#' @return A list containing the averaged end-state of each of the population
#' 
#' @export
sdcomp_project_end_avg <- function(sdcomp_model, sp1_0, sp2_0, timestep, avgtf = NULL, dstb_func = NULL, vrstoch_func = NULL) {
  n_stages <- sdcomp_model$n_stages
  
  if(is.null(avgtf)) { avgtf <- timestep %/% 2 } # Set avgtf as half of timestep if not provided
  
  # Get simulation trajectory using sdcomp_projects
  sim_traj <- NULL
  sim_traj <- sdcomp_project_dstb(sdcomp_model, sp1_0, sp2_0, timestep, dstb_func, vrstoch_func)
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
#' @param dstb_func Disturbance function
#' @param vrstoch_func Vital rate stochasticity function
#' 
#' @return A boolean value, TRUE for stable coexistence and FALSE for competitive exclusion
#' 
#' @export
does_coexist <- function(sdcomp_model, sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT, timestep = SIMEQ_TIMESTEP, proj_end = NULL, dstb_func = NULL, vrstoch_func = NULL) {
  if(is.null(proj_data)) { # If the projection data is not already given
    proj_end <- sdcomp_project_end(sdcomp_model, sp1_0, sp2_0, timestep, dstb_func, vrstoch_func) # Populate projection data
  }
  
  if(sum(proj_end$sp1_pop) == 0 | sum(proj_end$sp2_pop) == 0) { # Competitive exclusion
    return(FALSE)
  } else { # Stable coexistence
    return(TRUE)
  }
}