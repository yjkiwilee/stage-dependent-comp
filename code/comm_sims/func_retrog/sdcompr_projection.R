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

# ===== Functions handling stochastic simulations =====

# The functions in this file are dependent on functions and variables in sdcompr_det.R.

#' Project the community through time to get a time-series, with a time-dependent absolute 'disturbance function'
#'
#' This function projects the community forward in time and returns the trajectory.
#' At each timestep, the timestep is fed into the disturbance function, which specifies the amount by which
#' the sizes of each stage should be altered after projection.
#' This function also takes in vrstoch_func, which specifies the new stochastic density-independent vital rates.
#' NB: vrstoch_func sees the population structure BEFORE the dstb_func is applied.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_0 Vector with the initial population structure of species 1
#' @param sp2_0 Vector with the initial population structure of species 2
#' @param timestep Number of time steps to simulate
#' @param dstb_func Function that takes in the current timestep and population sizes and returns a list (w/ sp1 & sp2) of vectors that specify absolute 'disturbance'
#' @param vrstoch_func Function that takes in the current timestep, the model, and population sizes and returns a list (w/ sp1_vrs & sp2_vrs) of lists and specify the modified density-independent vital rate
#' 
#' @return A wide dataframe with columns for time and population structure of the two species
#'
#' @export
sdcomp_project_dstb <- function(sdcomp_model, sp1_0, sp2_0, timestep, dstb_func = NULL, vrstoch_func = NULL) {
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
    new_sdcomp <- NULL
    
    # Update density-independent vital rates unless vrstoch is null
    if(is.null(vrstoch_func)) { new_sdcomp <- sdcomp_model }
    else {
      new_vrs <- vrstoch_func(t, sdcomp_model, sp1_curr, sp2_curr)
      new_sdcomp <- modify_vrs(sdcomp_model, new_vrs$sp1_vrs, new_vrs$sp2_vrs)
    }
    
    # Apply disturbance unless dstb_func is null
    if(!is.null(dstb_func)) {
      # Get absolute disturbance
      dstb <- dstb_func(t, sp1_curr, sp2_curr)
      sp1_dstb <- dstb$sp1
      sp2_dstb <- dstb$sp2
      
      # Apply disturbance (set minimum after disturbance as 0)
      sp1_curr <- pmax(sp1_curr + sp1_dstb, 0)
      sp2_curr <- pmax(sp2_curr + sp2_dstb, 0)
    }
    
    # Project by one time step
    comm_curr <- sdcomp_project_step(new_sdcomp, sp1_curr, sp2_curr)
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
