#' Collapse the vital rates of a species into a named vector
#' 
#' This function collapses the vital rates of a species, given as a list,
#' into a named vector.
#'
#' @param vrs Vital rates
#' @param sp_prefix Optional species prefix
#'
#' @return A named vector of vital rate values where the names are the vital rate names
#' 
#' @export
collapse_vrs <- function(vrs, prefix = "", suffix = "") {
  # Get namd vectors for each vital rate
  vr_vecs <- lapply(names(vrs), function(vname) {
    vr_vals <- vrs[[vname]]
    names(vr_vals) <- paste0(prefix, vname, seq(1, length(vrs[[vname]])), suffix)
    
    return(vr_vals)
  })
  
  # Merge into single df
  return(do.call(c, vr_vecs))
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
#' @return A list of population structure of the two species
#'
#' @export
sdcomp_project_step <- function(sdcomp_model, sp1_pop, sp2_pop,
                                sp1_rnorm = NULL, sp2_rnorm = NULL) {
  n_stages <- sdcomp_model$n_stages
  
  # Get MPM for timestep
  comm_mpm <- sdcomp_to_mpm(sdcomp_model, sp1_pop, sp2_pop, sp1_rnorm, sp2_rnorm)
  
  # Merge population vectors into column vector
  comm_cv <- matrix(c(sp1_pop, sp2_pop), ncol = 1)
  
  # Project by one time step
  # Non-reproductive projection
  comm_cv_f_nr <- comm_mpm$matU %*% comm_cv
  # Reproductive projection
  # comm_cv_f_r <- apply(
  #   comm_mpm$matF %*% comm_cv,
  #   c(1, 2),
  #   function(x) { rpois(1, x) }
  # )
  comm_cv_f_r <- comm_mpm$matF %*% comm_cv
  comm_vec <- as.vector(comm_cv_f_nr + comm_cv_f_r)
  
  # Take floor of every element
  # comm_vec <- floor(comm_vec)
  
  # Split and return as list
  sp1_pop_f <- comm_vec[1:n_stages]
  sp2_pop_f <- comm_vec[(n_stages + 1):(2 * n_stages)]
  list(
    sp1_pop = sp1_pop_f,
    sp2_pop = sp2_pop_f
  )
}

#' Project the community one time step forward in time & return vital rates
#'
#' This function projects the community forward in time by one time step.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_vec Vector with the population structure of species 1
#' @param sp2_vec Vector with the population structure of species 2
#' @param min_size Minimum size of stage
#' 
#' @return A list containing population structure and the vital rates for the transition
#'
#' @export
sdcomp_project_step_vr <- function(sdcomp_model, sp1_pop, sp2_pop,
                                   sp1_rnorm = NULL, sp2_rnorm = NULL) {
  n_stages <- sdcomp_model$n_stages
  
  # Get MPM for timestep
  comm_mpm_vr <- sdcomp_to_mpm_vr(sdcomp_model, sp1_pop, sp2_pop, sp1_rnorm, sp2_rnorm)
  comm_mpm <- comm_mpm_vr$comm_mpm
  
  # Merge population vectors into column vector
  comm_cv <- matrix(c(sp1_pop, sp2_pop), ncol = 1)
  
  # Project by one time step
  # Non-reproductive projection
  comm_cv_f_nr <- comm_mpm$matU %*% comm_cv
  # Reproductive projection
  # comm_cv_f_r <- apply(
  #   comm_mpm$matF %*% comm_cv,
  #   c(1, 2),
  #   function(x) { rpois(1, x) }
  # )
  comm_cv_f_r <- comm_mpm$matF %*% comm_cv
  comm_vec <- as.vector(comm_cv_f_nr + comm_cv_f_r)
  
  # Take floor of every element
  # comm_vec <- floor(comm_vec)
  
  # Split and return as list
  sp1_pop_f <- comm_vec[1:n_stages]
  sp2_pop_f <- comm_vec[(n_stages + 1):(2 * n_stages)]
  list(
    sp1_pop = sp1_pop_f,
    sp2_pop = sp2_pop_f,
    sp1_vrs = comm_mpm_vr$sp1_vrs,
    sp2_vrs = comm_mpm_vr$sp2_vrs
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
sdcomp_project <- function(sdcomp_model, sp1_0, sp2_0, timestep,
                           sp1_rnorms = NULL, sp2_rnorms = NULL) {
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
    comm_curr <- sdcomp_project_step(sdcomp_model, sp1_curr, sp2_curr,
                                     if(is.null(sp1_rnorms)){NULL}else{sp1_rnorms[[t]]},
                                     if(is.null(sp2_rnorms)){NULL}else{sp2_rnorms[[t]]})
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

#' Project the community through time to get a time-series & also get vital rates
#'
#' This function projects the community forward in time and returns the trajectory.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_0 Vector with the initial population structure of species 1
#' @param sp2_0 Vector with the initial population structure of species 2
#' @param timestep Number of time steps to simulate
#' 
#' @return A wide dataframe with columns for time, population structure and vital rates of the two species
#'
#' @export
sdcomp_project_vr <- function(sdcomp_model, sp1_0, sp2_0, timestep,
                              sp1_rnorms = NULL, sp2_rnorms = NULL) {
  n_stages <- sdcomp_model$n_stages
  n_vrs <- sdcomp_model$sp1_n_vrs + sdcomp_model$sp2_n_vrs
  
  # Vector for storing time
  time_vec <- 0:timestep
  
  # Vector for storing current population structure
  sp1_curr <- sp1_0
  sp2_curr <- sp2_0
  
  # Matrix for storing the community structure through time
  comm_mat <- matrix(c(sp1_0, sp2_0), ncol = n_stages * 2)
  colnames(comm_mat) <- c(paste0("sp1_s", 1:n_stages), paste0("sp2_s", 1:n_stages))
  
  # Matrix for storing the vital rates through time
  vr_mat <- NULL
  
  # Simulate
  for(t in 1:timestep) {
    # Project by one time step
    comm_curr <- sdcomp_project_step_vr(sdcomp_model, sp1_curr, sp2_curr,
                                        if(is.null(sp1_rnorms)){NULL}else{sp1_rnorms[[t]]},
                                        if(is.null(sp2_rnorms)){NULL}else{sp2_rnorms[[t]]})
    
    sp1_curr <- comm_curr$sp1_pop
    sp2_curr <- comm_curr$sp2_pop
    
    # Store community structure
    comm_row <- c(
      comm_curr$sp1_pop,
      comm_curr$sp2_pop
    )
    comm_mat <- rbind(comm_mat, comm_row)
    
    # Extract species vital rates
    vr_row <- c(
      collapse_vrs(comm_curr$sp1_vrs, prefix = "vr_sp1_"),
      collapse_vrs(comm_curr$sp2_vrs, prefix = "vr_sp2_")
    )

    # Bind vital rates to matrix
    if(is.null(vr_mat)) {
      vr_mat <- matrix(vr_row, nrow = 1, dimnames = list(NULL, names(vr_row)))
    } else {
      vr_mat <- rbind(
        vr_mat,
        matrix(vr_row, nrow = 1, dimnames = list(NULL, names(vr_row)))
      )
    }
  }
  
  # Add an additional NA row to vital rates
  vr_mat <- rbind(
    vr_mat,
    matrix(rep(NA, times = n_vrs), nrow = 1)
  )
  
  # Convert time to matrix
  time_mat <- matrix(time_vec, ncol = 1)
  colnames(time_mat) <- c("time")
  
  # Bind matrices together
  res_mat <- cbind(
    time_mat, comm_mat, vr_mat
  )
  
  # Coerce into dataframe
  res_df <- data.frame(res_mat)
  rownames(res_df) <- NULL
  
  # Return
  res_df
}

#' Convert output of sdcomp_project into a longer format more suitable for plotting
#'
#' This function converts the output of sdcomp_project into a df more suitable for time series plotting.
#'
#' @param project_out
#' 
#' @return Processed df
#'
#' @export
long_comm_from_proj <- function(project_out) {
  out_comm <- tibble()
  
  # Make community structure longer
  out_comm <- project_out %>%
    select(time, starts_with("sp")) %>%
    distinct() %>%
    pivot_longer(
      cols = !time,
      names_to = c("species", "stage"),
      names_pattern = "^sp([0-9]+)_s([0-9]+)$",
      values_to = "size"
    )
  
  # Return df
  return(out_comm)
}

#' Extract vital rate data from sdcomp_project output in longer format
#'
#' This function extracts only vital rate data from sdcomp_project output.
#'
#' @param project_out
#' 
#' @return Processed df
#'
#' @export
long_vr_from_proj <- function(project_out) {
  out_vr <- project_out %>%
    select(time, starts_with("vr")) %>%
    pivot_longer(starts_with("vr"), names_to = "vr_name", values_to = "vr_val")
  
  # Process vr_name to get vital rate info
  vr_info <- str_match(out_vr$vr_name, "^vr_sp([0-9]+)_([a-z]+)([0-9]+)$")
  
  out_vr <- out_vr %>%
    mutate(
      vr_name = str_match(vr_name, "^vr_(.*)$")[,2],
      species = vr_info[,2],
      vr_type = vr_info[,3],
      vr_index = vr_info[,4]
    )
  
  return(out_vr)
}

#' Extract demographic variables from sdcomp_project output 
#'
#' This function extracts demographic variables from sdcomp_project output.
#'
#' @param project_out
#' 
#' @return DF
#'
#' @export
demo_from_proj <- function(project_out, n_stages) {
  # Get long vr
  long_vr <- long_vr_from_proj(project_out) %>%
    mutate(vr_index = as.numeric(vr_index))
  
  # Matrix for storing results
  res_mat <- matrix(nrow = 0, ncol = 6)
  colnames(res_mat) <- c("time", "species",
                         "lambda", "mean_life_expect",
                         "gen_time", "net_reprod_rate")
  
  # At each time and species, collate vital rates into a list
  for(t in unique(long_vr$time)) {
    for(sp in unique(long_vr$species)) {
      sp_vrs <- list()
      
      for(vr_type in unique(long_vr$vr_type)) {
        # Vector for subsetting data
        vr_subset <- long_vr$time == t &
          long_vr$species == sp &
          long_vr$vr_type == vr_type
        
        # Extract values
        vr_vals <- long_vr$vr_val[vr_subset][long_vr$vr_index[vr_subset]]
        
        # Add to species vr list
        sp_vrs[[vr_type]] <- vr_vals
      }
      
      # Convert vital rates to mpm
      sp_mpm <- vrs_to_mpm(n_stages, sp_vrs)
      
      # Is NA?
      is_mpm_na <- TRUE %in% is.na(sp_mpm$matA)
      
      if(is_mpm_na) {
        res_mat <- rbind(
          res_mat,
          matrix(c(
            t, as.numeric(sp),
            NA,
            NA,
            NA,
            NA
          ), ncol = 6)
        )
      } else {
        res_mat <- rbind(
          res_mat,
          matrix(c(
            t, as.numeric(sp),
            lambda(sp_mpm$matA),
            life_expect_mean(sp_mpm$matU),
            gen_time(sp_mpm$matU, sp_mpm$matF),
            net_repro_rate(sp_mpm$matU, sp_mpm$matF)
          ), ncol = 6)
        )
      }
    }
  }
  
  res_df <- as.data.frame(res_mat)
  res_df$species <- as.factor(res_df$species)
  
  return(res_df)
}
