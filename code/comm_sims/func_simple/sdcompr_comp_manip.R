#' Calculate average per-capita interaction strengths of a species given a sdcomp model (INCOMPATIBLE WITH CYCLES, STOCHASTIC AND CHAOS)
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
calc_avg_int_str <- function(sdcomp_model, vr_type, mat_row, species,
                             timestep = SIMEQ_TIMESTEP, sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT) {
  n_stages <- sdcomp_model$n_stages
  
  # Find 'equilibrium' of species with current parameters
  end_state <- sdcomp_project_end(sdcomp_model, sp1_0, sp2_0, timestep)
  
  end_state_sp <- end_state[[ifelse(species == 1, "sp1_pop", "sp2_pop")]] # End state for the specific species
  
   if(TRUE %in% is.na(end_state_sp)) { # Escape if the simulation 'blows up'
    warning("Population size exceeded calculable limit")
    return(NA)
   } else if(sum(end_state_sp) <= 0) { # Escape if the focal species goes extinct
     warning("Species extinction occurred")
     return(NA)
   }
  
  # Normalise end state
  end_state_sp_norm <- end_state_sp / sum(end_state_sp)
  
  # Determine the current per-capia interaction strength of the species on the specific vital rate
  mat_cols <- c((1:n_stages) + ((species - 1) * n_stages))
  k_mat_elems <- sdcomp_model$k_mat[[vr_type]][mat_row, mat_cols] # Extract k_mat elements for the species & vital rate
  avg_int_str <- sum(end_state_sp_norm * k_mat_elems) # Average per-capita interaction strength
  
  return(avg_int_str)
}

#' Calculate average per-capita interaction strengths of a species given a sdcomp model that's robust for cycles, stochastic and chaos
#' 
#' This function calculates the average per-capita competitive strength of a species on a specific vital rate at equilibrium.
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param vr_type Type of vital rate being affected
#' @param mat_row Row of matrix corresponding to vital rate being affected
#' @param species Species exerting the interactive effect. 1 = sp1, 2 = sp2
#' @param timestep Number of time steps to simulate to reach equilibrium (Default = SIMEQ_TIMESTEP)
#' @param sp1_0 Initial population structure for species 1 (Default = DEFAULT_POP_STRUCT)
#' @param sp2_0 Initial population structure for species 2 (Default = DEFAULT_POP_STRUCT)
#' @param avgtf Timeframe over which to average population structure (Default = SIMEQ_AVGTF)
#' 
#' @return The average per-capita competitive strength or NA if competitive exclusion occurs
#' 
#' @export
calc_avg_int_str_rb <- function(sdcomp_model, vr_type, mat_row, species, timestep = SIMEQ_TIMESTEP_ROBUST,
                                sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT, avgtf = SIMEQ_AVGTF) {
  n_stages <- sdcomp_model$n_stages
  
  # Find 'equilibrium' of species with current parameters
  end_state <- sdcomp_project_end_avg(sdcomp_model, sp1_0, sp2_0, timestep, avgtf)
  
  end_state_sp <- end_state[[ifelse(species == 1, "sp1_pop", "sp2_pop")]] # End state for the specific species
  
  if(TRUE %in% is.na(end_state_sp)) { # Escape if the simulation 'blows up'
    warning("Population size exceeded calculable limit")
    return(NA)
  } else if(sum(end_state_sp) <= 0) { # Escape if the focal species goes extinct
    warning("Species extinction occurred")
    return(NA)
  } 
  
  # Normalise end state
  end_state_sp_norm <- end_state_sp / sum(end_state_sp)
  
  # Determine the current per-capia interaction strength of the species on the specific vital rate
  mat_cols <- c((1:n_stages) + ((species - 1) * n_stages))
  k_mat_elems <- sdcomp_model$k_mat[[vr_type]][mat_row, mat_cols] # Extract k_mat elements for the species & vital rate
  avg_int_str <- sum(end_state_sp_norm * k_mat_elems) # Average per-capita interaction strength
  
  return(avg_int_str)
}

#' Calculate multiple average per-capita interaction strengths of a species given a sdcomp model that's robust for cycles, stochastic and chaos
#' 
#' This function calculates the average per-capita competitive strength of a species on a specific vital rate at equilibrium.
#' NOT COMPATIBLE WITH STOCHASTIC MODEL
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param vr_type Type of vital rate being affected
#' @param mat_rows Rows of matrix corresponding to vital rate being affected
#' @param species Species exerting the interactive effect. 1 = sp1, 2 = sp2 (NB: this is a vector)
#' @param timestep Number of time steps to simulate to reach equilibrium (Default = SIMEQ_TIMESTEP)
#' @param sp1_0 Initial population structure for species 1 (Default = DEFAULT_POP_STRUCT)
#' @param sp2_0 Initial population structure for species 2 (Default = DEFAULT_POP_STRUCT)
#' @param avgtf Timeframe over which to average population structure (Default = SIMEQ_AVGTF)
#' 
#' @return The average per-capita competitive strength or NA if competitive exclusion occurs
#' 
#' @export
calc_avg_int_str_mult <- function(sdcomp_model, vr_type, mat_rows, species, timestep = SIMEQ_TIMESTEP_ROBUST,
                                  sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT, avgtf = SIMEQ_AVGTF) {
  # Number of species-level average interactions strengths to compute
  n_eval <- length(mat_rows)
  # Number of stages
  n_stages <- sdcomp_model$n_stages
  
  # Find 'equilibrium' of species with current parameters
  end_state <- sdcomp_project_end_avg(sdcomp_model, sp1_0, sp2_0, timestep, avgtf)
  
  end_state_sp <- lapply(species, function(sp) { # Store end states relevant for each interaction strength being computed
    end_state[[ifelse(sp == 1, "sp1_pop", "sp2_pop")]] # End state for the specific species
  })
  
  if(TRUE %in% is.na(end_state_sp)) { # Escape if the simulation 'blows up'
    warning("Population size exceeded calculable limit")
    return(NA)
  } else if(TRUE %in% (sapply(end_state_sp, sum) <= 0)) { # Escape if the focal species goes extinct
    warning("Species extinction occurred")
    return(NA)
  } 
  
  # Normalise end state
  end_state_sp_norm <- lapply(end_state_sp, function(end_st) {
    end_st / sum(end_st)
  })
  
  # Determine the current per-capita interaction strength of the species on the specific vital rate
  avg_int_strs <- sapply(1:n_eval, function(i){
    mat_cols <- c((1:n_stages) + ((species[i] - 1) * n_stages))
    k_mat_elems <- sdcomp_model$k_mat[[vr_type]][mat_rows[i], mat_cols] # Extract k_mat elements for the species & vital rate
    avg_int_str <- sum(end_state_sp_norm[[i]] * k_mat_elems) # Average per-capita interaction strength
    # Return average per-capita interaction strength
    return(avg_int_str)
  })
  
  return(avg_int_strs)
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
alter_k_mat <- function(sdcomp_model, vr_type, mat_row, mat_col, alt_val, timestep = SIMEQ_TIMESTEP,
                        avgtf = SIMEQ_AVGTF, iter = 40, scal_fact_l = 0, scal_fact_r = 5,
                        sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT) {
  n_stages <- sdcomp_model$n_stages
  species <- (mat_col - 1) %/% n_stages + 1 # 1 if species 1, 2 if species 2
  
  # Determine the target per-capita interaction strength
  targ_int_str <- calc_avg_int_str(sdcomp_model, vr_type, mat_row, species, timestep = timestep, sp1_0 = sp1_0, sp2_0 = sp2_0)
  
  # Escape if sim returns NA
  if(is.na(targ_int_str)) {
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

#' Adjust a competitive coefficient while keeping species-level per-capita interaction strength equal (ROBUST TO STOCHASTIC, CYCLES AND CHAOS)
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
#' @param timestep Number of time steps to simulate to reach equilibrium (Default = SIMEQ_TIMESTEP_ROBUST)
#' @param avgtf Number of time steps to average over (Default = SIMEQ_AVGTF)
#' @param iter Number of iterations used to approsimate the solution
#' @param scal_fact_l Initial minimum scale factor to try for bisection
#' @param scal_fact_r Initial minimum scale factor to try for bisection
#' @param sp1_0 Initial population structure of species 1
#' @param sp2_0 Initial population structure of species 2
#' 
#' @return A sdcomp model with the adjustment controlling for average per-capita interaction strength
#' 
#' @export
alter_k_mat_rb <- function(sdcomp_model, vr_type, mat_row, mat_col, alt_val,
                           timestep = SIMEQ_TIMESTEP_ROBUST, avgtf = SIMEQ_AVGTF,
                           iter = 40, scal_fact_l = 0, scal_fact_r = 5,
                           sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT) {
  n_stages <- sdcomp_model$n_stages
  species <- (mat_col - 1) %/% n_stages + 1 # 1 if species 1, 2 if species 2
  
  # Determine the target per-capita interaction strength
  targ_int_str <- calc_avg_int_str_rb(sdcomp_model, vr_type, mat_row, species,
                                      timestep = timestep, avgtf = avgtf,
                                      sp1_0 = sp1_0, sp2_0 = sp2_0)
  
  # Escape if sim returns NA
  if(is.na(targ_int_str)) {
    return(NA)
  }
  
  # Identify the columns of elements to be adjusted to accommodate the altered value
  adj_cols <- (1:n_stages) + ((species - 1) * n_stages)
  adj_cols <- adj_cols[adj_cols != mat_col]
  
  # Variable to store the appropriate scaling factor for the elements to be adjusted
  scal_fact_m <- mean(c(scal_fact_l, scal_fact_r))
  # Variables to store previous scaling factor extremeties
  scal_fact_l_prev <- scal_fact_l
  scal_fact_r_prev <- scal_fact_r
  
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
    
    l_int_str <- calc_avg_int_str_rb(l_model, vr_type, mat_row, species, timestep, avgtf, sp1_0, sp2_0)
    r_int_str <- calc_avg_int_str_rb(r_model, vr_type, mat_row, species, timestep, avgtf, sp1_0, sp2_0)
    m_int_str <- calc_avg_int_str_rb(m_model, vr_type, mat_row, species, timestep, avgtf, sp1_0, sp2_0)
    
    if(l_int_str > targ_int_str | r_int_str < targ_int_str) { # If target is not within the extremeties
      warning("Failed to converge, backing up...") # NB: the function may still get stuck regardless
      scal_fact_l <- scal_fact_l_prev
      scal_fact_r <- scal_fact_r_prev
      next
    } else if(m_int_str < targ_int_str) {
      scal_fact_l_prev <- scal_fact_l
      scal_fact_l <- scal_fact_m
    } else {
      scal_fact_r_prev <- scal_fact_r
      scal_fact_r <- scal_fact_m
    }
    
    scal_fact_m <- mean(c(scal_fact_l, scal_fact_r))
  }
  
  # Return resulting model
  return(m_model)
}

#' Adjust one or more competitive coefficients while keeping species-level per-capita interaction strength equal (ROBUST TO STOCHASTIC, CYCLES AND CHAOS)
#' 
#' This function takes in a sdcomp model, and competitive coefficients to be adjusted.
#' It then returns a new sdcomp model where the per-capita interaction strength at the stable population structure is equal.
#' This function uses the bisection method to approximate the appropriate competitive coefficient values.
#' This function is a multi-parameter extension of alter_k_mat_rb().
#' This function does not support altering multiple coefficients related to multiple stages of the same species.
#' 
#' @param sdcomp_model Stage-dependent competition model
#' @param vr_type Type of vital rate being affected
#' @param mat_rows Rows of matrix elements to alter
#' @param mat_cols Columns of matrix elements to alter
#' @param alt_val Alternate values of the competitive coefficients
#' @param timestep Number of time steps to simulate to reach equilibrium (Default = SIMEQ_TIMESTEP_ROBUST)
#' @param avgtf Number of time steps to average over (Default = SIMEQ_AVGTF)
#' @param iter Number of iterations used to approsimate the solution
#' @param scal_fact_l Initial minimum scale factor to try for bisection, vector with equal length as mat_rows / mat_cols or single value
#' @param scal_fact_r Initial minimum scale factor to try for bisection, vector with equal length as mat_rows / mat_cols or single value
#' @param sp1_0 Initial population structure of species 1
#' @param sp2_0 Initial population structure of species 2
#' 
#' @return A sdcomp model with the adjustment controlling for average per-capita interaction strength
#' 
#' @export
alter_k_mat_mult <- function(sdcomp_model, vr_type, mat_rows, mat_cols, alt_vals,
                             timestep = SIMEQ_TIMESTEP_ROBUST, avgtf = SIMEQ_AVGTF,
                             iter = 40, scal_fact_l = 0, scal_fact_r = 5,
                             sp1_0 = DEFAULT_POP_STRUCT, sp2_0 = DEFAULT_POP_STRUCT) {
  # Number of coefficients to alter
  n_coeff <- length(mat_rows)
  # Number of stages in the model
  n_stages <- sdcomp_model$n_stages
  species <- (mat_cols - 1) %/% n_stages + 1 # 1 if species 1, 2 if species 2. NB: this variable is a vector.
  
  # Determine the target per-capita interaction strengths
  targ_int_strs <- calc_avg_int_str_mult(sdcomp_model, vr_type, mat_rows, species,
                                         timestep = timestep, avgtf = avgtf,
                                         sp1_0 = sp1_0, sp2_0 = sp2_0)
  
  # Escape if any sims return NA
  if(TRUE %in% is.na(targ_int_strs)) {
    return(NA)
  }
  
  # Identify the columns of elements to be adjusted to accommodate the altered value
  # List to store the columns to be adjusted
  adj_cols <- list()
  # For each coefficient to be altered, store columns to be adjusted
  for(i in 1:n_coeff) {
    adj_cols[[i]] <- (1:n_stages) + ((species[i] - 1) * n_stages)
    adj_cols[[i]] <- adj_cols[[i]][adj_cols[[i]] != mat_cols[i]] # Exclude column for the element being directly manipulated
  }
  
  # Throw error if the length of the scale factors are inconsistent with the number of coefficients being adjusted
  if(length(scal_fact_l) != 1 & length(scal_fact_l) != n_coeff) {
    stop("scal_fact_l must be of length 1 or have the same length as the number of coefficients being altered")
  } else if(length(scal_fact_r) != 1 & length(scal_fact_r) != n_coeff) {
    stop("scal_fact_r must be of length 1 or have the same length as the number of coefficients being altered")
  }
  # Create scal_fact_ls, scal_fact_rs to store left and right bounds for scaling factors corresponding to each of the coefficient groups being adjusted
  scal_fact_ls <- rep(scal_fact_l, ifelse(length(scal_fact_l) == 1, n_coeff, 1))
  scal_fact_rs <- rep(scal_fact_r, ifelse(length(scal_fact_r) == 1, n_coeff, 1))
  
  # Variable to store the appropriate scaling factor for the elements to be adjusted
  scal_fact_ms <- (scal_fact_ls + scal_fact_rs) / 2
  # Variables to store previous scaling factor extremeties
  scal_fact_ls_prev <- NA
  scal_fact_rs_prev <- NA
  
  # Alter k_mat elements
  sdcomp_model$k_mat[[vr_type]][cbind(mat_rows, mat_cols)] <- alt_vals
  
  # Variables to store sdcomp for left and right bounds
  l_model <- NA
  m_model <- NA
  r_model <- NA
  
  # Variables to store previous interaction strengths
  l_int_strs <- NA
  r_int_strs <- NA
  m_int_strs <- NA
  
  # Bisection method to find root
  for(i in 1:iter) {
    l_model <- sdcomp_model
    m_model <- sdcomp_model
    r_model <- sdcomp_model
    
    # Update models with scale factors
    for(j in 1:n_coeff) {
      mat_row <- mat_rows[j]
      adj_c <- adj_cols[[j]]
      l_model$k_mat[[vr_type]][mat_row, adj_c] <- l_model$k_mat[[vr_type]][mat_row, adj_c] * scal_fact_ls[j]
      r_model$k_mat[[vr_type]][mat_row, adj_c] <- r_model$k_mat[[vr_type]][mat_row, adj_c] * scal_fact_rs[j]
      m_model$k_mat[[vr_type]][mat_row, adj_c] <- m_model$k_mat[[vr_type]][mat_row, adj_c] * scal_fact_ms[j]
    }
    
    # Calculate average interaction strengths
    # Check if model simulation has been run before to avoid unnecessary calculations
    if((FALSE %in% (scal_fact_ls_prev == scal_fact_ls)) | (TRUE %in% is.na(scal_fact_ls_prev))) { # Calculate l_int_strs only if scal_fact_ls has been updated or hasn't been calculated before
      l_int_strs <- calc_avg_int_str_mult(l_model, vr_type, mat_rows, species,
                                          timestep = timestep, avgtf = avgtf,
                                          sp1_0 = sp1_0, sp2_0 = sp2_0)
    }
    if((FALSE %in% (scal_fact_rs_prev == scal_fact_rs)) | (TRUE %in% is.na(scal_fact_ls_prev))) { # Calculate r_int_strs only if scal_fact_rs has beeen updated or hasn't been calculated before
      r_int_strs <- calc_avg_int_str_mult(r_model, vr_type, mat_rows, species,
                                          timestep = timestep, avgtf = avgtf,
                                          sp1_0 = sp1_0, sp2_0 = sp2_0)
    }
    m_int_strs <- calc_avg_int_str_mult(m_model, vr_type, mat_rows, species,
                                        timestep = timestep, avgtf = avgtf,
                                        sp1_0 = sp1_0, sp2_0 = sp2_0)
    
    # Restrict extremeties
    for(j in 1:n_coeff) {
      if(l_int_strs[j] > targ_int_strs[j] | r_int_strs[j] < targ_int_strs[j]) { # If target is not within the extremeties, revert to the last value that did converge
        warning("Failed to converge, trying again...") # NB: the function may still get stuck regardless
        # scal_fact_ls[j] <- ifelse(is.na(scal_fact_ls_prev[j]), scal_fact_ls[j], scal_fact_ls_prev[j]) # scal_fact_ls_prev[j] will be NA in first iteration, in which case scal_fact_ls is used instead
        # scal_fact_rs[j] <- ifelse(is.na(scal_fact_rs_prev[j]), scal_fact_rs[j], scal_fact_rs_prev[j]) # Same as above
        next
      } else if(m_int_strs[j] < targ_int_strs[j]) {
        scal_fact_ls_prev[j] <- scal_fact_ls[j]
        scal_fact_ls[j] <- scal_fact_ms[j]
        scal_fact_rs_prev[j] <- scal_fact_rs[j]
      } else {
        scal_fact_rs_prev[j] <- scal_fact_rs[j]
        scal_fact_rs[j] <- scal_fact_ms[j]
        scal_fact_ls_prev[j] <- scal_fact_ls[j]
      }
    }
    
    # Update mean scale factors
    scal_fact_ms <- (scal_fact_ls + scal_fact_rs) / 2
  }
  
  # Return resulting model
  return(m_model)
}

# Calculate species-specific weights required to achieve a given equilibrium
# in deterministic simulation
find_sp_weights <- function(sdcomp_model, eq_size) {
  
}

# Adjust deterministic species-specific weights numerically
adjust_sp_weights <- function(sdcomp_model, sp_weights) {
  
}




