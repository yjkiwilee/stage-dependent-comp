#' Calculate a matrix population model from vital rates
#'
#' This function calculates an n-stage matrix population model from vital rates.
#' The vital rates are given as a list, which should include:
#' s (survival probabilities, n),
#' g (maturation probabilities, n-1), r (retrogression probabilities, n-1),
#' f (fecundity, n-1)
#' 
#' NB: Fecundity should include factors such as seed / seedling survival.
#' NB This formulation assumes that progression & retrogression occurs just
#' before the survey.
#' 
#' This is a pre-breeding matrix.
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
  retrog <- vrs$r
  
  # Survival and fecundity matrices
  matU <- matrix(rep(0, n_stages^2), nrow = n_stages, ncol = n_stages)
  matF <- matrix(rep(0, n_stages^2), nrow = n_stages, ncol = n_stages)
  
  # ===== Calculate matU =====
  # For each starting stage
  for(i in 1:n_stages) {
    if(i == n_stages) { # If final stage
      # Stasis
      matU[i, i] <- surv[i] * (1 - retrog[i-1])
      # Retrogression
      matU[i-1, i] <- surv[i] * retrog[i-1]
    } else if(i == 1) { # If first stage
      # Stasis
      matU[i, i] <- surv[i] * (1 - matur[i])
      # Maturation
      matU[i+1, i] <- surv[i] * matur[i]
    } else { # Otherwise
      # Stasis
      matU[i, i] <- surv[i] * (1 - matur[i] - retrog[i-1])
      # Maturation
      matU[i+1, i] <- surv[i] * matur[i]
      # Retrogression
      matU[i-1, i] <- surv[i] * retrog[i-1]
    }
  }
  
  # ===== Calculate matF =====
  # For each reproductive stage
  for(i in 2:n_stages) {
    # Reproduction
    matF[1, i] <- fec[i-1]
  }
  
  # Create & return mpm
  mpm <- list(
    matU = matU,
    matF = matF,
    matA = matU + matF
  )
  
  mpm
}

vrs_to_mpm_fast <- function(n_stages, vrs) {
  
  # Extract vital rates
  surv <- vrs$s
  matur <- c(vrs$g, 0)
  fec <- c(0, vrs$f)
  retrog <- c(0, vrs$r)
  
  # ===== Calculate matA =====
  matA <- matrix(rep(0, n_stages^2), nrow = n_stages, ncol = n_stages)
  # For each starting stage
  for(i in 1:n_stages) {
    if(i == 1) { # If first stage
      # Stasis
      matA[i, i] <- surv[i] * (1 - matur[i])
      # Maturation
      matA[i+1, i] <- surv[i] * matur[i]
    } else if(i == n_stages) { # If final stage
      # Stasis
      matA[i, i] <- surv[i] * (1 - retrog[i-1])
      # Retrogression
      matA[i-1, i] <- surv[i] * retrog[i-1]
      # Reproduction
      matA[1, i] <- matA[1, i] + fec[i-1]
    } else { # Otherwise
      # Stasis
      matA[i, i] <- surv[i] * (1 - matur[i] - retrog[i-1])
      # Maturation
      matA[i+1, i] <- surv[i] * matur[i]
      # Retrogression
      matA[i-1, i] <- surv[i] * retrog[i-1]
      # Reproduction
      matA[1, i] <- matA[1, i] + fec[i-1]
    }
  }
  
  # Return matA
  return(matA)
}

#' Calculate density-dependent vital rates from population structure
#'
#' This function calculates the density-dependent vital rates
#' given a community model and the population structure.
#' This function uses the vital rate variances defined at model definition to simulate stochasticity.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_vec Vector with the population structure of species 1
#' @param sp2_vec Vector with the population structure of species 2
#' 
#' @return Two species' vital rates, given as a list
#'
#' @export
sdcomp_to_vrs <- function(sdcomp_model, sp1_pop, sp2_pop, sp1_rnorm = NULL, sp2_rnorm = NULL) {
  
  # Number of stages
  n_stages <- sdcomp_model$n_stages
  
  # Lists to store vital rates after density-dependent scaling
  sp1_dd_vrs <- list()
  sp2_dd_vrs <- list()
  
  # Competitive coefficients + intercepts
  k_mat <- sdcomp_model$k_mat
  
  # For each vital rate type, calculate the version after density-dependent scaling
  # and stochastic vital rates
  for(vr_type in names(k_mat)) {
    # Number of vital rates
    n_vrs_sp1 <- length(sdcomp_model$sp1_vrs[[vr_type]])
    n_vrs_sp2 <- length(sdcomp_model$sp2_vrs[[vr_type]])
    
    # Survival, maturation, and retrogression are binary
    is_logistic <- (vr_type %in% c("s", "g", "r"))
    
    # Convert sp1/2_pop to column vector
    comm_cv <- matrix(c(sp1_pop, sp2_pop, 1), ncol = 1)
    
    # Variable to store scaling factor for vital rates
    dd_vrs <- c()
    
    # Linear sum to plug into transformation function
    dd_sums <- k_mat[[vr_type]] %*% comm_cv
    
    # Apply stochastic vital rates
    if(is.null(sp1_rnorm)) {
      dd_sums[1:(length(dd_sums) / 2)] <- rnorm(n_vrs_sp1) *
        sdcomp_model$sp1_vr_sd[[vr_type]] +
        dd_sums[1:(length(dd_sums) / 2)]
    } else {
      dd_sums[1:(length(dd_sums) / 2)] <- sp1_rnorm[[vr_type]] *
        sdcomp_model$sp1_vr_sd[[vr_type]] +
        dd_sums[1:(length(dd_sums) / 2)]
    }
    
    if(is.null(sp2_rnorm)) {
      dd_sums[(length(dd_sums) / 2 + 1):length(dd_sums)] <- rnorm(n_vrs_sp2) *
        sdcomp_model$sp2_vr_sd[[vr_type]] +
        dd_sums[(length(dd_sums) / 2 + 1):length(dd_sums)]
    } else {
      dd_sums[(length(dd_sums) / 2 + 1):length(dd_sums)] <- sp2_rnorm[[vr_type]] *
        sdcomp_model$sp2_vr_sd[[vr_type]] +
        dd_sums[(length(dd_sums) / 2 + 1):length(dd_sums)]
    }
    
    # Calculate final density-dependent vital rates
    if(is_logistic) { # If survival or maturation prob.
      dd_vrs <- as.vector(mod_logistic(dd_sums)) # Bound linear scaling
    } else { # If fecundity
      dd_vrs <- as.vector(mod_exp(dd_sums)) # Identity scaling
    }
    
    # Store density-dependent vital rates
    sp1_dd_vrs[[vr_type]] <- dd_vrs[1:(length(dd_vrs) / 2)]
    sp2_dd_vrs[[vr_type]] <- dd_vrs[(length(dd_vrs) / 2 + 1):length(dd_vrs)]
  }
  
  # Return density-dependent vital rates
  return(list(sp1_vrs = sp1_dd_vrs, sp2_vrs = sp2_dd_vrs))
}

#' Calculate a two-species projection matrix from population structure
#'
#' This function calculates a two-species projection matrix from a stage-dependent competition model,
#' given some population structure. Probabilities use a logistic density-dependent function,
#' while fecundity & recruitment factor uses a negative exponential function.
#' This function uses the vital rate variances defined at model definition to simulate stochasticity.
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_vec Vector with the population structure of species 1
#' @param sp2_vec Vector with the population structure of species 2
#' 
#' @return An Rage two-species matrix, given as a list
#'
#' @export
sdcomp_to_mpm <- function(sdcomp_model, sp1_pop, sp2_pop, sp1_rnorm = NULL, sp2_rnorm = NULL) {
  # Get density-dependent vital rates
  dd_vrs <- sdcomp_to_vrs(sdcomp_model, sp1_pop, sp2_pop, sp1_rnorm, sp2_rnorm)
  
  # Construct matrix from vital rates while passing variance in vital rates
  sp1_mpm <- vrs_to_mpm(n_stages, dd_vrs$sp1_vrs)
  sp2_mpm <- vrs_to_mpm(n_stages, dd_vrs$sp2_vrs)
  
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

#' Calculate a two-species projection matrix from population structure, but also return vital rates
#'
#' This function is identical to sdcomp_to_mpm but also returns the density-dependent vital rates
#'
#' @param sdcomp_model Stage-dependent competition model
#' @param sp1_vec Vector with the population structure of species 1
#' @param sp2_vec Vector with the population structure of species 2
#' 
#' @return A list containing the two species' density-dependent vital rates & Rage MPM
#'
#' @export
sdcomp_to_mpm_vr <- function(sdcomp_model, sp1_pop, sp2_pop, sp1_rnorm = NULL, sp2_rnorm = NULL) {
  # Get density-dependent vital rates
  dd_vrs <- sdcomp_to_vrs(sdcomp_model, sp1_pop, sp2_pop, sp1_rnorm, sp2_rnorm)
  
  # Construct matrix from vital rates while passing variance in vital rates
  sp1_mpm <- vrs_to_mpm(n_stages, dd_vrs$sp1_vrs)
  sp2_mpm <- vrs_to_mpm(n_stages, dd_vrs$sp2_vrs)
  
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
  
  # Return vital rates & MPM
  return(list(
    sp1_vrs = dd_vrs$sp1_vrs,
    sp2_vrs = dd_vrs$sp2_vrs,
    comm_mpm = comm_mpm
  ))
}

#' Get the density-independent deterministic MPM for a species from a sdcomp model
#'
#' This function takes in a sdcomp model and outputs an MPM for either one of the species
#' in a format that is compatible with Rage.
#'
#' @param sdcomp_model The sdcomp model, a list containing the vital rates of species 1 & 2 as well as competitive coefficients
#' @param species The species to extract the MPM of
#' 
#' @return A Rage MPM
#' 
#' @export
get_sp_mpm <- function(sdcomp_model, species) {
  # Get vital rates of the specified species
  vrs <- if(species == 1) {sdcomp_model$sp1_vrs}
    else {sdcomp_model$sp2_vrs}
  
  # Calculate & return MPM
  mpm <- vrs_to_mpm(sdcomp_model$n_stages, vrs)
  mpm
}
