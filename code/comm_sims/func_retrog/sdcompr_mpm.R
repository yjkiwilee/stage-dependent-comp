#' Calculate a matrix population model from vital rates
#'
#' This function calculates an n-stage matrix population model from vital rates.
#' The vital rates are given as a list, which should include:
#' s (survival probabilities, n),
#' g (maturation probabilities, n-1), r (retrogression probabilities, n-1),
#' f (fecundity, n-1), and p (recruitment factor; recruits per unit fecundity, n-1)
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
  rec <- vrs$p
  
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
    # Reproduction (with consideration of seed survival)
    matF[1, i] <- rec[i-1] * fec[i-1]
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
#' while fecundity & recruitment factor uses a negative exponential function.
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
  
  # Lists to store vital rates after density-dependent scaling
  sp1_dd_vrs <- list()
  sp2_dd_vrs <- list()
  
  # Competitive coefficients + intercepts
  k_mat <- sdcomp_model$k_mat
  
  # For each vital rate type, calculate the version after density-dependent scaling
  for(vr_type in names(k_mat)) {
    # Survival, maturation, and retrogression are binary
    is_logistic <- (vr_type %in% c("s", "g", "r"))
    
    # Convert sp1/2_pop to column vector
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

#' Get the density-independent MPM for a species from a sdcomp model
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
