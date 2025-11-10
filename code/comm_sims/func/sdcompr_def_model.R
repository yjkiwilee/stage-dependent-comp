#' Insert density-independent vital rates into k_mat intercepts
#'
#' This function converts the density-independent vital rates of two species into a k_mat.
#' 
#' @param n_stages Number of stages (n) in the model
#' @param sp1_vrs List containing the density-independent vital rates of species 1
#' @param sp2_vrs List containing the density-independent vital rates of species 2
#' @param prev_k_mat Previous k_mat
#' 
#' @return Updated k_mat
#' 
#' @export
vrs_to_kmat <- function(n_stages, sp1_vrs, sp2_vrs, prev_k_mat) {
  # List to store the new k_mat
  k_mat <- prev_k_mat
  
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
    
    # Insert intercepts into the last column of k_mat
    k_mat[[vr_type]][,ncol(k_mat[[vr_type]])] <- intrcpt
  }
  
  # Return k_mat
  k_mat
}

#' Calculate density-independent vital rates from k_mat
#'
#' This function calculates the density-independent vital rates from k_mat intercepts.
#' 
#' @param n_stages Number of stages (n) in the model
#' @param k_mat k_mat, with intercepts (last column) for density-independent vital rates of the two species
#' 
#' @return List containing vital rates of sp1 (sp1_vrs) and sp2 (sp2_vrs)
#' 
#' @export
kmat_to_vrs <- function(n_stages, k_mat) {
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
  
  # Return list
  list(sp1_vrs = sp1_vrs, sp2_vrs = sp2_vrs)
}

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
#' 
#' @return The pre-breeding stage-dependent competition model, given as a list
#' 
#' @export
def_sdcomp_model <- function(n_stages, sp1_vrs, sp2_vrs, k_mat) {
  # Insert extra column to k_mats
  k_mat$s <- cbind(k_mat$s, matrix(rep(0, nrow(k_mat$s)), ncol = 1))
  k_mat$g <- cbind(k_mat$g, matrix(rep(0, nrow(k_mat$g)), ncol = 1))
  k_mat$f <- cbind(k_mat$f, matrix(rep(0, nrow(k_mat$f)), ncol = 1))
  
  # Update k_mat to reflect density-independent vital rates
  k_mat <- vrs_to_kmat(n_stages, sp1_vrs, sp2_vrs, k_mat)
  
  # Store model
  sdcomp_model <- list(
    n_stages = n_stages,
    sp1_vrs = sp1_vrs,
    sp2_vrs = sp2_vrs,
    k_mat = k_mat
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
  # Convert k_mat intercepts into species density-independent vital rates
  vrs <- kmat_to_vrs(n_stages, k_mat)
  sp1_vrs <- vrs$sp1_vrs
  sp2_vrs <- vrs$sp2_vrs
  
  # List to store the model in
  sdcomp_model <- list(
    n_stages = n_stages,
    sp1_vrs = sp1_vrs,
    sp2_vrs = sp2_vrs,
    k_mat = k_mat
  )
  
  # Return model
  sdcomp_model
}

#' Modify the vital rates of species in a model
#'
#' This function modifies the vital rates of species in a pre-defined sdcompr model.
#' 
#' @param sdcomp_model Model to be modified
#' @param sp1_vrs List containing the vital rates of species 1 (must have s, g, and f components)
#' @param sp2_vrs List containing the vital rates of species 2 (must have s, g, and f components)
#' 
#' @return The modified sdcompr model, given as a list
#' 
#' @export
modify_vrs <- function(sdcomp_model, sp1_vrs = NULL, sp2_vrs = NULL) {
  new_sp1_vrs <- if(is.null(sp1_vrs)) sdcomp_model$sp1_vrs else sp1_vrs
  new_sp2_vrs <- if(is.null(sp2_vrs)) sdcomp_model$sp2_vrs else sp2_vrs
  
  # Convert density-independent vital rates into k_mat elements
  k_mat_temp <- NULL
  k_mat_temp <- vrs_to_kmat(sdcomp_model$n_stages, new_sp1_vrs, new_sp2_vrs, sdcomp_model$k_mat)
  
  # Update species vital rates
  sdcomp_model$sp1_vrs <- new_sp1_vrs
  sdcomp_model$sp2_vrs <- new_sp2_vrs
  
  # Update k_mat
  sdcomp_model$k_mat <- k_mat_temp
  
  # Return model
  sdcomp_model
}