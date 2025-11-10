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
    # Survival, progression, retrogression are binary
    is_logistic <- (vr_type %in% c("s", "g", "r"))
    
    comm_vrs <- matrix(c(sp1_vrs[[vr_type]], sp2_vrs[[vr_type]]), ncol = 1)
    
    intrcpt <- matrix()
    
    if(is_logistic) {
      intrcpt <- mod_logit(comm_vrs) # Inverse function of logistic
    } else {
      intrcpt <- mod_log(comm_vrs) # Inverse function of exponential
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
    is_logistic <- (vr_type %in% c("s", "g", "r"))
    
    intrcpt <- k_mat[[vr_type]][, 2 * n_stages + 1]
    comm_vrs <- c()
    
    if(is_logistic) {
      comm_vrs <- mod_logistic(intrcpt) # Logistic function of intercept
    } else {
      comm_vrs <- mod_exp(intrcpt) # Exponential function of intercept
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
#' With n stages, the vital rates of each species should include: s (survival probabilities, n),
#' g (maturation probabilities, n-1), r (retrogression probabilities, n-1), and f (fecundity, n-1).
#' With n stages, the competitive coefficient matrices is a list that contains s, g, r, and f.
#' Each contain the competitive coefficients controlling the vital rates.
#' Each are given as a 2m x 2n matrix, where m is the number of corresponding vital rates per species.
#' The final column in the survival and maturation matrices are the intercepts.
#' The vital rates are ordered by species and then by stage. Positive value denotes negative effect.
#'
#' @param n_stages Number of stages (n) in the model
#' @param sp1_vrs List containing the density-independent vital rates of species 1
#' @param sp2_vrs List containing the density-independent vital rates of species 2
#' @param k_mat Matrix containing competitive coefficients
#' @param sp1_vrs_var List containing the variances associated with vital rates of sp. 1
#' @param sp2_vrs_var List containing the variances associated with vital rates of sp. 2
#' 
#' @return The pre-breeding stage-dependent competition model, given as a list
#' 
#' @export
def_sdcomp_model <- function(n_stages, sp1_vrs, sp2_vrs, k_mat,
                             sp1_vr_sd,
                             sp2_vr_sd) {
  # Insert extra column to k_mats for the intercept
  k_mat$s <- cbind(k_mat$s, matrix(rep(0, nrow(k_mat$s)), ncol = 1))
  k_mat$g <- cbind(k_mat$g, matrix(rep(0, nrow(k_mat$g)), ncol = 1))
  k_mat$r <- cbind(k_mat$r, matrix(rep(0, nrow(k_mat$r)), ncol = 1))
  k_mat$f <- cbind(k_mat$f, matrix(rep(0, nrow(k_mat$f)), ncol = 1))
  
  # Update k_mat to reflect density-independent vital rates
  k_mat <- vrs_to_kmat(n_stages, sp1_vrs, sp2_vrs, k_mat)
  
  # Get number of vital rates
  sp1_n_vrs <- sum(sapply(sp1_vrs, function(i) { length(i) }))
  sp2_n_vrs <- sum(sapply(sp1_vrs, function(i) { length(i) }))
  
  # Store model
  sdcomp_model <- list(
    n_stages = n_stages,
    sp1_vrs = sp1_vrs,
    sp2_vrs = sp2_vrs,
    sp1_n_vrs = sp1_n_vrs,
    sp2_n_vrs = sp2_n_vrs,
    k_mat = k_mat,
    sp1_vr_sd = sp1_vr_sd,
    sp2_vr_sd = sp2_vr_sd
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
def_sdcomp_model_intrcpt <- function(n_stages, k_mat,
                                     sp1_vr_sd = list(s = 0, g = 0, r = 0, f = 0),
                                     sp2_vr_sd = list(s = 0, g = 0, r = 0, f = 0)) {
  # Convert k_mat intercepts into species density-independent vital rates
  vrs <- kmat_to_vrs(n_stages, k_mat)
  sp1_vrs <- vrs$sp1_vrs
  sp2_vrs <- vrs$sp2_vrs
  
  # Get number of vital rates
  sp1_n_vrs <- sum(sapply(sp1_vrs, function(i) { length(i) }))
  sp2_n_vrs <- sum(sapply(sp1_vrs, function(i) { length(i) }))
  
  # List to store the model in
  sdcomp_model <- list(
    n_stages = n_stages,
    sp1_vrs = sp1_vrs,
    sp2_vrs = sp2_vrs,
    sp1_n_vrs = sp1_n_vrs,
    sp2_n_vrs = sp2_n_vrs,
    k_mat = k_mat,
    sp1_vr_sd = sp1_vr_sd,
    sp2_vr_sd = sp2_vr_sd
  )
  
  # Return model
  sdcomp_model
}

#' Modify the vital rates of species in a model
#'
#' This function modifies the vital rates of species in a pre-defined sdcompr model.
#' 
#' @param sdcomp_model Model to be modified
#' @param sp1_vrs List containing the vital rates of species 1 (must have s, g, r, f components)
#' @param sp2_vrs List containing the vital rates of species 2 (must have s, g, r, f components)
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

#' Convert a matrix to string
mat_to_str <- function(mat) {
  r_collapsed <- apply(mat, 1, function(row) {
    paste(row, collapse = ",")
  })
  return(paste(r_collapsed, collapse = ";"))
}

#' Convert a string to matrix
str_to_mat <- function(mat_str) {
  rows <- strsplit(mat_str, ";")[[1]]
  mat_list <- lapply(rows, function(row) { as.numeric(strsplit(row, ",")[[1]]) })
  mat <- do.call(rbind, mat_list)
  rownames(mat) <- NULL
  return(mat)
}

#' Convert a sdcomp model into standardised string for export
sdcomp_to_str <- function(sdcomp_model, model_name = "") {
  paste0(
    "model\t", model_name, "\n",
    "\tn_stages\t", sdcomp_model$n_stages, "\n",
    do.call(paste0, lapply(names(sdcomp_model$sp1_vrs), function(vr_name) {
      paste0("\tsp1_vr_", vr_name, "\t", paste(sdcomp_model$sp1_vrs[[vr_name]], collapse = ","), "\n")
    })),
    do.call(paste0, lapply(names(sdcomp_model$sp2_vrs), function(vr_name) {
      paste0("\tsp2_vr_", vr_name, "\t", paste(sdcomp_model$sp1_vrs[[vr_name]], collapse = ","), "\n")
    })),
    "\tsp1_n_vrs\t", sdcomp_model$sp1_n_vrs, "\n",
    "\tsp2_n_vrs\t", sdcomp_model$sp2_n_vrs, "\n",
    do.call(paste0, lapply(names(sdcomp_model$k_mat), function(vr_name) {
      paste0("\tk_mat_", vr_name, "\t", mat_to_str(sdcomp_model$k_mat[[vr_name]]), "\n")
    })),
    do.call(paste0, lapply(names(sdcomp_model$sp1_vr_sd), function(vr_name) {
      paste0("\tsp1_vr_sd_", vr_name, "\t", paste(sdcomp_model$sp1_vr_sd[[vr_name]], collapse = ","), "\n")
    })),
    do.call(paste0, lapply(names(sdcomp_model$sp2_vr_sd), function(vr_name) {
      paste0("\tsp2_vr_sd_", vr_name, "\t", paste(sdcomp_model$sp2_vr_sd[[vr_name]], collapse = ","), "\n")
    }))
  )
}

#' Convert string representation of vectors to list
str_to_list <- function(vec_str) {
  lapply(strsplit(vec_str, ","), as.numeric)
}

#' Convert string to sdcomp model
str_to_sdcomp <- function(sdcomp_str) {
  # Define sdcomp model
  model_name <- str_match(sdcomp_str, "model\\t([a-zA-Z0-9._]*)\\n")[1,2]
  n_stages <- as.numeric(str_match(sdcomp_str, "\\tn_stages\\t([0-9]+)\\n")[1,2])
  sp1_vrs_match <- str_match_all(sdcomp_str, "\\tsp1_vr_([a-z])\\t(.+)\\n")[[1]]
  sp1_vrs <- setNames(str_to_list(sp1_vrs_match[,3]), sp1_vrs_match[,2])
  sp2_vrs_match <- str_match_all(sdcomp_str, "\\tsp2_vr_([a-z])\\t(.+)\\n")[[1]]
  sp2_vrs <- setNames(str_to_list(sp2_vrs_match[,3]), sp2_vrs_match[,2])
  k_mat_match <- str_match_all(sdcomp_str, "\\tk_mat_([a-z])\\t(.+)\\n")[[1]]
  k_mat <- setNames(lapply(k_mat_match[,3],
                           function(km_str) {
                             km <- str_to_mat(km_str)
                             km[,-ncol(km)]
                           }), k_mat_match[,2])
  sp1_vr_sd_match <- str_match_all(sdcomp_str, "\\tsp1_vr_sd_([a-z])\\t(.+)\\n")[[1]]
  sp1_vr_sd <- setNames(str_to_list(sp1_vr_sd_match[,3]), sp1_vr_sd_match[,2])
  sp2_vr_sd_match <- str_match_all(sdcomp_str, "\\tsp2_vr_sd_([a-z])\\t(.+)\\n")[[1]]
  sp2_vr_sd <- setNames(str_to_list(sp2_vr_sd_match[,3]), sp2_vr_sd_match[,2])
  
  def_sdcomp_model(
    n_stages = n_stages,
    sp1_vrs = sp1_vrs,
    sp2_vrs = sp2_vrs,
    k_mat = k_mat,
    sp1_vr_sd = sp1_vr_sd,
    sp2_vr_sd = sp2_vr_sd
  )
}

