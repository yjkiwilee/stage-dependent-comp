################################################################################
#' Script for manipulating stage-dependence from base model
#'
#' Written by Young Jun Lee
#' Mar 2025

# Range of stage-dependence to impose
sd_manip_scale <- seq(-0.4, 0.4, 0.05)

# List for storing the resulting models
sd_manip_models <- list()

# Information about the focal vital rate whose stage-dependence is being manipulated
focal_vr_type <- "s"
focal_kmat_rows <- c(1, 1, 3, 3)
focal_kmat_cols <- c(1, 3, 1, 3)
manip_coeff <- c(1, 1, 1, 1) # Scaling factors for manipulating each kmat element

# Calculate model for each level of stage-dependence
# This will take a long time!
for(i in 1:length(sd_manip_scale)) {
  cat(paste("Calculating model for stage-dependence manipulation of", sd_manip_scale[i]))
  
  # Get focal kmat element values
  km_elems <- base_model$k_mat[[focal_vr_type]][cbind(focal_kmat_rows, focal_kmat_cols)]
  
  # Calculate new desired values for these elements
  km_alt_vals <- km_elems * (1 + manip_coeff * sd_manip_scale[i])
  
  # Calculate model with altered stage-dependence
  temp_model <- alter_k_mat_mult(base_model, vr_type = focal_vr_type,
                                 mat_rows = focal_kmat_rows,
                                 mat_cols = focal_kmat_cols,
                                 alt_vals = km_alt_vals,
                                 iter = 20, scal_fact_l = -100, scal_fact_r = 100)
  
  # Save model
  sd_manip_models[[i]] <- temp_model
  
  cat(" Done!\n")
}

# Validate average interaction strength for each manipulations
avg_int_strs <- lapply(sd_manip_models, function(manip_model) {
  calc_avg_int_str_mult(manip_model, "s", mat_rows = c(1,1,3,3), species = c(1,2,1,2))
})

do.call(rbind, avg_int_strs)
