# Code for fitting model to ddvr and getting intercepts & slopes

# Timestep to start sampling
t_start <- 600
# Timestep to end sampling
t_end <- 2000
# Maximum timestep in simulations
t_max <- 2000

# Iterate through density-dependent vital rates, species combinations and stage-dependence
# vr_names <- c("sj", "sa", "g", "r", "f")
# sp_names <- c("F2", "F1", "M", "S1", "S2")
# sd_levels <- seq(-1, 1, 0.2)
# sd_levels_str <- as.character(sd_levels)
rand_seeds <- c(1:10, -1 * (1:10))

# DF to store model fit summaries
fit_mod_df <- NULL

# Convert vital rates to df
spp_vrs_df <- lapply(names(spp_vrs), function(ddvr) {
  lapply(names(spp_vrs[[ddvr]]), function(sp_name) {
    sp_vr <- spp_vrs[[ddvr]][[sp_name]]
    
    spp_vrs_row <- tibble(
      dd_vr = ddvr,
      sp_name = sp_name,
      sj = sp_vr$s[1],
      sa = sp_vr$s[2],
      g = sp_vr$g[1],
      r = sp_vr$r[1],
      f = sp_vr$f[1]
    )
    
    return(spp_vrs_row)
  })
})
spp_vrs_df <- bind_rows(spp_vrs_df)


for(ddvr in vr_names) {
  cat(paste("Density-dependent vital rate:", ddvr, "\n"))
  cat(paste0(format(Sys.time(), "%c"), "\n"))
  
  spp_vrs
  
  cat("\tLoading data... ")
  
  complete_tseries <- read_csv(
    file.path("result_data",
         paste0("sim_tseries_", ddvr, ".csv")),
    col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
  )
  
  complete_tseries <- complete_tseries %>%
    filter(time >= t_start & time < t_end)
  
  cat("done!\n\tInserting vital rates... ")
  
  # Calculate ddvr from stage-specific abundances and stochastic term
  if(ddvr == "f") {
    complete_tseries <- complete_tseries %>%
      left_join(
        kmat_df,
        by = join_by(dd_vr, sp1_name, sp2_name, st_dep_level)
      ) %>%
      mutate(
        sp1_dd_vr = mod_exp(
          sp1_cj * sp1_s1 + sp1_ca * sp1_s2 + sp1_hj * sp2_s1 + sp1_ha * sp2_s2 +
            sp1_intercept + sp1_stoch_ddvr
        ),
        sp2_dd_vr = mod_exp(
          sp2_cj * sp2_s1 + sp2_ca * sp2_s2 + sp2_hj * sp1_s1 + sp2_ha * sp1_s2 +
            sp2_intercept + sp2_stoch_ddvr
        )
      ) %>%
      select(
        !ends_with(c("_cj", "_ca", "_hj", "_ha", "_intercept"))
      )
  } else {
    complete_tseries <- complete_tseries %>%
      left_join(
        kmat_df,
        by = join_by(dd_vr, sp1_name, sp2_name, st_dep_level)
      ) %>%
      mutate(
        sp1_dd_vr = mod_logistic(
          sp1_cj * sp1_s1 + sp1_ca * sp1_s2 + sp1_hj * sp2_s1 + sp1_ha * sp2_s2 +
            sp1_intercept + sp1_stoch_ddvr
        ),
        sp2_dd_vr = mod_logistic(
          sp2_cj * sp2_s1 + sp2_ca * sp2_s2 + sp2_hj * sp1_s1 + sp2_ha * sp1_s2 +
            sp2_intercept + sp2_stoch_ddvr
        )
      ) %>%
      select(
        !ends_with(c("_cj", "_ca", "_hj", "_ha", "_intercept"))
      )
  }
  
  # Merge vital rates
  complete_tseries <- complete_tseries %>%
    left_join(
      spp_vrs_df,
      by = join_by(dd_vr, sp1_name == sp_name)
    ) %>%
    rename_with(
      function(vn) { paste0("sp1_vr_", vn) },
      c(sj, sa, g, r, f)
    ) %>%
    left_join(
      spp_vrs_df,
      by = join_by(dd_vr, sp2_name == sp_name)
    ) %>%
    rename_with(
      function(vn) { paste0("sp2_vr_", vn) },
      c(sj, sa, g, r, f)
    )
  complete_tseries[[paste0("sp1_vr_", ddvr)]] <- complete_tseries$sp1_dd_vr
  complete_tseries[[paste0("sp2_vr_", ddvr)]] <- complete_tseries$sp2_dd_vr
  
  cat("done!\n\tCalculating variables... ")
  
  tseries_part <- complete_tseries %>%
    mutate( # Calculate additional variables
      sp1_n = sp1_s1 + sp1_s2,
      sp2_n = sp2_s1 + sp2_s2
    )
  
  tseries_part$sp1_dd_vr <- tseries_part[[
    paste0("sp1_vr_", ddvr)
  ]]
  tseries_part$sp2_dd_vr <- tseries_part[[
    paste0("sp2_vr_", ddvr)
  ]]
  
  # Calculate realised lambda
  sp1_n_future <- c(tseries_part$sp1_n, NA)
  sp1_n_future <- sp1_n_future[2:length(sp1_n_future)]
  sp1_n_future[seq(t_max+1, length(sp1_n_future), t_max+1)] <- NA
  tseries_part$sp1_real_lambda <- sp1_n_future / tseries_part$sp1_n
  
  sp2_n_future <- c(tseries_part$sp2_n, NA)
  sp2_n_future <- sp2_n_future[2:length(sp2_n_future)]
  sp2_n_future[seq(t_max+1, length(sp2_n_future), t_max+1)] <- NA
  tseries_part$sp2_real_lambda <- sp2_n_future / tseries_part$sp2_n
  
  cat("done!\n")
  
  for(sp1_nam_i in 1:length(sp_names)) {
    sp1_nam <- sp_names[sp1_nam_i]
    for(sp2_nam_i in sp1_nam_i:length(sp_names)) {
      sp2_nam <- sp_names[sp2_nam_i]
      cat(paste0("\t", sp1_nam, " ", sp2_nam, "\n"))
      cat(paste0("\t", format(Sys.time(), "%c"), "\n"))
      
      cat("\t\tCalculating statistical summaries... ")
      
      fit_mod_df_part <- tseries_part %>%
        filter(sp1_name == sp1_nam & sp2_name == sp2_nam) %>%
        group_by(dd_vr, sp1_name, sp2_name, st_dep_level) %>%
        group_modify(calc_ddvr_lm)
      
      if(is.null(fit_mod_df)) {
        fit_mod_df <- fit_mod_df_part
      } else {
        fit_mod_df <- bind_rows(fit_mod_df, fit_mod_df_part)
      }
      
      write_csv(
        fit_mod_df,
        file.path("result_data","ddvr_fit.csv")
      )
      
      cat(" done!\n")
    }
  }
}

complete_tseries <- NULL