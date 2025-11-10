################################################################################
#' Run script for calculating mean & variances and statistical summaries
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

library("here")
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_def_models.R"))

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

# Variables to calculate means and variances of
meanvar_names <- c(
  "n", "s1", "s2", "struct",
  "dd_vr", "lambda", "real_lambda", "dl",
  "damping_r"
)
# meanvar_names <- c(
#   "dd_vr_slope", "lambda_slope"
# )
meanvar_df <- NULL

# DF to store statistical summaries
stat_summ_df <- NULL

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
    here("StageDependentComp","result_data","wp3",
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
      sp2_n = sp2_s1 + sp2_s2,
      sp1_struct = sp1_s1 / (sp1_s1 + sp1_s2),
      sp2_struct = sp2_s1 / (sp2_s1 + sp2_s2),
      sp1_lambda = lambda_from_vr(
        sp1_vr_sj,
        sp1_vr_sa,
        sp1_vr_g,
        sp1_vr_r,
        sp1_vr_f
      ),
      sp2_lambda = lambda_from_vr(
        sp2_vr_sj,
        sp2_vr_sa,
        sp2_vr_g,
        sp2_vr_r,
        sp2_vr_f
      ),
      sp1_dl = d_lambda_vec(
        ddvr,
        sp1_vr_sj,
        sp1_vr_sa,
        sp1_vr_g,
        sp1_vr_r,
        sp1_vr_f
      ),
      sp2_dl = d_lambda_vec(
        ddvr,
        sp2_vr_sj,
        sp2_vr_sa,
        sp2_vr_g,
        sp2_vr_r,
        sp2_vr_f
      ),
      sp1_damping_r = damping_ratio_vr(
        sp1_vr_sj,
        sp1_vr_sa,
        sp1_vr_g,
        sp1_vr_r,
        sp1_vr_f
      ),
      sp2_damping_r = damping_ratio_vr(
        sp2_vr_sj,
        sp2_vr_sa,
        sp2_vr_g,
        sp2_vr_r,
        sp2_vr_f
      )
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
  
  cat("done!\n\tCalculating means and variances... ")
  
  meanvar_df_part <- tseries_part %>%
    group_by(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed) %>%
    summarise(
      across(
        all_of(c(paste0("sp1_",meanvar_names), paste0("sp2_",meanvar_names))),
        list(
          mean = function(x) { mean(x, na.rm = TRUE) },
          var = function(x) { var(x, na.rm = TRUE) }
        )
      ),
      .groups = "keep"
    )
  
  cat("done!\n\tSaving results... ")
  
  if(is.null(meanvar_df)) {
    meanvar_df <- meanvar_df_part
  } else {
    meanvar_df <- bind_rows(meanvar_df, meanvar_df_part)
  }
  
  # Store results
  write_csv(
    meanvar_df,
    here("StageDependentComp","result_data","wp3","meanvar.csv")
  )
  
  cat("done!\n")
  
  for(sp1_nam_i in 1:length(sp_names)) {
    sp1_nam <- sp_names[sp1_nam_i]
    for(sp2_nam_i in sp1_nam_i:length(sp_names)) {
      sp2_nam <- sp_names[sp2_nam_i]
      cat(paste0("\t", sp1_nam, " ", sp2_nam, "\n"))
      cat(paste0("\t", format(Sys.time(), "%c"), "\n"))
      
      cat("\t\tCalculating statistical summaries... ")

      stat_summ_df_part <- tseries_part %>%
        filter(sp1_name == sp1_nam & sp2_name == sp2_nam) %>%
        group_by(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed) %>%
        group_modify(calc_stat_summs)

      if(is.null(stat_summ_df)) {
        stat_summ_df <- stat_summ_df_part
      } else {
        stat_summ_df <- bind_rows(stat_summ_df, stat_summ_df_part)
      }

      write_csv(
        stat_summ_df,
        here("StageDependentComp","result_data","wp3","stat_summ.csv")
      )

      cat(" done!\n")
    }
  }
}

# 
# # Process meanvar to create a column for the difference in CV of each stage
# meanvar_df <- read_csv(
#   here("StageDependentComp","result_data","wp3","meanvar.csv")
# )
# 
# meanvar_df %>%
#   mutate(
#     sp1_s_cv_diff = sp1_s1_mean / sqrt(sp1_s1_var) - 
#       sp1_s2_mean / sqrt(sp1_s2_var),
#     sp2_s_cv_diff = sp2_s1_mean / sqrt(sp2_s1_var) -
#       sp2_s2_mean / sqrt(sp2_s2_var)
#   ) %>%
#   select(ends_with("s_cv_diff"))


