################################################################################
#' Run script for calculating mean & variances and statistical summaries
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

library("here")
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_init_setup.R"))

dd_vrs <- c("sj", "sa", "g", "r", "f")

# Timestep to start sampling
t_start <- 600
# Timestep to end sampling
t_end <- 2000

# Iterate through density-dependent vital rates, species combinations and stage-dependence
dd_vrs <- c("sj", "sa", "g", "r", "f")
sp_names <- c("F2", "F1", "M", "S1", "S2")
sd_levels <- seq(-1, 1, 0.2)
rand_seeds <- c(1:10, -1 * (1:10))
vr_sd <- 0.2

# Load interaction coefficients
kmat_df <- read_csv(
  here("StageDependentComp","result_data","wp3","inter_coeffs.csv")
)

# Variables to calculate means and variances of
meanvar_names <- c(
  "n", "s1", "s2", "struct",
  "dd_vr", "dd_vr_slope", "lambda", "lambda_slope", "real_lambda", "dl"
)
# meanvar_names <- c(
#   "dd_vr_slope", "lambda_slope"
# )
meanvar_df <- NULL

# DF to store statistical summaries
stat_summ_df <- NULL


for(ddvr in dd_vrs) {
  cat(paste("Density-dependent vital rate:", ddvr, "\n"))
  cat(paste0(format(Sys.time(), "%c"), "\n"))
  
  cat("Loading data... ")
  
  complete_tseries <- read_csv(
    here("StageDependentComp","result_data","wp3",
         paste0("sim_tseries_", ddvr, ".csv")),
    col_types = cols(dd_vr = col_character())
  )
  
  cat("done!\n")
  
  cat("\tCalculating variables... ")
  
  tseries_part <- complete_tseries %>%
    group_by(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed) %>%
    mutate( # Calculate additional variables
      sp1_n = sp1_s1 + sp1_s2,
      sp2_n = sp2_s1 + sp2_s2,
      sp1_struct = sp1_s1 / (sp1_s1 + sp1_s2),
      sp2_struct = sp2_s1 / (sp2_s1 + sp2_s2),
      sp1_dd_vr = .data[[
        paste0("sp1_vr_", ddvr)
      ]],
      sp2_dd_vr = .data[[
        paste0("sp2_vr_", ddvr)
      ]],
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
      sp1_damping_r = damping_ratio_list(vrs_to_mpms_vec(
        sp1_vr_sj,
        sp1_vr_sa,
        sp1_vr_g,
        sp1_vr_r,
        sp1_vr_f
      )),
      sp2_damping_r = damping_ratio_list(vrs_to_mpms_vec(
        sp2_vr_sj,
        sp2_vr_sa,
        sp2_vr_g,
        sp2_vr_r,
        sp2_vr_f
      )),
      sp1_dd_vr_intercept = .data[[
        paste0("sp1_stoch_", ddvr)
      ]] * vr_sd + kmat_df[
        kmat_df$dd_vr == ddvr &
          kmat_df$sp1_name == sp1_name &
          kmat_df$sp2_name == sp2_name,
      ]$sp1_intercept,
      sp2_dd_vr_intercept = .data[[
        paste0("sp2_stoch_", ddvr)
      ]] * vr_sd + kmat_df[
        kmat_df$dd_vr == ddvr &
          kmat_df$sp1_name == sp1_name &
          kmat_df$sp2_name == sp2_name,
      ]$sp2_intercept
    ) %>%
    mutate(
      sp1_dd_vr_slope = ifelse(
        dd_vr == "f",
        exp_slope(sp1_dd_vr_intercept),
        logistic_slope(sp1_dd_vr_intercept)
      ),
      sp2_dd_vr_slope = ifelse(
        dd_vr == "f",
        exp_slope(sp2_dd_vr_intercept),
        logistic_slope(sp2_dd_vr_intercept)
      )
    ) %>%
    mutate(
      sp1_lambda_slope = sp1_dl * sp1_dd_vr_slope,
      sp2_lambda_slope = sp2_dl * sp2_dd_vr_slope
    )
  
  tseries_part_future <- tseries_part %>%
    mutate(time = time - 1) %>%
    select(
      c(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time, sp1_n, sp2_n)
    ) %>%
    rename_with(
      function(vn) { paste0(vn, "_future") },
      c(sp1_n, sp2_n)
    )
  
  tseries_part <- tseries_part %>%
    left_join(tseries_part_future, by = join_by(
      dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time
    )) %>%
    mutate(
      sp1_real_lambda = sp1_n_future / sp1_n,
      sp2_real_lambda = sp2_n_future / sp2_n
    )
  
  tseries_part_future <- NULL
  
  cat("done!\n\tCalculating means and variances... ")
  
  meanvar_df_part <- tseries_part %>%
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
    here("StageDependentComp","result_data","wp3","meanvar_temp.csv")
  )
  
  cat("done!\n")
  
  for(sp1_nam_i in 1:length(sp_names)) {
    sp1_nam <- sp_names[sp1_nam_i]
    for(sp2_nam_i in sp1_nam_i:length(sp_names)) {
      sp2_nam <- sp_names[sp2_nam_i]
      cat(paste("\t",sp1_nam,sp2_nam,"\n"))
      cat(paste0("\t", format(Sys.time(), "%c"), "\n"))
      
      cat("\tCalculating statistical summaries...\n")

      stat_summ_df_part <- tseries_part %>%
        filter(sp1_name == sp1_nam & sp2_name == sp2_nam) %>%
        group_modify(calc_stat_summs)

      if(is.null(stat_summ_df)) {
        stat_summ_df <- stat_summ_df_part
      } else {
        stat_summ_df <- bind_rows(stat_summ_df, stat_summ_df_part)
      }

      write_csv(
        stat_summ_df,
        here("StageDependentComp","result_data","wp3","stat_summ_temp.csv")
      )

      cat("\tDone!\n")
    }
  }
}






