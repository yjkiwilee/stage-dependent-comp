# Script for projecting populations based on parameterised models

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
exp_names <- c("ab", "st_ab")

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

# Get ddvr fit parameters
fit_mod_df <- read_csv(
  here("StageDependentComp","result_data","wp3",
       "ddvr_fit.csv"),
  col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
)

fit_mod_df_wider <- fit_mod_df %>%
  mutate(
    foc_sp = str_match(resp_name, "(sp[0-9]+)_dd_vr")[,2]
  ) %>%
  select(
    !c(resp_name, n) & !starts_with("var_") | c(var_name, var_est)
  ) %>%
  pivot_wider(
    names_from = var_name,
    values_from = var_est
  ) %>%
  mutate(
    stoch_ddvr = ifelse(is.na(sp1_stoch_ddvr), sp2_stoch_ddvr, sp1_stoch_ddvr)
  )

# DF for storing mean squared errors
mse_df <- NULL

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
  
  cat("done!\n\tCalculating estimates & MSEs...\n")
  
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
  
  # Get future n & s1 & s2
  tseries_future <- tseries_part %>%
    mutate(
      time = time - 1
    ) %>%
    select(
      ends_with(c("_n", "_s1", "_s2")) | c(dd_vr, sp1_name, sp2_name, st_dep_level,
                                        rand_seed, time)
    )
  tseries_part <- tseries_part %>%
    left_join(
      tseries_future,
      by = join_by(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time),
      suffix = c("", "_future")
    )
  tseries_future <- NULL
  
  # Estimate ddvr & calculate projections & mean squared error
  for(sp1_nam_i in 1:length(sp_names)) {
    sp1_nam <- sp_names[sp1_nam_i]
    for(sp2_nam_i in sp1_nam_i:length(sp_names)) {
      sp2_nam <- sp_names[sp2_nam_i]
      cat(paste0("\t\t", sp1_nam, " ", sp2_nam, "\n"))
      cat(paste0("\t\t", format(Sys.time(), "%c"), "\n"))
      
      est_df_org <- tseries_part %>%
        filter(sp1_name == sp1_nam & sp2_name == sp2_nam) %>%
        select(!ends_with("_dd_vr")) %>%
        pivot_longer(
          !c(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time),
          names_to = c("foc_sp", ".value"),
          names_pattern = "(sp[0-9]+)_(.+)"
        ) %>%
        left_join(
          fit_mod_df_wider,
          by = join_by(dd_vr, sp1_name, sp2_name, st_dep_level, foc_sp)
        ) %>%
        left_join(
          tseries_part %>% select(
            c(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time) |
              ends_with(c("_s1", "_s2", "_n"))
          ),
          by = join_by(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time),
          suffix = c("_slope", "")
        )
      
      for(exp_name in exp_names) {
        est_df <- est_df_org
        
        if(exp_name == "ab") {
          est_df <- est_df %>% mutate(
            dd_vr_est = `(Intercept)` +
              sp1_n_slope * sp1_n +
              sp2_n_slope * sp2_n
          )
        } else if(exp_name == "st_ab") {
          est_df <- est_df %>% mutate(
            dd_vr_est = `(Intercept)` +
              sp1_s1_slope * sp1_s1 +
              sp1_s2_slope * sp1_s2 +
              sp2_s1_slope * sp2_s1 +
              sp2_s2_slope * sp2_s2
          )
        }
        
        est_df[[paste0("vr_", ddvr)]] <-
          est_df$dd_vr_est
        
        est_df <- est_df %>%
          mutate(
            s1_est = project_juv_vr(
              s1, s2, vr_sj, vr_sa, vr_g, vr_r, vr_f
            ),
            s2_est = project_ad_vr(
              s1, s2, vr_sj, vr_sa, vr_g, vr_r, vr_f
            )
          ) %>%
          mutate(
            n_est = s1_est + s2_est
          ) %>%
          mutate(
            s1_err = s1_future - s1_est,
            s2_err = s2_future - s2_est,
            n_err = n_future - n_est
          )
        
        est_summ_df <- est_df %>%
          group_by(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, foc_sp) %>%
          summarise(
            s1_mse = sum(s1_err * s1_err, na.rm = TRUE) / sum(!is.na(s1_err)),
            s2_mse = sum(s2_err * s2_err, na.rm = TRUE) / sum(!is.na(s2_err)),
            n_mse = sum(n_err * n_err, na.rm = TRUE) / sum(!is.na(n_err)),
            n_mape = mean(abs(n_err / n_future), na.rm = TRUE) * 100,
            s1_mape = mean(abs(s1_err / s1_future), na.rm = TRUE) * 100,
            s2_mape = mean(abs(s2_err / s2_future), na.rm = TRUE) * 100
          ) %>%
          mutate(
            exp_name = exp_name
          ) %>%
          relocate(
            exp_name, .after = rand_seed
          )
        
        est_df <- NULL
        
        if(is.null(mse_df)) {
          mse_df <- est_summ_df
        } else {
          mse_df <- bind_rows(mse_df, est_summ_df)
        }
      }
      
      write_csv(
        mse_df,
        here("StageDependentComp","result_data","wp3","mse_summ.csv")
      )
    }
  }
}

complete_tseries <- NULL
tseries_part <- NULL


