# Script for projecting populations based on parameterised models

library("here")
source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_def_models.R"))


# Timestep to start sampling
t_start <- 600
# Timestep to simulate
t_sim <- 100

# Iterate through density-dependent vital rates, species combinations and stage-dependence
# vr_names <- c("sj", "sa", "g", "r", "f")
# sp_names <- c("F2", "F1", "M", "S1", "S2")
# sd_levels <- seq(-1, 1, 0.2)
# sd_levels_str <- as.character(sd_levels)
rand_seeds <- c(1:10, -1 * (1:10))
exp_names <- c("ab", "st_ab", "env_only", "ab_only", "st_ab_only")
n_exps <- length(exp_names)

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
  )

# DF for storing mean squared errors
mse_df <- tibble()

for(ddvr in vr_names) {
# for(ddvr in c("g", "r", "f")) {
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
    filter(time >= t_start & time < t_start + t_sim)
  
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
  
  cat("done!\n\tSimulating estimates & MSEs...\n")
  
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
  
  est_df_org <- tseries_part %>%
    select(!ends_with("_dd_vr")) %>%
    pivot_longer(
      !c(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time),
      names_to = c("foc_sp", ".value"),
      names_pattern = "(sp[0-9]+)_(.+)"
    ) %>%
    left_join(
      fit_mod_df_wider,
      by = join_by(dd_vr, sp1_name, sp2_name, st_dep_level, foc_sp)
    )
  
  est_df_org <- est_df_org %>%
    left_join(
      tseries_part %>% select(
        c(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time) |
          ends_with(c("_s1", "_s2", "_n", "_stoch_ddvr"))
      ),
      by = join_by(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, time),
      suffix = c("_slope", "")
    ) %>%
    mutate(
      stoch_ddvr = ifelse(foc_sp == "sp1", sp1_stoch_ddvr, sp2_stoch_ddvr),
      stoch_ddvr_slope = ifelse(foc_sp == "sp1",
                                sp1_stoch_ddvr_slope,
                                sp2_stoch_ddvr_slope)
    ) %>%
    select(!c("sp1_stoch_ddvr", "sp2_stoch_ddvr",
              "sp1_stoch_ddvr_slope", "sp2_stoch_ddvr_slope"))
  
  # Estimate ddvr & calculate projections & mean squared error
  for(sp1_nam_i in 1:length(sp_names)) {
    sp1_nam <- sp_names[sp1_nam_i]
    for(sp2_nam_i in sp1_nam_i:length(sp_names)) {
      sp2_nam <- sp_names[sp2_nam_i]
      
      cat(paste0("\t\t", sp1_nam, " ", sp2_nam, "\n"))
      cat(paste0("\t\t", format(Sys.time(), "%c"), "\n"))
      
      for(sd_level in sd_levels_str) {
        for(rand_sd in rand_seeds) {
          tseries_small <- est_df_org %>%
            filter(
              dd_vr == ddvr &
                sp1_name == sp1_nam &
                sp2_name == sp2_nam &
                st_dep_level == sd_level &
                rand_seed == rand_sd &
                time >= t_start &
                time < t_start + t_sim
            )
          
          tseries_real <- tseries_small %>%
            select(
              c(n, s1, s2)
            ) %>%
            rename(
              n_real = n,
              s1_real = s1,
              s2_real = s2
            )
          
          tseries_start <- tseries_small %>%
            filter(time == t_start)
          
          spp_n <- tseries_start$n
          spp_s1 <- tseries_start$s1
          spp_s2 <- tseries_start$s2
          
          sim_res_df <- tibble(
            time = t_start,
            foc_sp = c(rep("sp1", n_exps), rep("sp2", n_exps)),
            exp_name = rep(exp_names, 2),
            n = spp_n,
            s1 = spp_s1,
            s2 = spp_s2
          )
          
          exp_vars_mat <- matrix(c(
            tseries_start$sp1_n,
            tseries_start$sp1_s1,
            tseries_start$sp1_s2,
            tseries_start$sp2_n,
            tseries_start$sp2_s1,
            tseries_start$sp2_s2,
            tseries_start$stoch_ddvr,
            rep(1, 2 * n_exps)
          ), ncol = 8)
          
          slope_mat <- matrix(filter_na(c(
            tseries_start$sp1_n_slope,
            tseries_start$sp1_s1_slope,
            tseries_start$sp1_s2_slope,
            tseries_start$sp2_n_slope,
            tseries_start$sp2_s1_slope,
            tseries_start$sp2_s2_slope,
            tseries_start$stoch_ddvr_slope,
            tseries_start[["(Intercept)"]]
          ), 0), ncol = 8)
          
          spp_vrs_mat <- matrix(
            c(
              tseries_start$vr_sj,
              tseries_start$vr_sa,
              tseries_start$vr_g,
              tseries_start$vr_r,
              tseries_start$vr_f
            ),
            ncol = 5
          )
          colnames(spp_vrs_mat) <- vr_names
          
          for(t_i in 1:(t_sim-1)) {
            dd_vr_est <- rowSums(exp_vars_mat * slope_mat)
            
            spp_vrs_mat[,ddvr] <- dd_vr_est
            
            spp_s1_est <- project_juv_vr(
              spp_s1, spp_s2,
              spp_vrs_mat[,"sj"],
              spp_vrs_mat[,"sa"],
              spp_vrs_mat[,"g"],
              spp_vrs_mat[,"r"],
              spp_vrs_mat[,"f"]
            )
            spp_s2_est <- project_ad_vr(
              spp_s1, spp_s2,
              spp_vrs_mat[,"sj"],
              spp_vrs_mat[,"sa"],
              spp_vrs_mat[,"g"],
              spp_vrs_mat[,"r"],
              spp_vrs_mat[,"f"]
            )
            spp_n_est <- spp_s1_est + spp_s2_est
            spp_s1 <- spp_s1_est
            spp_s2 <- spp_s2_est
            spp_n <- spp_n_est
            
            sim_res_df <- sim_res_df %>%
              bind_rows(
                tibble(
                  time = t_start + t_i,
                  foc_sp = c(rep("sp1", n_exps), rep("sp2", n_exps)),
                  exp_name = rep(exp_names, 2),
                  n = spp_n,
                  s1 = spp_s1,
                  s2 = spp_s2
                )
              )
            
            exp_vars_mat <- matrix(c(
              spp_n[rep(1:n_exps, 2)],
              spp_s1[rep(1:n_exps, 2)],
              spp_s2[rep(1:n_exps, 2)],
              spp_n[rep(1:n_exps, 2) + n_exps],
              spp_s1[rep(1:n_exps, 2) + n_exps],
              spp_s2[rep(1:n_exps, 2) + n_exps],
              tseries_small[(1:(n_exps*2)) + t_i * n_exps * 2,]$stoch_ddvr,
              rep(1, 2 * n_exps)
            ), ncol = 8)
          }
          
          sim_res_df <- bind_cols(sim_res_df, tseries_real)
          
          sim_summ_df <- sim_res_df %>%
            group_by(foc_sp, exp_name) %>%
            summarise(
              n_mse = mean((n_real - n)^2),
              s1_mse = mean((s1_real - s1)^2),
              s2_mse = mean((s2_real - s2)^2),
              n_mape = mean(abs((n_real - n) / n_real)) * 100,
              s1_mape = mean(abs((s1_real - s1) / s1_real)) * 100,
              s2_mape = mean(abs((s2_real - s2) / s2_real)) * 100
            ) %>%
            mutate(
              dd_vr = ddvr,
              sp1_name = sp1_nam,
              sp2_name = sp2_nam,
              st_dep_level = sd_level,
              rand_seed = rand_sd
            ) %>%
            relocate(
              dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed,
              .before = everything()
            )
          
          mse_df <- bind_rows(
            mse_df,
            sim_summ_df
          )
        }
      }
      
      write_csv(
        mse_df,
        here("StageDependentComp","result_data","wp3","mse_long_summ.csv")
      )
    }
  }
}

complete_tseries <- NULL
tseries_part <- NULL
est_df_org <- NULL


