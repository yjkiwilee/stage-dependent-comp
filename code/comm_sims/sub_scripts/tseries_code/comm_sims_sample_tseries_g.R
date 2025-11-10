################################################################################
#' Run script for collecting time series data
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

require("here")

source(here("StageDependentComp","code","comm_sims","sub_scripts",
            "comm_sims_def_models.R"))

# === Calculate stochastic terms ===

# Timestep to end sampling
t_end <- 2000

# Number of replicates
n_reps <- 20
rand_seeds <- c(1:(n_reps / 2), -1 * (1:(n_reps / 2)))

stoch_rnorms_seeds <- list()

for(rand_seed in rand_seeds) {
  set.seed(abs(rand_seed))
  stoch_rnorms <- matrix(
    rnorm(t_end * 2 * (2+1+1+1)), ncol = 2 * (2+1+1+1)
  )
  stoch_rnorms <- rbind(stoch_rnorms, rep(NA, 2 * (2+1+1+1)))
  
  # Negative random seed indicates swapping stochastic terms between spp
  if(rand_seed < 0) {
    stoch_temp <- stoch_rnorms[,1:(2+1+1+1)]
    stoch_rnorms[,1:(2+1+1+1)] <- stoch_rnorms[,(2+1+1+1+1):(2*(2+1+1+1))]
    stoch_rnorms[,(2+1+1+1+1):(2*(2+1+1+1))] <- stoch_temp
  }
  
  colnames(stoch_rnorms) <- c(
    "stoch_sj_1",
    "stoch_sa_1",
    "stoch_g_1",
    "stoch_r_1",
    "stoch_f_1",
    "stoch_sj_2",
    "stoch_sa_2",
    "stoch_g_2",
    "stoch_r_2",
    "stoch_f_2"
  )
  stoch_rnorms_seeds[[length(stoch_rnorms_seeds) + 1]] <- stoch_rnorms
}

# === Calculate time series data ===

# DF to store the time series data
tseries_df <- NULL

# Test each density-dependent vital rate
# for(dd_vr_i in 1:length(vr_names)) {
for(dd_vr_i in c(3)) {
  dd_vr <- vr_names[dd_vr_i]
  
  cat(paste0("Testing ", dd_vr, "\n"))
  
  # Calculate & extract stochasticity term of ddvr
  ddvr_stoch_rnorms <-
    stoch_rnorms[,c(dd_vr_i, length(vr_names) + dd_vr_i)] *
    vr_sd[dd_vr]
  
  colnames(ddvr_stoch_rnorms) <- c(
    "sp1_stoch_ddvr",
    "sp2_stoch_ddvr"
  )
  
  # Test each species pair
  for(sp_comb_i in 1:length(sp_combs)) {
    sp_comb <- sp_combs[sp_comb_i]
    
    cat(paste0("\tTesting ", dd_vr, " ", sp_comb, "\n"))
    
    # Get the names of two species
    spp_names <- str_match(sp_comb, "^([A-Z][0-9]*)([A-Z][0-9]*)$")
    sp1_name <- spp_names[1,2]
    sp2_name <- spp_names[1,3]
    
    # Get equilibrium size
    sp1_eq_pop <- as.numeric(pop_eq_avg[
      pop_eq_avg$dd_vr == ddvr &
        pop_eq_avg$sp_name == sp1_name,
    ][1, c("mean_j", "mean_a")])
    sp2_eq_pop <- as.numeric(pop_eq_avg[
      pop_eq_avg$dd_vr == ddvr &
        pop_eq_avg$sp_name == sp2_name,
    ][1, c("mean_j", "mean_a")])
    
    # Replicates
    for(rep_i in 1:n_reps) {
      cat(paste0("\t\tReplicate no. ", rep_i, "\n"))
      cat(paste0("\t\t", format(Sys.time(), "%c"), "\n"))
        
      # Generate stochasticity terms
      rand_seed <- rand_seeds[rep_i]
      stoch_rnorms <- stoch_rnorms_seeds[[rep_i]]
      cat(paste0("\t\tRandom seed: ", rand_seed, "\n"))
      
      # Organise stochasticity terms for each species
      stoch_rnorms_sp1 <- lapply(1:t_end, function(t) {
        list(
          s = c(stoch_rnorms[t, c(1,2)]), # Small & large survival
          g = stoch_rnorms[t, 3], # Maturation rate
          r = stoch_rnorms[t, 4], # Retrogression rate
          f = stoch_rnorms[t, 5] # Fecundity
        )
      })
      stoch_rnorms_sp2 <- lapply(1:t_end, function(t) {
        list(
          s = c(stoch_rnorms[t, c(6,7)]), # Small & large survival
          g = stoch_rnorms[t, 8], # Maturation rate
          r = stoch_rnorms[t, 9], # Retrogression rate
          f = stoch_rnorms[t, 10] # Fecundity
        )
      })
      
      # Test all levels of stage-dependence
      for(sd_level_i in 1:length(sd_levels)) {
        sd_level_str <- sd_levels_str[sd_level_i]
        
        cat(paste0("\t\t\tTesting stage-dependence ", sd_level_str, "\n"))
        
        comm_model <- sd_comm_models[[dd_vr]][[sp_comb]][[sd_level_str]]
        
        # Project community
        comm_project <- sdcomp_project(
          comm_model,
          sp1_eq_pop,
          sp2_eq_pop,
          timestep = t_end,
          stoch_rnorms_sp1, stoch_rnorms_sp2
        )
        
        # Merge DF with stochasticity terms
        project_res_df <- bind_cols(comm_project, ddvr_stoch_rnorms) %>%
          mutate( # Insert dd vital rate, species identity & stage dependence
            dd_vr = dd_vr,
            sp1_name = sp1_name,
            sp2_name = sp2_name,
            st_dep_level = sd_level_str,
            rand_seed = rand_seed
          ) %>%
          select(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed,
                 time, everything())
        
        # Merge
        if(is.null(tseries_df)) {
          tseries_df <- project_res_df
        } else {
          tseries_df <- rbind(tseries_df, project_res_df)
        }
        
        # Purge DF
        project_res_df <- NULL
      }
      
      # Save DF so far into a file
      if(file.exists(here("StageDependentComp","result_data","wp3",
                          paste0("sim_tseries_", dd_vr, ".csv")))) {
        write_csv(
          tseries_df,
          here("StageDependentComp","result_data","wp3",
               paste0("sim_tseries_", dd_vr, ".csv")),
          append = TRUE
        )
      } else {
        write_csv(
          tseries_df,
          here("StageDependentComp","result_data","wp3",
               paste0("sim_tseries_", dd_vr, ".csv"))
        )
      }
      
      # Purge DF
      tseries_df <- NULL
    }
  }
}
