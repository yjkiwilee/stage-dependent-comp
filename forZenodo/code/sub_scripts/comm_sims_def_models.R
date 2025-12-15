################################################################################
#' Script for defining the virtual species
#'
#' Written by Young Jun Lee
#' Mar 2025

# ======= Define virtual species =======

source(file.path("code","sub_scripts","comm_sims_init_setup.R"))

# Constant parameters
## Number of life cycle stages
n_stages <- 2
## Names of vital rates
vr_names <- c("sj", "sa", "g", "r", "f")
## Names of species
sp_names <- c("F2", "F1", "M", "S1", "S2")
## Baseline survival, growth and reproduction of each species
spp_vrs_base <- lapply(seq(0, 1, length.out = 5), function(x) {
  vrs <- list(
    s = c( # Survival
      lerp(x, 0.4, 0.9),
      lerp(x, 0.6, 0.95)
    ),
    g = c( # Progression
      lerp(x, 0.9, 0.3)
    ),
    r = c( # Retrogression
      lerp(x, 0.05, 0.1)
    )
  )
  # Return vital rates
  return(vrs)
})
# spp_vrs_base <- lapply(seq(0, 1, length.out = 5), function(x) {
#   vrs <- list(
#     s = c( # Survival
#       lerp(x, 0.3, 0.9),
#       lerp(x, 0.6, 0.95)
#     ),
#     g = c( # Progression
#       lerp(x, 0.12, 0.02)
#     ),
#     r = c( # Retrogression
#       lerp(x, 0.006, 0.06)
#     )
#   )
#   # Return vital rates
#   return(vrs)
# })
names(spp_vrs_base) <- sp_names
## Minimum values of baseline vital rates
min_vrs <- c(
  sj = 0.4,
  sa = 0.6,
  g = 0.3,
  r = 0.05,
  f = 0.1
)
## Standard deviation of stochasticity in vital rate
vr_sd <- min_vrs * 0.1
## Target equilibrium population size
n_eq <- 500
## Intra/interspecific interaction weights
w_intra <- min_vrs * 0.3 / (n_eq * (1 + 1/4.5))
w_inter <- w_intra / 4.5
w_intra
## Density-dependence sign coefficients
dd_sign <- c(
  sj = -1,
  sa = -1,
  g = -1,
  r = 1,
  f = -1
)
## Combinations of species
sp_combs <- c()
for(sp1_i in 1:length(sp_names)) {
  for(sp2_i in sp1_i:length(sp_names)) {
    sp_combs <- c(sp_combs, paste0(sp_names[sp1_i], sp_names[sp2_i]))
  }
}
## Stage dependence levels
sd_levels <- seq(-1, 1, 0.2)
sd_levels_str <- as.character(sd_levels)

# Update fecundity so that lambda = 1 at the desired equilibrium population sizes
spp_vrs <- lapply(vr_names, function(ddvr_name) {
  lapply(spp_vrs_base, function(base_vrs) {
    update_f_for_eq(
      base_vrs, n_eq, ddvr_name, w_intra[ddvr_name], w_inter[ddvr_name]
    )
  })
})
names(spp_vrs) <- vr_names

spp_vrs

spp_vrs_list <- lapply(names(spp_vrs), function(ddvr) {
  lapply(names(spp_vrs[[ddvr]]), function(sp_name) {
    vrs <- spp_vrs[[ddvr]][[sp_name]]
    
    sp_mpm <- vrs_to_mpm(n_stages, vrs)
    
    spp_vrs_row <- tibble(
      dd_vr = ddvr,
      sp_name = sp_name,
      sj = vrs$s[1],
      sa = vrs$s[2],
      g = vrs$g[1],
      r = vrs$r[1],
      f = vrs$f[1],
      m_1_1 = sp_mpm$matA[1,1],
      m_1_2 = sp_mpm$matA[1,2],
      m_2_1 = sp_mpm$matA[2,1],
      m_2_2 = sp_mpm$matA[2,2],
      gen_time = generation.time(sp_mpm$matA),
      net_rep_rate = net.reproductive.rate(sp_mpm$matA),
      life_expect = life_expect_mean(sp_mpm$matU),
      rep_val_j = reproductive.value(sp_mpm$matA)[1],
      rep_val_a = reproductive.value(sp_mpm$matA)[2],
      sens_sj = d_lambda("sj", vrs),
      sens_sa = d_lambda("sa", vrs),
      sens_g = d_lambda("g", vrs),
      sens_r = d_lambda("r", vrs),
      sens_f = d_lambda("f", vrs)
    )
  })
})

write_csv(
  bind_rows(spp_vrs_list),
  file.path("result_data","spp_vrs.csv")
)

# ======= Define stage-independent community models =======

# List of lists to store base community models
comm_models <- list(
  sj = NULL,
  sa = NULL,
  g = NULL,
  r = NULL,
  f = NULL
)

# For each density-dependent vital rate:
for(ddvr in vr_names) {
  comm_models[[ddvr]] <- list()
  
  # Define stochasticity magnitudes
  vr_sds <- list(
    s = c(0, 0), # Small & large survival
    g = c(0), # Maturation rate
    r = c(0), # Retrogression rate
    f = c(0) # Fertility
  )
  
  if(ddvr == "sj") {
    vr_sds$s[1] <- vr_sd[ddvr]
  } else if(ddvr == "sa") {
    vr_sds$s[2] <- vr_sd[ddvr]
  } else {
    vr_sds[[ddvr]][1] <- vr_sd[ddvr]
  }
  
  # Iterate through pairs of virtual species
  for(sp1_i in 1:length(sp_names)) {
    for(sp2_i in sp1_i:length(sp_names)) {
      sp1_name <- sp_names[sp1_i]
      sp2_name <- sp_names[sp2_i]
      
      # Initialise interaction matrix
      k_mat <- list(
        s = matrix(c(
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0
        ), ncol = 4, byrow = TRUE),
        g = matrix(c( 
          0, 0, 0, 0,
          0, 0, 0, 0
        ), ncol = 4, byrow = TRUE),
        r = matrix(c(
          0, 0, 0, 0,
          0, 0, 0, 0
        ), ncol = 4, byrow = TRUE),
        f = matrix(c(
          0, 0, 0, 0,
          0, 0, 0, 0
        ), ncol = 4, byrow = TRUE)
      )
      
      w_c_base <- w_intra[ddvr] * dd_sign[[ddvr]]
      w_h_base <- w_inter[ddvr] * dd_sign[[ddvr]]
      
      if(ddvr == "sj") {
        k_mat$s[1,c(1,2)] <- w_c_base
        k_mat$s[1,c(3,4)] <- w_h_base
        k_mat$s[3,c(1,2)] <- w_h_base
        k_mat$s[3,c(3,4)] <- w_c_base
      } else if(ddvr == "sa") {
        k_mat$s[2,c(1,2)] <- w_c_base
        k_mat$s[2,c(3,4)] <- w_h_base
        k_mat$s[4,c(1,2)] <- w_h_base
        k_mat$s[4,c(3,4)] <- w_c_base
      } else {
        k_mat[[ddvr]][1,c(1,2)] <- w_c_base
        k_mat[[ddvr]][1,c(3,4)] <- w_h_base
        k_mat[[ddvr]][2,c(1,2)] <- w_h_base
        k_mat[[ddvr]][2,c(3,4)] <- w_c_base
      }
      
      # Define community model
      comm_models[[ddvr]][[paste0(sp1_name, sp2_name)]] <-
        def_sdcomp_model(
          n_stages,
          spp_vrs[[ddvr]][[sp1_name]],
          spp_vrs[[ddvr]][[sp2_name]],
          k_mat,
          vr_sds, vr_sds
        )
    }
  }
}

# ====== Calculate equilibrium population structures ======

eq_struct_tmin <- 1000
eq_struct_tmax <- 2000
eq_struct_nrep <- 10

set.seed(100)

pop_eq_structs_df <- tibble(
  dd_vr = character(),
  sp1_name = character(),
  sp2_name = character(),
  sp1_juvenile = numeric(),
  sp1_adult = numeric(),
  sp1_n = numeric(),
  sp1_j_norm = numeric(),
  sp1_a_norm = numeric(),
  sp2_juvenile = numeric(),
  sp2_adult = numeric(),
  sp2_n = numeric(),
  sp2_j_norm = numeric(),
  sp2_a_norm = numeric()
)

if(file.exists(file.path("result_data","pop_eq_structs.csv"))) {
  pop_eq_structs_df <- read_csv(
    file.path("result_data","pop_eq_structs.csv")
  )
} else {
  for(dd_vr_name in vr_names) {
    for(sp1_i in 1:length(sp_names)) {
      sp1_name <- sp_names[sp1_i]
      for(sp2_i in sp1_i:length(sp_names)) {
        sp2_name <- sp_names[sp2_i]
        cat(paste("Processing",dd_vr_name,sp1_name,sp2_name,"\n"))
        cat(paste0(format(Sys.time(), "%c"), "\n"))
        
        eq_sp1_j <- 0
        eq_sp1_a <- 0
        eq_sp2_j <- 0
        eq_sp2_a <- 0
        
        # Repeat simulations nrep times
        for(i in 1:eq_struct_nrep) {
          temp_comm_proj <- sdcomp_project(
            comm_models[[dd_vr_name]][[paste0(sp1_name, sp2_name)]],
            c(n_eq / 2, n_eq / 2), c(n_eq / 2, n_eq / 2),
            timestep = eq_struct_tmax
          )
          
          # Add average stage-specific abundances
          eq_sp1_j <- eq_sp1_j + mean(temp_comm_proj$sp1_s1[eq_struct_tmin:eq_struct_tmax])
          eq_sp1_a <- eq_sp1_a + mean(temp_comm_proj$sp1_s2[eq_struct_tmin:eq_struct_tmax])
          eq_sp2_j <- eq_sp2_j + mean(temp_comm_proj$sp2_s1[eq_struct_tmin:eq_struct_tmax])
          eq_sp2_a <- eq_sp2_a + mean(temp_comm_proj$sp2_s2[eq_struct_tmin:eq_struct_tmax])
        }
        
        eq_sp1_j <- eq_sp1_j / eq_struct_nrep
        eq_sp1_a <- eq_sp1_a / eq_struct_nrep
        eq_sp2_j <- eq_sp2_j / eq_struct_nrep
        eq_sp2_a <- eq_sp2_a / eq_struct_nrep
        
        
        pop_eq_structs_row <- tibble(
          dd_vr = dd_vr_name,
          sp1_name = sp1_name,
          sp2_name = sp2_name,
          sp1_juvenile = eq_sp1_j,
          sp1_adult = eq_sp1_a,
          sp1_n = eq_sp1_j + eq_sp1_a,
          sp1_j_norm = eq_sp1_j / (eq_sp1_j + eq_sp1_a),
          sp1_a_norm = eq_sp1_a / (eq_sp1_j + eq_sp1_a),
          sp2_juvenile = eq_sp2_j,
          sp2_adult = eq_sp2_a,
          sp2_n = eq_sp1_j + eq_sp2_a,
          sp2_j_norm = eq_sp2_j / (eq_sp2_j + eq_sp2_a),
          sp2_a_norm = eq_sp2_a / (eq_sp2_j + eq_sp2_a)
        )
        
        pop_eq_structs_df <- bind_rows(pop_eq_structs_df, pop_eq_structs_row)
        
        write_csv(
          pop_eq_structs_df,
          file.path("result_data","pop_eq_structs.csv")
        )
      }
    }
  }
}

# Calculate average equilibrium structures from runs

pop_eq_avg <- tibble(
  dd_vr = rep(pop_eq_structs_df$dd_vr, times = 2),
  sp_name = c(pop_eq_structs_df$sp1_name, pop_eq_structs_df$sp2_name),
  j = c(pop_eq_structs_df$sp1_juvenile, pop_eq_structs_df$sp2_juvenile),
  j_norm = c(pop_eq_structs_df$sp1_j_norm, pop_eq_structs_df$sp2_j_norm),
  a = c(pop_eq_structs_df$sp1_adult, pop_eq_structs_df$sp2_adult),
  a_norm = c(pop_eq_structs_df$sp1_a_norm, pop_eq_structs_df$sp2_a_norm),
  n = c(pop_eq_structs_df$sp1_n, pop_eq_structs_df$sp2_n)
) %>%
  mutate(
    sp_name = factor(sp_name, levels = sp_names)
  ) %>%
  group_by(dd_vr, sp_name) %>%
  summarise(
    mean_j = mean(j),
    mean_a = mean(a),
    mean_j_norm = mean(j_norm),
    mean_a_norm = mean(a_norm),
    mean_n = mean(n),
    .groups='keep'
  )

write_csv(
  pop_eq_avg,
  file.path("result_data","pop_eq_structs_avg.csv")
)

# pop_eq_avg_an <- NULL
# 
# for(ddvr in vr_names) {
#   for(sp_name in sp_names) {
#     # Calculate density dependent vital rate for stationary
#     stat_ddvr <- dd_vr_for_stationary(spp_vrs[[ddvr]][[sp_name]], ddvr)
#     
#     stat_vrs <- spp_vrs[[ddvr]][[sp_name]]
#     if(ddvr == "sj") {
#       stat_vrs$s[1] <- stat_ddvr
#     } else if(ddvr == "sa") {
#       stat_vrs$s[2] <- stat_ddvr
#     } else {
#       stat_vrs[[ddvr]][1] <- stat_ddvr
#     }
#     
#     stat_mpm <- vrs_to_mpm(n_stages, stat_vrs)
#     
#     eq_struct <- stable.stage(stat_mpm$matA)
#     
#     pop_eq_avg_an_row <- tibble(
#       dd_vr = ddvr,
#       sp_name = sp_name,
#       mean_j = eq_struct[1] * n_eq,
#       mean_a = eq_struct[2] * n_eq,
#       mean_j_norm = eq_struct[1],
#       mean_a_norm = eq_struct[2],
#       mean_n = n_eq
#     )
#     
#     if(is.null(pop_eq_avg_an)) {
#       pop_eq_avg_an <- pop_eq_avg_an_row
#     } else {
#       pop_eq_avg_an <- bind_rows(pop_eq_avg_an, pop_eq_avg_an_row)
#     }
#   }
# }
# 
# pop_eq_avg_an

# ====== Calculate stage-dependent models ======

# List of lists to store stage-dependent community models
sd_comm_models <- list(
  sj = NULL,
  sa = NULL,
  g = NULL,
  r = NULL,
  f = NULL
)

# DF for storing interaction coefficients
kmat_df <- NULL

# For each density-dependent vital rate:
for(ddvr in vr_names) {
  sd_comm_models[[ddvr]] <- list()
  # Iterate through pairs of virtual species
  for(sp1_i in 1:length(sp_names)) {
    for(sp2_i in sp1_i:length(sp_names)) {
      sp1_name <- sp_names[sp1_i]
      sp2_name <- sp_names[sp2_i]
      sd_comm_models[[ddvr]][[paste0(sp1_name, sp2_name)]] <- list()
      
      # Store stage-dependent models
      for(sd_level_i in 1:length(sd_levels)) {
        sd_level <- sd_levels[sd_level_i]
        sd_level_str <- sd_levels_str[sd_level_i]
        
        sp1_eq_struct <- as.numeric(pop_eq_structs_df[
          pop_eq_structs_df$dd_vr == ddvr &
            pop_eq_structs_df$sp1_name == sp1_name &
            pop_eq_structs_df$sp2_name == sp2_name,
        ][1, c("sp1_j_norm", "sp1_a_norm")])
        sp2_eq_struct <- as.numeric(pop_eq_structs_df[
          pop_eq_structs_df$dd_vr == ddvr &
            pop_eq_structs_df$sp1_name == sp1_name &
            pop_eq_structs_df$sp2_name == sp2_name,
        ][1, c("sp2_j_norm", "sp2_a_norm")])
        
        # sp1_eq_struct <- as.numeric(pop_eq_avg_an[
        #   pop_eq_avg_an$dd_vr == ddvr &
        #     pop_eq_avg_an$sp_name == sp1_name,
        # ][1, c("mean_j_norm", "mean_a_norm")])
        # sp2_eq_struct <- as.numeric(pop_eq_avg_an[
        #   pop_eq_avg_an$dd_vr == ddvr &
        #     pop_eq_avg_an$sp_name == sp2_name,
        # ][1, c("mean_j_norm", "mean_a_norm")])
        
        sd_comm_model <-
          alter_stdep(
            comm_models[[ddvr]][[paste0(sp1_name, sp2_name)]],
            sd_level,
            ddvr,
            sp1_eq_struct,
            sp2_eq_struct
          )
        
        sd_comm_models[[ddvr]][[paste0(sp1_name, sp2_name)]][[sd_level_str]] <-
          sd_comm_model
        
        # Extract interaction coefficients
        kmat_section <- NULL
        if(ddvr == "sj") {
          kmat_section <- get_kmat_section(sd_comm_model, "s", 1)
        } else if(ddvr == "sa") {
          kmat_section <- get_kmat_section(sd_comm_model, "s", 2)
        } else {
          kmat_section <- get_kmat_section(sd_comm_model, ddvr, 1)
        }
        kmat_row <- c(kmat_section[1,], kmat_section[2,])
        kmat_row_temp <- kmat_row
        # Swap to make conspecific & heterospecific effects consistens
        kmat_row[6:7] <- kmat_row_temp[8:9]
        kmat_row[8:9] <- kmat_row_temp[6:7]
        kmat_row <- as.list(kmat_row)
        names(kmat_row) <- c(
          "sp1_cj", "sp1_ca", "sp1_hj", "sp1_ha", "sp1_intercept",
          "sp2_cj", "sp2_ca", "sp2_hj", "sp2_ha", "sp2_intercept"
        )
        kmat_row <- as.data.frame(kmat_row)
        kmat_row <- bind_cols(
          tibble(
            dd_vr = ddvr,
            sp1_name = sp1_name,
            sp2_name = sp2_name,
            st_dep_level = sd_level_str
          ),
          kmat_row
        )
        
        if(is.null(kmat_df)) {
          kmat_df <- kmat_row
        } else {
          kmat_df <- bind_rows(kmat_df, kmat_row)
        }
      }
    }
  }
}

kmat_df <- kmat_df %>%
  mutate(
    st_dep_level = factor(st_dep_level, sd_levels_str)
  )

write_csv(
  kmat_df,
  file.path("result_data","inter_coeffs.csv")
)
