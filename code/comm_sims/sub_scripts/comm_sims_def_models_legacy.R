################################################################################
#' Script for defining the virtual species
#'
#' Written by Young Jun Lee
#' Mar 2025

# ======= Define virtual species =======

require("here")

source(here("StageDependentComp","code","comm_sims","sub_scripts","comm_sims_init_setup.R"))


# Constant parameters
## Number of life cycle stages
n_stages <- 2
## Standard deviation of stochasticity to be applied to weighted sum
vr_sd <- 0.2
## Intra/interspecific interaction coefficient
w_intra <- 4.5
w_inter <- 1
## Target equilibrium population size
n_eq <- 500
## Density-dependence sign coefficients
dd_sign <- list(
  sj = -1,
  sa = -1,
  g = -1,
  r = 1,
  f = -1
)


# Density-independent vital rates of each species
spp_vrs <- lapply(seq(0, 1, length.out = 5), function(x) {
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
  # Calculate fertility so that lambda is equal to target
  vrs$f = c(
    fert_for_lambda(
      vrs$s[1], vrs$s[2],
      vrs$g[1],
      vrs$r[1],
      target_lambda
    )
  )
  # Return vital rates
  return(vrs)
})

names(spp_vrs) <- c("F2", "F1", "M", "S1", "S2")

spp_vrs_list <- lapply(names(spp_vrs), function(sp_name) {
  vrs <- spp_vrs[[sp_name]]
  sp_mpm <- vrs_to_mpm(n_stages, vrs)
  
  spp_vrs_row <- tibble(
    sp_name = sp_name,
    sj = vrs$s[1],
    sa = vrs$s[2],
    g = vrs$g,
    r = vrs$r,
    f = vrs$f,
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

write_csv(
  bind_rows(spp_vrs_list),
  here("StageDependentComp","result_data","wp3","spp_vrs.csv")
)

# ======= Define stage-independent community models =======


# Calculate population size coefficients
pop_coeffs_df <- tibble(
  dd_vr = character(),
  sp1_name = character(),
  sp2_name = character(),
  sp1_pop_coeff = numeric(),
  sp2_pop_coeff = numeric()
)

# Check if population coefficients have already been calculated
if(file.exists(here("StageDependentComp","result_data","wp3","pop_coeffs.csv"))) {
  pop_coeffs_df <- read_csv(
    here("StageDependentComp","result_data","wp3","pop_coeffs.csv")
  )
} else {
  set.seed(1)
  for(dd_vr_name in c("sj", "sa", "g", "r", "f")) {
    vr_sds <- list(
      s = c(0, 0), # Small & large survival
      g = c(0), # Maturation rate
      r = c(0), # Retrogression rate
      f = c(0) # Fecundity
    )
    
    if(dd_vr_name =="sj") {
      vr_sds$s[1] <- vr_sd
    } else if(dd_vr_name == "sa") {
      vr_sds$s[2] <- vr_sd
    } else {
      vr_sds[[dd_vr_name]][1] <- vr_sd
    }
    
    for(sp1_i in 1:length(names(spp_vrs))) {
      sp1_name <- names(spp_vrs)[sp1_i]
      for(sp2_i in sp1_i:length(names(spp_vrs))) {
        sp2_name <- names(spp_vrs)[sp2_i]
        cat(paste("Processing",dd_vr_name,sp1_name,sp2_name,"\n"))
        cat(paste0(format(Sys.time(), "%c"), "\n"))
        
        pop_coeffs_temp <-
          calc_pop_coeff_num(
            spp_vrs[[sp1_name]],
            spp_vrs[[sp2_name]],
            dd_vr_name,
            dd_sign[[dd_vr_name]],
            n_eq,
            w_intra,
            w_inter,
            vr_sds
          )
        
        pop_coeffs_row <- NULL
        
        # Average two population coefficients if two populations are of the same species
        if(sp1_name == sp2_name) {
          pop_coeffs_row <- tibble(
            dd_vr = dd_vr_name,
            sp1_name = sp1_name,
            sp2_name = sp2_name,
            sp1_pop_coeff = mean(pop_coeffs_temp),
            sp2_pop_coeff = mean(pop_coeffs_temp)
          )
        } else {
          pop_coeffs_row <- tibble(
            dd_vr = dd_vr_name,
            sp1_name = sp1_name,
            sp2_name = sp2_name,
            sp1_pop_coeff = pop_coeffs_temp[1],
            sp2_pop_coeff = pop_coeffs_temp[2]
          )
        }
        
        pop_coeffs_df <- bind_rows(pop_coeffs_df, pop_coeffs_row)
        
        write_csv(
          pop_coeffs_df,
          here("StageDependentComp","result_data","wp3","pop_coeffs.csv")
        )
      }
    }
  }
}

# Test correlation between population coefficient and competitor identity
# sp_i <- c(
#   "F2" = 1,
#   "F1" = 2,
#   "M" = 3,
#   "S1" = 4,
#   "S2" = 5
# )
# pop_coeffs_test <- pop_coeffs_df %>%
#   mutate(
#     sp1_i = sp_i[sp1_name],
#     sp2_i = sp_i[sp2_name]
#   )
# pop_coeffs_test_long <- tibble(
#   dd_vr = rep(pop_coeffs_test$dd_vr, times = 2),
#   foc_sp = c(pop_coeffs_test$sp1_name, pop_coeffs_test$sp2_name),
#   foc_sp_i = c(pop_coeffs_test$sp1_i, pop_coeffs_test$sp2_i),
#   comp_sp_i = c(pop_coeffs_test$sp2_i, pop_coeffs_test$sp1_i),
#   comp_sp = c(pop_coeffs_test$sp2_name, pop_coeffs_test$sp1_name),
#   foc_pop_coeff = c(pop_coeffs_test$sp1_pop_coeff, pop_coeffs_test$sp2_pop_coeff)
# )
# test_lm <- lm(
#   foc_pop_coeff ~ dd_vr + foc_sp_i:dd_vr + comp_sp_i:dd_vr,
#   pop_coeffs_test_long
# )
# summary(test_lm)
# test_lm_inter <- lm(
#   foc_pop_coeff ~ dd_vr + foc_sp_i:dd_vr + comp_sp_i:dd_vr + foc_sp_i:comp_sp_i:dd_vr,
#   pop_coeffs_test_long
# )
# summary(test_lm_inter)
# # Population coefficients are not significantly correlated with competitor identity
# # Therefore: average population coefficients
# 
# pop_coeffs_avg <- pop_coeffs_test_long %>%
#   group_by(dd_vr, foc_sp) %>%
#   summarise(
#     avg_pop_coeff = mean(foc_pop_coeff)
#   )
# 
# write_csv(
#   pop_coeffs_avg,
#   here("StageDependentComp","result_data","wp3","pop_coeffs_avg.csv")
# )

# Calculate equilibrium population structures

eq_struct_tmin <- 500
eq_struct_tmax <- 1500
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


if(file.exists(here("StageDependentComp","result_data","wp3","pop_eq_structs.csv"))) {
  pop_eq_structs_df <- read_csv(
    here("StageDependentComp","result_data","wp3","pop_eq_structs.csv")
  )
} else {
  for(dd_vr_name in c("sj", "sa", "g", "r", "f")) {
    for(sp1_i in 1:length(names(spp_vrs))) {
      sp1_name <- names(spp_vrs)[sp1_i]
      for(sp2_i in sp1_i:length(names(spp_vrs))) {
        sp2_name <- names(spp_vrs)[sp2_i]
        cat(paste("Processing",dd_vr_name,sp1_name,sp2_name,"\n"))
        cat(paste0(format(Sys.time(), "%c"), "\n"))
        
        vr_sds <- list(
          s = c(0, 0), # Small & large survival
          g = c(0), # Maturation rate
          r = c(0), # Retrogression rate
          f = c(0) # Fecundity
        )
        
        if(dd_vr_name =="sj") {
          vr_sds$s[1] <- vr_sd
        } else if(dd_vr_name == "sa") {
          vr_sds$s[2] <- vr_sd
        } else {
          vr_sds[[dd_vr_name]][1] <- vr_sd
        }
        
        # Build community model
        temp_comm_model <- build_base_model(
          spp_vrs[[sp1_name]],
          spp_vrs[[sp2_name]],
          dd_vr_name,
          dd_sign[[dd_vr_name]],
          w_intra,
          w_inter,
          c(
            pop_coeffs_df[pop_coeffs_df$dd_vr == dd_vr_name &
                             pop_coeffs_df$sp1_name == sp1_name &
                            pop_coeffs_df$sp2_name == sp2_name,]$sp1_pop_coeff,
            pop_coeffs_df[pop_coeffs_df$dd_vr == dd_vr_name &
                            pop_coeffs_df$sp1_name == sp1_name &
                            pop_coeffs_df$sp2_name == sp2_name,]$sp2_pop_coeff
          ),
          vr_sds
        )
        
        eq_sp1_j <- 0
        eq_sp1_a <- 0
        eq_sp2_j <- 0
        eq_sp2_a <- 0
        
        # Repeat simulations nrep times
        for(i in 1:eq_struct_nrep) {
          temp_comm_proj <- sdcomp_project(
            temp_comm_model,
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
          here("StageDependentComp","result_data","wp3","pop_eq_structs.csv")
        )
      }
    }
  }
}

# # Test correlation between population structure and competitor identity
# sp_i <- c(
#   "F2" = 1,
#   "F1" = 2,
#   "M" = 3,
#   "S1" = 4,
#   "S2" = 5
# )
# eq_struct_test <- pop_eq_structs_df %>%
#   mutate(
#     sp1_i = sp_i[sp1_name],
#     sp2_i = sp_i[sp2_name]
#   )
# eq_struct_test_long <- tibble(
#   dd_vr = rep(eq_struct_test$dd_vr, times = 2),
#   foc_sp = c(eq_struct_test$sp1_name, eq_struct_test$sp2_name),
#   foc_sp_i = c(eq_struct_test$sp1_i, eq_struct_test$sp2_i),
#   comp_sp_i = c(eq_struct_test$sp2_i, eq_struct_test$sp1_i),
#   comp_sp = c(eq_struct_test$sp2_name, eq_struct_test$sp1_name),
#   foc_j_norm = c(eq_struct_test$sp1_j_norm, eq_struct_test$sp2_j_norm),
#   foc_n = c(eq_struct_test$sp1_n, eq_struct_test$sp2_n)
# )
# test_lm <- lm(
#   foc_j_norm ~ dd_vr + foc_sp_i:dd_vr + comp_sp_i:dd_vr,
#   eq_struct_test_long
# )
# summary(test_lm)
# test_lm_inter <- lm(
#   foc_j_norm ~ dd_vr + foc_sp_i:dd_vr + comp_sp_i:dd_vr + foc_sp_i:comp_sp_i:dd_vr,
#   eq_struct_test_long
# )
# summary(test_lm_inter)
# # Competitor identity doesn't affect population structure;
# # Therefore, average population structure across simulations
# 
# # Test correlation between population size and competitor identity
# test_lm <- lm(
#   foc_n ~ dd_vr + foc_sp_i:dd_vr + comp_sp_i:dd_vr,
#   eq_struct_test_long
# )
# summary(test_lm)
# test_lm_inter <- lm(
#   foc_n ~ dd_vr + foc_sp_i:dd_vr + comp_sp_i:dd_vr + foc_sp_i:comp_sp_i:dd_vr,
#   eq_struct_test_long
# )
# summary(test_lm_inter)
# # Unfortunately, this procedure didn't effectively control population size...

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
  group_by(dd_vr, sp_name) %>%
  summarise(
    mean_j = mean(j),
    mean_a = mean(a),
    mean_j_norm = mean(j_norm),
    mean_a_norm = mean(a_norm),
    mean_n = mean(n)
  )

write_csv(
  pop_eq_avg,
  here("StageDependentComp","result_data","wp3","pop_eq_structs_avg.csv")
)


# List of lists to store base community models
base_comm_models <- list(
  sj = NULL,
  sa = NULL,
  g = NULL,
  r = NULL,
  f = NULL
)

# DF for storing interaction coefficients
kmat_df <- NULL
# Levels of stage dependence
sd_levels <- seq(-1, 1, 0.2)

# For each density-dependent vital rate:
for(ddvr in names(base_comm_models)) {
  base_comm_models[[ddvr]] <- list()
  # Iterate through pairs of virtual species
  for(sp1_i in 1:length(names(spp_vrs))) {
    for(sp2_i in sp1_i:length(names(spp_vrs))) {
      sp1_name <- names(spp_vrs)[sp1_i]
      sp2_name <- names(spp_vrs)[sp2_i]
      
      vr_sds <- list(
        s = c(0, 0), # Small & large survival
        g = c(0), # Maturation rate
        r = c(0), # Retrogression rate
        f = c(0) # Fecundity
      )
      
      if(ddvr == "sj") {
        vr_sds$s[1] <- vr_sd
      } else if(ddvr == "sa") {
        vr_sds$s[2] <- vr_sd
      } else {
        vr_sds[[ddvr]][1] <- vr_sd
      }
      
      # Make community model
      base_comm_models[[ddvr]][[paste0(sp1_name, sp2_name)]] <-
        build_base_model(
          spp_vrs[[sp1_name]],
          spp_vrs[[sp2_name]],
          ddvr,
          dd_sign[[ddvr]],
          w_intra,
          w_inter,
          c(
            # pop_coeffs_df[pop_coeffs_df$dd_vr == ddvr &
            #                 pop_coeffs_df$sp1_name == sp1_name &
            #                 pop_coeffs_df$sp2_name == sp2_name,]$sp1_pop_coeff,
            # pop_coeffs_df[pop_coeffs_df$dd_vr == ddvr &
            #                 pop_coeffs_df$sp1_name == sp1_name &
            #                 pop_coeffs_df$sp2_name == sp2_name,]$sp2_pop_coeff
            0.00001,
            0.00001
          ),
          vr_sds
        )
      
      # Store stage-dependent coefficients
      for(sd_level in sd_levels) {
        sp1_j_norm <- pop_eq_avg[
          pop_eq_avg$dd_vr == ddvr &
            pop_eq_avg$sp_name == sp1_name,
        ]$mean_j_norm
        sp2_j_norm <- pop_eq_avg[
          pop_eq_avg$dd_vr == ddvr &
            pop_eq_avg$sp_name == sp2_name,
        ]$mean_j_norm
        
        sp1_eq_struct <- c(sp1_j_norm, 1 - sp1_j_norm)
        sp2_eq_struct <- c(sp2_j_norm, 1 - sp2_j_norm)
        
        sd_comm_model <- alter_stdep(
          base_comm_models[[ddvr]][[paste0(sp1_name, sp2_name)]],
          sd_level,
          ddvr,
          sp1_eq_struct,
          sp2_eq_struct
        )
        
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
            st_dep_level = sd_level
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

# write_csv(
#   kmat_df,
#   here("StageDependentComp","result_data","wp3","inter_coeffs.csv")
# )


stoch_rnorms <- matrix(
  rnorm(2000 * 2 * (2+1+1+1)), ncol = 2 * (2+1+1+1)
)
stoch_rnorms <- rbind(stoch_rnorms, rep(NA, 2 * (2+1+1+1)))

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

# Organise stochasticity terms for each species
stoch_rnorms_sp1 <- lapply(1:2000, function(t) {
  list(
    s = c(stoch_rnorms[t, c(1,2)]), # Small & large survival
    g = stoch_rnorms[t, 3], # Maturation rate
    r = stoch_rnorms[t, 4], # Retrogression rate
    f = stoch_rnorms[t, 5] # Fecundity
  )
})

stoch_rnorms_sp2 <- lapply(1:2000, function(t) {
  list(
    s = c(stoch_rnorms[t, c(6,7)]), # Small & large survival
    g = stoch_rnorms[t, 8], # Maturation rate
    r = stoch_rnorms[t, 9], # Retrogression rate
    f = stoch_rnorms[t, 10] # Fecundity
  )
})

# comm_project <- sdcomp_project_vr(
#   test_model,
#   c(250, 250),
#   c(250, 250),
#   timestep = 2000
# )

calc_f_for_dd <- function(base_vrs, target_dn) {
  sj <- logistic(logit(base_vrs$s[1]) - target_dn)
  sa <- base_vrs$s[2]
  g <- base_vrs$g
  r <- base_vrs$r
  
  return(fert_for_lambda(sj, sa, g, r, 1))
}

spp_vrs_test <- spp_vrs

spp_vrs_test$F2$f <- calc_f_for_dd(
  spp_vrs$F2,
  1000 * 0.00045 + 1000 * 0.0001
)

spp_vrs_test$S2$f <- calc_f_for_dd(
  spp_vrs$S2,
  1000 * 0.00045 + 1000 * 0.0001
)

vr_sds <- list(
  s = c(0.2, 0), # Small & large survival
  g = c(0), # Maturation rate
  r = c(0), # Retrogression rate
  f = c(0) # Fecundity
)

test_base_model <- def_sdcomp_model(
  2,
  spp_vrs_test$F2,
  spp_vrs_test$S2,
  list(
    s = matrix(c(
      -0.00045, -0.00045, -0.0001, -0.0001,
      0, 0, 0, 0,
      -0.0001, -0.0001, -0.00045, -0.00045,
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
  ),
  vr_sds, vr_sds
)

comm_project <- sdcomp_project_vr(
  test_base_model,
  c(500, 500),
  c(500, 500),
  timestep = 2000,
  stoch_rnorms_sp1,
  stoch_rnorms_sp2
)

comm_project <- comm_project %>%
  bind_cols(stoch_rnorms)

comm_project <- comm_project %>%
  mutate(
    sp1_n = sp1_s1 + sp1_s2,
    sp2_n = sp2_s1 + sp2_s2,
    sp1_struct = sp1_s1 / (sp1_s1 + sp1_s2),
    sp2_struct = sp2_s1 / (sp2_s1 + sp2_s2),
    vr_sp1_s1_trans = logit(vr_sp1_s1),
    vr_sp2_s1_trans = logit(vr_sp2_s1),
    sp1_lambda = lambda_from_vr(
      vr_sp1_s1,
      vr_sp1_s2,
      vr_sp1_g1,
      vr_sp1_r1,
      vr_sp1_f1
    ),
    sp2_lambda = lambda_from_vr(
      vr_sp2_s1,
      vr_sp2_s2,
      vr_sp2_g1,
      vr_sp2_r1,
      vr_sp2_f1
    ),
  )

ggplot(comm_project[1:2000,], aes(x = time)) +
  geom_line(aes(y = sp1_n), color = "red") +
  geom_line(aes(y = sp2_n), color = "blue")

sp1_pop_vec <- matrix(c(
  mean(comm_project[500:2000,]$sp1_s1),
  mean(comm_project[500:2000,]$sp1_s2)
), ncol = 1)

sp2_pop_vec <- matrix(c(
  mean(comm_project[500:2000,]$sp2_s1),
  mean(comm_project[500:2000,]$sp2_s2)
), ncol = 1)

summary(lm(
  vr_sp1_s1 ~ sp1_n + sp2_n,
  comm_project
))

summary(lm(
  vr_sp1_s1 ~ sp1_n + sp2_n + sp1_struct + sp2_struct,
  comm_project
))

summary(lm(
  vr_sp1_s1 ~ sp1_n + sp2_n + sp1_struct + sp2_struct +
    stoch_sj_1,
  comm_project
))

summary(lm(
  vr_sp2_s1_trans ~ sp1_n + sp2_n,
  comm_project
))

summary(lm(
  vr_sp2_s1_trans ~ sp1_n + sp2_n + sp1_struct + sp2_struct,
  comm_project
))

summary(lm(
  vr_sp2_s1 ~ sp1_n + sp2_n + sp1_struct + sp2_struct +
    stoch_sj_2,
  comm_project
))

test_sd_model <- alter_stdep(
  test_base_model,
  1,
  "sj",
  sp1_pop_vec / sum(sp1_pop_vec),
  sp2_pop_vec / sum(sp2_pop_vec)
)

comm_project <- sdcomp_project_vr(
  test_sd_model,
  c(500, 500),
  c(500, 500),
  timestep = 2000,
  stoch_rnorms_sp1,
  stoch_rnorms_sp2
)

comm_project <- comm_project %>%
  bind_cols(stoch_rnorms)

comm_project <- comm_project %>%
  mutate(
    sp1_n = sp1_s1 + sp1_s2,
    sp2_n = sp2_s1 + sp2_s2,
    sp1_struct = sp1_s1 / (sp1_s1 + sp1_s2),
    sp2_struct = sp2_s1 / (sp2_s1 + sp2_s2),
    vr_sp1_s1_trans = logit(vr_sp1_s1),
    vr_sp2_s1_trans = logit(vr_sp2_s1),
    sp1_lambda = lambda_from_vr(
      vr_sp1_s1,
      vr_sp1_s2,
      vr_sp1_g1,
      vr_sp1_r1,
      vr_sp1_f1
    ),
    sp2_lambda = lambda_from_vr(
      vr_sp2_s1,
      vr_sp2_s2,
      vr_sp2_g1,
      vr_sp2_r1,
      vr_sp2_f1
    ),
  )

ggplot(comm_project[1:2000,], aes(x = time)) +
  geom_line(aes(y = sp1_n), color = "red") +
  geom_line(aes(y = sp2_n), color = "blue")

summary(lm(
  vr_sp1_s1 ~ sp1_n + sp2_n,
  comm_project
))

summary(lm(
  sp1_lambda ~ sp1_n + sp2_n + sp1_struct + sp2_struct,
  comm_project
))

summary(lm(
  vr_sp1_s1 ~ sp1_n + sp2_n + sp1_struct + sp2_struct +
    stoch_sj_1,
  comm_project
))

summary(lm(
  vr_sp2_s1 ~ sp1_n + sp2_n,
  comm_project
))

summary(lm(
  sp2_lambda ~ sp1_n + sp2_n + sp1_struct + sp2_struct,
  comm_project
))

summary(lm(
  vr_sp2_s1 ~ sp1_n + sp2_n + sp1_struct + sp2_struct +
    stoch_sj_2,
  comm_project
))
