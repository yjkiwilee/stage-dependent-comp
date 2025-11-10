# ===== File for preliminary model simulations for stochastic models, fast life history species =====

# ==== 1. Test run of the stochastic model =====

# Load required packages (Rage is not used here but NB the functions were written to generate MPMs that are compatible with Rage)
pacman::p_load("tidyverse", "Rage", "here", "scales", "RColorBrewer")

# Exported figure dimensions
fig_def_width <- 6
fig_def_height <- 4

# All paths are relative to the resilienceNERC master folder

# Load R file with deterministic model functions
source(here("StageDependentComp", "code", "comm_sims", "func", "sdcompr.R"))

# Number of life cycle stages
n_stages <- 2

# Density-independent vital rates of species 1 (~Annual species)
sp1_vrs <- list(
  s = c(1, 0.7, 0.1), # Seed, juvenile, adult survival
  g = c(0.9), # Maturation rate
  f = c(2)
)

# Density-independent vital rates of species 2 (~Annual species)
sp2_vrs <- list(
  s = c(1, 0.8, 0.1),
  g = c(0.9),
  f = c(2)
)

# Initial population structures for simulations
sp1_0 <- c(100, 0)
sp2_0 <- c(100, 0)

# Simulation length for the figures
sim_timestep <- 100

# Simulation length for the comparison of mean and standard deviation
sim_lat_timestep <- 1000
# Number of timesteps before the end to sample to calculate mean and sd
sim_lat_sample <- 700

# ====== 1.1. Stochastic stage-independent competition model ======

# K matrices describing stage-dependent competition for each of the vital rates
# Coffeicients are stage-independent
k_mat <- list(
  s = matrix(c(
    0, 0, 0, 0,
    0.001, 0.001, 0.0003, 0.0003, # Only juvenile survival is stage-dependent, with greater intraspecific competition
    0, 0, 0, 0,
    0, 0, 0, 0,
    0.0003, 0.0003, 0.001, 0.001,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  g = matrix(c( # No density-dependence on maturation
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  f = matrix(c( # Or fecundity
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE)
)

# Function defining the absolute disturbance to be applied to the population at each timestep
dstb_func <- function(t, sp1_pop, sp2_pop) {
  # Proportional disturbance
  sp1_prop_dstb <- -abs(rnorm(2, mean = 0, sd = 0.07))
  sp2_prop_dstb <- -abs(rnorm(2, mean = 0, sd = 0.07))
  
  # Calculate absolute disturbance
  sp1_abs_dstb <- sp1_pop * sp1_prop_dstb
  sp2_abs_dstb <- sp2_pop * sp2_prop_dstb
  
  # Return absolute disturbance as list
  dstb <- list(
    sp1 = sp1_abs_dstb,
    sp2 = sp2_abs_dstb
  )
  return(dstb)
}

# Function defining the stochastic density-indepdentend vital rates at each timestep
vrstoch_func <- function(t, sdcomp_model, sp1_pop, sp2_pop) {
  sp1_vrs <- sdcomp_model$sp1_vrs
  sp2_vrs <- sdcomp_model$sp2_vrs
  
  sp1_vrs$s <- sp1_vrs$s + rnorm(length(sp1_vrs$s), mean = 0, sd = 0.01)
  sp1_vrs$g <- sp1_vrs$g + rnorm(length(sp1_vrs$g), mean = 0, sd = 0.01)
  sp1_vrs$f <- sp1_vrs$f + rnorm(length(sp1_vrs$f), mean = 0, sd = 0.05)
  sp1_vrs$s <- pmax(pmin(sp1_vrs$s, 1), 0)
  sp1_vrs$g <- pmax(pmin(sp1_vrs$g, 1), 0)
  sp1_vrs$f <- pmax(sp1_vrs$f, 0)
  
  sp2_vrs$s <- sp2_vrs$s + rnorm(length(sp2_vrs$s), mean = 0, sd = 0.01)
  sp2_vrs$g <- sp2_vrs$g + rnorm(length(sp2_vrs$g), mean = 0, sd = 0.01)
  sp2_vrs$f <- sp2_vrs$f + rnorm(length(sp2_vrs$f), mean = 0, sd = 0.05)
  sp2_vrs$s <- pmax(pmin(sp2_vrs$s, 1), 0)
  sp2_vrs$g <- pmax(pmin(sp2_vrs$g, 1), 0)
  sp2_vrs$f <- pmax(sp2_vrs$f, 0)
  
  return(
    list(sp1_vrs = sp1_vrs, sp2_vrs = sp2_vrs)
  )
}

# Define model with these parameters
st_indep_m <- def_sdcomp_model(n_stages, sp1_vrs, sp2_vrs, k_mat)
# Set seed for reproducibility
set.seed(1)
# Calculate trajectory from initial population structures
st_indep_traj <- sdcomp_project_dstb(st_indep_m, sp1_0, sp2_0, timestep = sim_timestep, vrstoch_func = vrstoch_func)
# Convert trajectory data to long dataframe
st_indep_traj_long <- st_indep_traj %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")
# Plot trajectory
ggplot(st_indep_traj_long, aes(x = time, y = density, color = species, linetype = stage)) +
  geom_line() +
  theme_jun1() +
  scale_color_manual(breaks = c("sp1", "sp2"), values = brewer.pal(n=2,"Set2"), labels = c("Species 1", "Species 2")) +
  scale_linetype_manual(breaks = c("s1", "s2"), values = c("solid", "dashed"), labels = c("Small", "Large")) +
  labs(
    x = "Time",
    y = "Density",
    color = "Species",
    linetype = "Stage"
  )
# Save figure
ggsave(here("StageDependentComp", "figures", "prelim", "stoch_stindep_ann.png"), width = fig_def_width - 1, height = fig_def_height)

# Calculate adult / juevnile ratio & drop raw stage sizes
st_indep_ajr <- st_indep_traj %>%
  mutate(
    sp1 = sp1_s2 / sp1_s1,
    sp2 = sp2_s2 / sp2_s1
  ) %>%
  select(c(time, sp1, sp2))
# Convert to long dataframe
st_indep_ajr_long <- st_indep_ajr %>%
  pivot_longer(!time, names_to = "species", values_to = "aj_ratio")
# Plot data
ggplot(st_indep_ajr_long, aes(x = time, y = aj_ratio, color = species)) +
  geom_line() +
  theme_bw()
# Save figure
ggsave(here("StageDependentComp/figures/prelim/stoch_stindep_ann_ajr.png"), width = fig_def_width, height = fig_def_height)

# ====== 1.2. Stochastic stage-dependent competition model ======

# Alter elements of k_mat for stage-dependent competition while keeping the species-average per-capita competitive strength
# at equilibrium population structure constant
# Here, the juvenile survival is being set to be stage-dependent
st_dep_m <- alter_k_mat(st_indep_m,
                        vr_type = "s", # Alter survival component of k_mat
                        mat_row = 2,
                        mat_col = 1, # Altering row 2, column 1 of the survival k_mat
                        alt_val = 0.0007, # Set it so that the adults' effect on juvenile survival is now greater than that of juveniles
                        iter = 70, # Number of iterations to numerically determine the appropriate k_mat configuration
                        scal_fact_l = -100,
                        scal_fact_r = 100) # Minimum and maximum scaling factors to try in the first instance; the range should contain the desired scaling factor.
st_dep_m <- alter_k_mat(st_dep_m,
                        vr_type = "s", # Alter survival component of k_mat
                        mat_row = 2,
                        mat_col = 3, # Altering row 2, column 4 of the survival k_mat
                        alt_val = 0.00025, # Set it so that the adults' effect on juvenile survival is now greater than that of juveniles
                        iter = 70, # Number of iterations to numerically determine the appropriate k_mat configuration
                        scal_fact_l = -100,
                        scal_fact_r = 100) # Minimum and maximum scaling factors to try in the first instance; the range should contain the desired scaling factor.
# Do the same with juvenile survival of species 2
st_dep_m <- alter_k_mat(st_dep_m, vr_type = "s", mat_row = 5, mat_col = 1, alt_val = 0.00025, iter = 70, scal_fact_l = -100, scal_fact_r = 100)
st_dep_m <- alter_k_mat(st_dep_m, vr_type = "s", mat_row = 5, mat_col = 3, alt_val = 0.0007, iter = 70, scal_fact_l = -100, scal_fact_r = 100)

# Set seed for reproducibility
set.seed(1)
# Calculate trajectory from initial population structures
st_dep_traj <- sdcomp_project_dstb(st_dep_m, sp1_0, sp2_0, timestep = sim_timestep, dstb_func = dstb_func)
# Convert trajectory data to long dataframe
st_dep_traj_long <- st_dep_traj %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")
# Plot trajectory
ggplot(st_dep_traj_long, aes(x = time, y = density, color = species, linetype = stage)) +
  geom_line() +
  theme_bw()
ggsave(here("StageDependentComp/figures/prelim/stoch_stdep_ann.png"), width = fig_def_width, height = fig_def_height)

# Calculate adult / juevnile ratio & drop raw stage sizes
st_dep_ajr <- st_dep_traj %>%
  mutate(
    sp1 = sp1_s2 / sp1_s1,
    sp2 = sp2_s2 / sp2_s1
  ) %>%
  select(c(time, sp1, sp2))
# Convert to long dataframe
st_dep_ajr_long <- st_dep_ajr %>%
  pivot_longer(!time, names_to = "species", values_to = "aj_ratio")
# Plot data
ggplot(st_dep_ajr_long, aes(x = time, y = aj_ratio, color = species)) +
  geom_line() +
  theme_bw()
# Save figure
ggsave(here("StageDependentComp/figures/prelim/stoch_stdep_ann_ajr.png"), width = fig_def_width, height = fig_def_height)


# ===== 1.3. Compare dispersion in the stage-dependent vs. stage-independent models =====

# Simulate models for longer, setting seed for reproducibility
set.seed(1)
st_indep_traj <- sdcomp_project_dstb(st_indep_m, sp1_0, sp2_0, timestep = sim_lat_timestep, dstb_func = dstb_func)
set.seed(1)
st_dep_traj <- sdcomp_project_dstb(st_dep_m, sp1_0, sp2_0, timestep = sim_lat_timestep, dstb_func = dstb_func)

# Extract sample of population structure & convert to long format
st_indep_samp <- st_indep_traj[(nrow(st_indep_traj) - sim_lat_sample + 1):nrow(st_indep_traj),] %>%
  mutate(stage_dep = FALSE) %>%
  pivot_longer(!c(time, stage_dep), names_to = c("group"), values_to = "density")
st_dep_samp <- st_dep_traj[(nrow(st_dep_traj) - sim_lat_sample + 1):nrow(st_dep_traj),] %>%
  mutate(stage_dep = TRUE) %>%
  pivot_longer(!c(time, stage_dep), names_to = c("group"), values_to = "density")

# Merge datasets together
sims_samp <- rbind(st_indep_samp, st_dep_samp)

# Produce summary plot
ggplot(sims_samp, aes(x = group, y = density, color = stage_dep)) +
  geom_boxplot() +
  geom_point(size = 0.3, alpha = 0.1, position = position_jitterdodge(jitter.width = 0.2)) +
  theme_bw()
# Save figure
ggsave(here("StageDependentComp/figures/prelim/stoch_summ.png"), width = fig_def_width, height = fig_def_height)
