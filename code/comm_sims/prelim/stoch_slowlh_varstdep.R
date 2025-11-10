# ====== Testing stochastic models across a range of stage-dependence in ~annual species ======

# ====== 1. Define model parameters & setup ======

# Load required packages (Rage is not used here but NB the functions were written to generate MPMs that are compatible with Rage)
pacman::p_load("tidyverse", "Rage", "here", "ggpubr")

# Exported figure dimensions
fig_def_width <- 6
fig_def_height <- 4

# All paths are relative to the resilienceNERC master folder

# Load R file with model functions
source(here("StageDependentComp/code/sdcompr/sdcompr.R"))
# Load R file with custom themes
source(here("StageDependentComp/code/themes/custom_themes.R"))

# Number of life cycle stages
n_stages <- 2

# Density-independent vital rates of species 1 (Fast lh species)
sp1_vrs <- list(
  s = c(1, 0.2, 0.9), # Seed, juvenile, adult survival
  g = c(0.6), # Maturation rate
  f = c(2)
)

# Density-independent vital rates of species 2 (Fast lh species)
sp2_vrs <- list(
  s = c(1, 0.2, 0.9), # Seed, juvenile, adult survival
  g = c(0.6), # Maturation rate
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

# Original strength of intraspecific & interspecific competition
intrasp_comp <- 0.001
intersp_comp <- 0.0003

# K matrices describing stage-dependent competition for each of the vital rates
# Coffeicients are stage-independent
k_mat <- list(
  s = matrix(c(
    0, 0, 0, 0,
    intrasp_comp, intrasp_comp, intersp_comp, intersp_comp, # Only juvenile survival is stage-dependent, with greater intraspecific competition
    0, 0, 0, 0,
    0, 0, 0, 0,
    intersp_comp, intersp_comp, intrasp_comp, intrasp_comp,
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
  sp1_prop_dstb <- -abs(rnorm(2, mean = 0, sd = 0.01))
  sp2_prop_dstb <- -abs(rnorm(2, mean = 0, sd = 0.01))
  
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

# ====== 2. Calculate models with varying levels of stage-dependence ======

# Range of stage-dependence to impose
sds_scale <- seq(0, 1, 0.1)
intrasp_sds <- sds_scale * intrasp_comp # On intraspecific competition
intersp_sds <- sds_scale * intersp_comp # On interspecific competition

# List for storing the resulting models
sd_models <- list()

# Calculate model for each level of stage-dependence
# This will take a long time!
for(i in 1:length(sds_scale)) {
  cat(paste("Calculating model for stage-dependence of", sds_scale[i]))
  
  # New values for competitive effect of intraspecific and interspecific juveniles
  intra_jcomp <- intrasp_comp - intrasp_sds[i]
  inter_jcomp <- intersp_comp - intersp_sds[i]
  
  # Calculate stage-dependent model
  temp_model <- alter_k_mat_mult(st_indep_m, vr_type = "s",
                                 mat_rows = c(2, 2, 5, 5),
                                 mat_cols = c(1, 3, 1, 3),
                                 alt_vals = c(intra_jcomp, inter_jcomp, inter_jcomp, intra_jcomp),
                                 iter = 20, scal_fact_l = -100, scal_fact_r = 100, vrstoch_func = vrstoch_func)
  
  # Save model
  sd_models[[i]] <- temp_model
  
  cat(" Done! 🦊\n")
}

# ====== 3. Calculate trajectories for each model ======

# Dataframe for storing the resulting trajectories
sd_models_traj <- data.frame(
  st_dep = numeric(),
  intra_jcomp = numeric(),
  inter_jcomp = numeric(),
  time = numeric(),
  sp1_s1 = numeric(),
  sp1_s2 = numeric(),
  sp2_s1 = numeric(),
  sp2_s2 = numeric()
)

# Run models
for(i in 1:length(sd_models)) {
  # Calculate trajectory
  temp_traj <- sdcomp_project_dstb(sd_models[[i]], sp1_0, sp2_0, timestep = sim_lat_timestep, vrstoch_func = vrstoch_func) %>%
    mutate(st_dep = sds_scale[i],
           intra_jcomp = intrasp_comp - intrasp_sds[i],
           inter_jcomp = intersp_comp - intersp_sds[i]) # Save degree of stage-dependence
  
  # Merge calculated trajectory into the dataframe
  sd_models_traj <- rbind(sd_models_traj, temp_traj)
}

# ====== 4. Plot means ======

# Extract standard deviations for each level of stage-dependence from the data
sd_models_slh_m <- sd_models_traj %>%
  filter(time >= sim_lat_timestep - sim_lat_sample, time <= sim_lat_timestep + 1) %>% # Extract latter part of the simulation
  group_by(st_dep) %>%
  summarise_all(mean)

# Plot standard deviations for each stage
sp1_s1_mplot <- ggline(sd_models_slh_m, x = "st_dep", y = "sp1_s1", ylab = "mean_sp1_s1")
sp1_s2_mplot <- ggline(sd_models_slh_m, x = "st_dep", y = "sp1_s2", ylab = "mean_sp1_s2")
sp2_s1_mplot <- ggline(sd_models_slh_m, x = "st_dep", y = "sp2_s1", ylab = "mean_sp2_s1")
sp2_s2_mplot <- ggline(sd_models_slh_m, x = "st_dep", y = "sp2_s2", ylab = "mean_sp2_s2")

# Arrange plots
ggarrange(sp1_s1_mplot, sp1_s2_mplot, sp2_s1_mplot, sp2_s2_mplot, ncol = 2, nrow = 2)
ggsave(here("StageDependentComp/figures/prelim/stoch_slowlh_varstdep_mean_mult.png"), width = 8, height = 5)

# ====== 5. Plot standard deviations ======

# Extract standard deviations for each level of stage-dependence from the data
sd_models_slh_sd <- sd_models_traj %>%
  filter(time >= sim_lat_timestep - sim_lat_sample, time <= sim_lat_timestep + 1) %>% # Extract latter part of the simulation
  group_by(st_dep) %>%
  summarise_all(sd)

# Plot standard deviations for each stage
sp1_s1_sdplot <- ggline(sd_models_slh_sd, x = "st_dep", y = "sp1_s1", ylab = "sd_sp1_s1")
sp1_s2_sdplot <- ggline(sd_models_slh_sd, x = "st_dep", y = "sp1_s2", ylab = "sd_sp1_s2")
sp2_s1_sdplot <- ggline(sd_models_slh_sd, x = "st_dep", y = "sp2_s1", ylab = "sd_sp2_s1")
sp2_s2_sdplot <- ggline(sd_models_slh_sd, x = "st_dep", y = "sp2_s2", ylab = "sd_sp2_s2")

# Arrange plots
ggarrange(sp1_s1_sdplot, sp1_s2_sdplot, sp2_s1_sdplot, sp2_s2_sdplot, ncol = 2, nrow = 2)
ggsave(here("StageDependentComp/figures/prelim/stoch_slowlh_varstdep_sd_mult.png"), width = 8, height = 5)

# ====== B. Repeat analysis but focusing on the range that doesn't result in cyclic dynamics ======

# ====== B1. Calculate models with varying levels of stage-dependence ======

# Range of stage-dependence to impose
sds_scale <- seq(0, 0.4, 0.01)
intrasp_sds <- sds_scale * intrasp_comp # On intraspecific competition
intersp_sds <- sds_scale * intersp_comp # On interspecific competition

# List for storing the resulting models
sd_models_slh_zoom <- list()

# Calculate model for each level of stage-dependence
# This will take a long time!
for(i in 1:length(sds_scale)) {
  cat(paste("Calculating model for stage-dependence of", sds_scale[i]))
  
  # New values for competitive effect of intraspecific and interspecific juveniles
  intra_jcomp <- intrasp_comp - intrasp_sds[i]
  inter_jcomp <- intersp_comp - intersp_sds[i]
  
  # Calculate stage-dependent model
  temp_model <- alter_k_mat_mult(st_indep_m, vr_type = "s",
                                 mat_rows = c(2, 2, 5, 5),
                                 mat_cols = c(1, 3, 1, 3),
                                 alt_vals = c(intra_jcomp, inter_jcomp, inter_jcomp, intra_jcomp),
                                 iter = 20, scal_fact_l = -100, scal_fact_r = 100, vrstoch_func = vrstoch_func)
  
  # Save model
  sd_models_slh_zoom[[i]] <- temp_model
  
  cat(" Done! 🦊\n")
}

# ====== B2. Calculate trajectories for each model ======

# Dataframe for storing the resulting trajectories
sd_models_traj <- data.frame(
  st_dep = numeric(),
  intra_jcomp = numeric(),
  inter_jcomp = numeric(),
  time = numeric(),
  sp1_s1 = numeric(),
  sp1_s2 = numeric(),
  sp2_s1 = numeric(),
  sp2_s2 = numeric()
)

# Run models
for(i in 1:length(sd_models_slh_zoom)) {
  # Calculate trajectory
  temp_traj <- sdcomp_project_dstb(sd_models_slh_zoom[[i]], sp1_0, sp2_0, timestep = sim_lat_timestep, vrstoch_func = vrstoch_func) %>%
    mutate(st_dep = sds_scale[i]) # Save degree of stage-dependence
  
  # Merge calculated trajectory into the dataframe
  sd_models_traj <- rbind(sd_models_traj, temp_traj)
}

# ====== B3. Plot means ======

# Extract standard deviations for each level of stage-dependence from the data
sd_models_slh_m <- sd_models_traj %>%
  filter(time >= sim_lat_timestep - sim_lat_sample, time <= sim_lat_timestep + 1) %>% # Extract latter part of the simulation
  group_by(st_dep) %>%
  summarise_all(mean)
sd_models_slh_m_long <- sd_models_slh_m %>%
  pivot_longer(
    !c(st_dep, time),
    names_to = c("species", "stage"),
    names_sep = "_",
    values_to = "avg"
  ) %>%
  group_by(species)

# Plot means for each stage
species_labs <- c("Species 1", "Species 2")
names(species_labs) <- c("sp1", "sp2")
stage_labs <- c("Juvenile", "Adult")
names(stage_labs) <- c("s1", "s2")

ggplot(sd_models_slh_m_long, aes(x = st_dep, y = avg)) +
  geom_point() +
  theme_jun2() +
  labs(
    x = "Degree of stage-dependence in competition",
    y = "Long-term mean of density"
  ) +
  facet_grid(stage ~ species, scales = "free_y", labeller = labeller(species = species_labs, stage = stage_labs))

ggsave(here("StageDependentComp/figures/prelim/stoch_slowlh_varstdep_mean_zoom_mult.png"), width = 8, height = 5)

# ====== B4. Plot standard deviations ======

# Extract standard deviations for each level of stage-dependence from the data
sd_models_slh_sd <- sd_models_traj %>%
  filter(time >= sim_lat_timestep - sim_lat_sample, time <= sim_lat_timestep + 1) %>% # Extract latter part of the simulation
  group_by(st_dep) %>%
  summarise_all(sd)

sd_models_slh_sd_long <- sd_models_slh_sd %>%
  pivot_longer(
    !c(st_dep, time),
    names_to = c("species", "stage"),
    names_sep = "_",
    values_to = "sd"
  ) %>%
  group_by(species)

# Plot standard deviations for each stage
species_labs <- c("Species 1", "Species 2")
names(species_labs) <- c("sp1", "sp2")
stage_labs <- c("Juvenile", "Adult")
names(stage_labs) <- c("s1", "s2")

# y limit setter
ylim_setter <- tibble(
  st_dep = c(0.2, 0.2, 0.2, 0.2),
  sd = c(16, 53, 5, 25),
  species = c("sp1", "sp1", "sp1", "sp1"),
  stage = c("s1", "s1", "s2", "s2")
)

ggplot(sd_models_slh_sd_long, aes(x = st_dep, y = sd, color = species)) +
  geom_point() +
  geom_point(data = ylim_setter, alpha = 0) +
  theme_jun2() +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = "Degree of stage-dependence in competition",
    y = "Long-term standard deviation of density"
  ) +
  guides(color = "none") +
  facet_grid(stage ~ species, scales = "free_y", labeller = labeller(species = species_labs, stage = stage_labs))

# Save plot
ggsave(here("StageDependentComp/figures/prelim/stoch_slowlh_varstdep_sd_zoom_mult.png"), width = 5, height = 4)

# ====== B5. Test correlation between stage dependence and standard deviation =====

cor.test(sd_models_slh_sd$st_dep, sd_models_slh_sd$sp1_s1, method = "pearson")
cor.test(sd_models_slh_sd$st_dep, sd_models_slh_sd$sp1_s2, method = "pearson")
cor.test(sd_models_slh_sd$st_dep, sd_models_slh_sd$sp2_s1, method = "pearson")
cor.test(sd_models_slh_sd$st_dep, sd_models_slh_sd$sp2_s2, method = "pearson")
cor.test(sd_models_slh_sd$sp1_s2, sd_models_slh_sd$sp2_s2, method = "pearson")
cor.test(sd_models_slh_sd$sp1_s1, sd_models_slh_sd$sp2_s1, method = "pearson")

# ====== C. Inspect k_mat elements ======

# Variable for storing the k_mat elements for each level of stage-dependence as a matrix
sd_kmat_elems <- NULL

# Extract relevant rows in the k_mat and insert into the matrix (Only those that weren't directly manipulated)
for(i in 1:length(sd_models_slh_zoom)) {
  # Get row for sp1 juvenile survival
  sp1_sj_row <- sd_models_slh_zoom[[i]]$k_mat$s[2,c(2,4)]
  # Get row for sp2 juvenile survival
  sp2_sj_row <- sd_models_slh_zoom[[i]]$k_mat$s[5,c(2,4)]
  # Merge into single vector
  sj_rows <- c(sp1_sj_row, sp2_sj_row)
  # Insert into matrix
  if(is.null(sd_kmat_elems)) {
    sd_kmat_elems <- matrix(sj_rows, ncol = length(sj_rows))
  } else {
    sd_kmat_elems <- rbind(sd_kmat_elems, matrix(sj_rows, ncol = length(sj_rows)))
  }
}

# Convert into dataframe
sd_kmat_elems <- as.data.frame(sd_kmat_elems)
sd_kmat_elems$st_dep <- sds_scale
sd_kmat_elems
# Convert to long
sd_kmat_elems_l <- sd_kmat_elems %>%
  pivot_longer(!st_dep, names_to = "elem", values_to = "value")

# Plot
ggplot(sd_kmat_elems_l, aes(x = st_dep, y = value, color = elem)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(
    values = c("red", "blue", "green", "black"),
    breaks = c("V1", "V2", "V3", "V4"),
    labels = c("A1->J1", "A2->J1", "A1->J2", "A2->J2")
  )
ggsave(here("StageDependentComp/figures/prelim/stoch_slowlh_elemvals_mult.png"), width = 8, height = 5)

# ===== Same but with sd_models =====

# Variable for storing the k_mat elements for each level of stage-dependence as a matrix
sd_kmat_elems <- NULL
# Extract relevant rows in the k_mat and insert into the matrix (Only those that weren't directly manipulated)
for(i in 1:length(sd_models)) {
  # Get row for sp1 juvenile survival
  sp1_sj_row <- sd_models[[i]]$k_mat$s[2,c(2,4)]
  # Get row for sp2 juvenile survival
  sp2_sj_row <- sd_models[[i]]$k_mat$s[5,c(2,4)]
  # Merge into single vector
  sj_rows <- c(sp1_sj_row, sp2_sj_row)
  # Insert into matrix
  if(is.null(sd_kmat_elems)) {
    sd_kmat_elems <- matrix(sj_rows, ncol = length(sj_rows))
  } else {
    sd_kmat_elems <- rbind(sd_kmat_elems, matrix(sj_rows, ncol = length(sj_rows)))
  }
}

# Convert into dataframe
sd_kmat_elems <- as.data.frame(sd_kmat_elems)
sd_kmat_elems$st_dep <- sds_scale
sd_kmat_elems
# Convert to long
sd_kmat_elems_l <- sd_kmat_elems %>%
  pivot_longer(!st_dep, names_to = "elem", values_to = "value")

# Plot
ggplot(sd_kmat_elems_l, aes(x = st_dep, y = value, color = elem)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(
    values = c("red", "blue", "green", "black"),
    breaks = c("V1", "V2", "V3", "V4"),
    labels = c("A1->J1", "A2->J1", "A1->J2", "A2->J2")
  )

# ====== Extras ======

# # Generate trajectory figure for stage-dependence of 0.5
# 
# # Calculate trajectory from initial population structures
# temp_traj <- sdcomp_project_dstb(sd_models[[6]], sp1_0, sp2_0, timestep = 400, dstb_func = dstb_func)
# # Convert trajectory data to long dataframe
# temp_traj_long <- temp_traj %>%
#   pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")
# # Plot trajectory
# ggplot(temp_traj_long, aes(x = time, y = density, color = species, linetype = stage)) +
#   geom_line() +
#   theme_bw()
# ggsave(here("StageDependentComp/figures/prelim/stoch_ann_0.5stdep_sim.png"), width = fig_def_width, height = fig_def_height)
# 
# # Generate trajectory figure for stage-dependence of 0.7
# 
# # Calculate trajectory from initial population structures
# temp_traj <- sdcomp_project_dstb(sd_models[[8]], sp1_0, sp2_0, timestep = 400, dstb_func = dstb_func)
# # Convert trajectory data to long dataframe
# temp_traj_long <- temp_traj %>%
#   pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")
# # Plot trajectory
# ggplot(temp_traj_long, aes(x = time, y = density, color = species, linetype = stage)) +
#   geom_line() +
#   theme_bw()
# ggsave(here("StageDependentComp/figures/prelim/stoch_ann_0.7stdep_sim.png"), width = fig_def_width, height = fig_def_height)
