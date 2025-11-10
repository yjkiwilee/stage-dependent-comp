# ===== File for preliminary model simulations for deterministic models, ~annual species =====

# ==== 1. Test run of the deterministic model =====

# Load required packages (Rage is not used here but NB the functions were written to generate MPMs that are compatible with Rage)
pacman::p_load("tidyverse", "Rage", "here")

# Exported figure dimensions
fig_def_width <- 4.8
fig_def_height <- 4

# All paths are relative to the resilienceNERC master folder

# Load R file with model functions
source(here("StageDependentComp/code/sdcompr/sdcompr.R"))
# Load R file with custom ggplot themes
source(here("StageDependentComp/code/themes/custom_themes.R"))

# Number of life cycle stages
n_stages <- 2

# Density-independent vital rates of species 1 (Fast life history)
sp1_vrs <- list(
  s = c(1, 0.8, 0.1), # Seed, juvenile, adult survival
  g = c(0.9), # Maturation rate
  f = c(2)
)

# Density-independent vital rates of species 2 (Fast life history)
sp2_vrs <- list(
  s = c(1, 0.8, 0.1),
  g = c(0.9),
  f = c(2)
)

# Initial population structures for simulations
sp1_0 <- c(500, 0)
sp2_0 <- c(100, 0)

# Simulation lengths for the figures
sim_timestep <- 50

# ====== 1.1. Deterministic stage-independent competition model ======

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

# Define model with these parameters
st_indep_m <- def_sdcomp_model(n_stages, sp1_vrs, sp2_vrs, k_mat)
# Calculate trajectory from initial population structures
st_indep_traj <- sdcomp_project(st_indep_m, sp1_0, sp2_0, timestep = sim_timestep)
# Convert trajectory data to long dataframe
st_indep_traj_long <- st_indep_traj %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")
# Plot trajectory
ggplot(st_indep_traj_long, aes(x = time, y = density, color = species, linetype = stage)) +
  geom_line() +
  labs(
    x = "Time",
    y = "Density",
    linetype = "Stage",
    color = "Species"
  ) +
  scale_linetype_discrete(
    labels = c("Adults", "Juveniles")
  ) +
  scale_color_discrete(
    labels = c("Sp. 1", "Sp. 2")
  ) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 0)) +
  theme_jun1()

# Save figure
ggsave(here("StageDependentComp/figures/prelim/det_stindep_fastlh.png"), width = fig_def_width, height = fig_def_height)

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
  labs(
    x = "Time",
    y = "Adult-to-juvenile ratio",
    color = "Species"
  ) +
  theme_jun1()
# Save figure
ggsave(here("StageDependentComp/figures/prelim/det_stindep_fastlh_ajr.png"), width = fig_def_width, height = fig_def_height)

# ====== 1.2. Deterministic stage-dependent competition model ======

# Alter elements of k_mat for stage-dependent competition while keeping the species-average per-capita competitive strength
# at equilibrium population structure constant
# Here, the juvenile survival is being set to be stage-dependent
st_dep_m <- alter_k_mat_mult(st_indep_m,
                        vr_type = "s", # Alter survival component of k_mat
                        mat_rows = c(2,2,5,5),
                        mat_cols = c(1,3,1,3),
                        alt_vals = c(0.0006, 0.0002, 0.0002, 0.0006), # Set it so that the adults' effect on juvenile survival is now greater than that of juveniles
                        iter = 70, # Number of iterations to numerically determine the appropriate k_mat configuration
                        scal_fact_l = -100,
                        scal_fact_r = 100) # Minimum and maximum scaling factors to try in the first instance; the range should contain the desired scaling factor.

# Calculate trajectory from initial population structures
st_dep_traj <- sdcomp_project(st_dep_m, sp1_0, sp2_0, timestep = sim_timestep)
# Convert trajectory data to long dataframe
st_dep_traj_long <- st_dep_traj %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")
# Plot trajectory
ggplot(st_dep_traj_long, aes(x = time, y = density, color = species, linetype = stage)) +
  geom_line() +
  labs(
    x = "Time",
    y = "Density",
    linetype = "Stage",
    color = "Species"
  ) +
  scale_linetype_discrete(
    labels = c("Adults", "Juveniles")
  ) +
  scale_color_discrete(
    labels = c("Sp. 1", "Sp. 2")
  ) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 0)) +
  theme_jun1()

ggsave(here("StageDependentComp/figures/prelim/det_stdep_fastlh.png"), width = fig_def_width, height = fig_def_height)

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
ggsave(here("StageDependentComp/figures/prelim/det_stdep_fastlh_ajr.png"), width = fig_def_width, height = fig_def_height)

# ===== 1.3. Generate plots to compare stage-dependent and stage-independent simulations =====

# Merge trajectories & convert to long
merged_traj <- tibble(
  sp1_stindep = st_indep_traj$sp1_s1 + st_indep_traj$sp1_s2,
  sp2_stindep = st_indep_traj$sp2_s1 + st_indep_traj$sp2_s2,
  sp1_stdep = st_dep_traj$sp1_s1 + st_dep_traj$sp1_s2,
  sp2_stdep = st_dep_traj$sp2_s1 + st_dep_traj$sp2_s2,
  time = st_indep_traj$time
)
merged_traj_long <- merged_traj %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage_dep"), values_to = "density")

# Plot trajectory
ggplot(merged_traj_long, aes(x = time, y = density, color = species, linetype = stage_dep)) +
  geom_line() +
  labs(
    x = "Time",
    y = "Density",
    linetype = "Competition",
    color = "Species"
  ) +
  scale_linetype_manual(
    labels = c("Stage-dependent", "Stage-independent"),
    breaks = c("stdep", "stindep"),
    values = c(5, 1)
  ) +
  scale_color_brewer(palette = "Set2", labels = c("Sp. 1", "Sp. 2")) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 0)) +
  theme_jun2()

ggsave(here("StageDependentComp/figures/prelim/det_fastlh_comparison.png"), width = 5.5, height = 3.3)
