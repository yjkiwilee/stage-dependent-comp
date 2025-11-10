# ====== Testing stochastic models across a range of stage-dependence in ~annual species ======

# ====== 1. Define model parameters & setup ======

pacman::p_unload("MASS")

# Exported figure dimensions
fig_def_width <- 6
fig_def_height <- 4

# All paths are relative to the resilienceNERC master folder

# Number of life cycle stages
n_stages <- 2

# Density-independent vital rates of species 1 (Fast lh species) & associated st dev
sp1_vrs <- list(
  s = c(0.9, 0.95), # Small & large survival
  g = c(0.3), # Maturation rate
  r = c(0.1), # Retrogression rate
  f = c(0.3314815) # Fertility
)

sp1_eq_vrs <- list(
  s = c(0.9, 0.95), # Small & large survival
  g = c(0.3), # Maturation rate
  r = c(0.1), # Retrogression rate
  f = c(0.103704) # Fertility
)

stable.stage(vrs_to_mpm(2, sp1_eq_vrs)$matA)

sp1_vr_sd <- list(
  s = c(0.2, 0.2), # Small & large survival
  g = c(0.2), # Maturation rate
  r = c(0.2), # Retrogression rate
  f = c(0.1) # Fecundity
)

# Density-independent vital rates of species 1 (Fast lh species) & associated st dev
sp2_vrs <- list(
  s = c(0.4, 0.6), # Small & large survival
  g = c(0.9), # Maturation rate
  r = c(0.05), # Retrogression rate
  f = c(1.530556) # Fertility
)

sp2_vr_sd <- list(
  s = c(0.2, 0.2), # Small & large survival
  g = c(0.2), # Maturation rate
  r = c(0.2), # Retrogression rate
  f = c(0.1) # Fecundity
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
intrasp_comp <- -0.0003
intersp_comp <- -0.0001

# K matrices describing stage-dependent competition for each of the vital rates
# Coffeicients are stage-independent
# k_mat <- list(
#   s = matrix(c(
#     -0.001, -0.001, -0.0003, -0.0003, # Only juvenile survival is stage-dependent, with greater intraspecific competition
#     0, 0, 0, 0,
#     -0.0003, -0.0003, -0.001, -0.001,
#     0, 0, 0, 0
#   ), ncol = 4, byrow = TRUE),
#   g = matrix(c( # No density-dependence on maturation
#     0, 0, 0, 0,
#     0, 0, 0, 0
#   ), ncol = 4, byrow = TRUE),
#   r = matrix(c( # No density-dependence on retrogression
#     0, 0, 0, 0,
#     0, 0, 0, 0
#   ), ncol = 4, byrow = TRUE),
#   f = matrix(c( # Or fecundity
#     0, 0, 0, 0,
#     0, 0, 0, 0
#   ), ncol = 4, byrow = TRUE)
# )

# k_mat <- list(
#   s = matrix(c(
#     0, 0, 0, 0,
#     0, 0, 0, 0,
#     0, 0, 0, 0,
#     0, 0, 0, 0
#   ), ncol = 4, byrow = TRUE),
#   g = matrix(c( # No density-dependence on retrogression
#     -0.003, -0.003, -0.0001, -0.0001,
#     -0.00001, -0.0001, -0.003, -0.003
#   ), ncol = 4, byrow = TRUE),
#   r = matrix(c(
#     0, 0, 0, 0,
#     0, 0, 0, 0
#   ), ncol = 4, byrow = TRUE),
#   f = matrix(c( # Or fecundity
#     0, 0, 0, 0,
#     0, 0, 0, 0
#   ), ncol = 4, byrow = TRUE)
# )

k_mat <- list(
  s = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  g = matrix(c( # No density-dependence on retrogression
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  r = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  f = matrix(c( # Or fecundity
    -0.001, -0.001, -0.00023, -0.00023,
    -0.00023, -0.00023, -0.001, -0.001
  ), ncol = 4, byrow = TRUE)
)

alpha <- 1

k_mat_stdep <- list(
  s = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  g = matrix(c( # No density-dependence on retrogression
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  r = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  f = matrix(c( # Or fecundity
    -0.001 * (1 - alpha), -0.001 * (1 + alpha * 0.54), -0.00023 * (1 - alpha), -0.00023 * (1 + alpha * 1.21),
    -0.00023 * (1 - alpha), -0.00023 * (1 + alpha * 0.54), -0.001 * (1 - alpha), -0.001 * (1 + alpha * 1.21)
  ), ncol = 4, byrow = TRUE)
)

# k_mat <- list(
#   s = matrix(c(
#     intrasp_comp, intrasp_comp, intersp_comp, intersp_comp,
#     intrasp_comp, intrasp_comp, intersp_comp, intersp_comp,
#     intersp_comp, intersp_comp, intrasp_comp, intrasp_comp,
#     intersp_comp, intersp_comp, intrasp_comp, intrasp_comp
#   ), ncol = 4, byrow = TRUE),
#   g = matrix(c(
#     intrasp_comp, intrasp_comp, intersp_comp, intersp_comp,
#     intersp_comp, intersp_comp, intrasp_comp, intrasp_comp
#   ), ncol = 4, byrow = TRUE),
#   r = matrix(c(
#     -intrasp_comp, -intrasp_comp, -intersp_comp, -intersp_comp,
#     -intersp_comp, -intersp_comp, -intrasp_comp, -intrasp_comp
#   ), ncol = 4, byrow = TRUE),
#   f = matrix(c(
#     intrasp_comp/2, intrasp_comp/2, intersp_comp/2, intersp_comp/2,
#     intersp_comp/2, intersp_comp/2, intrasp_comp/2, intrasp_comp/2
#   ), ncol = 4, byrow = TRUE)
# )

# Define model with these parameters
st_indep_m <- def_sdcomp_model(n_stages, sp1_vrs, sp2_vrs, k_mat,
                               sp1_vr_sd, sp2_vr_sd)

# Get density-independent matrices
sp1_mpm <- get_sp_mpm(st_indep_m, 1)
sp2_mpm <- get_sp_mpm(st_indep_m, 2)
sp1_mpm

# Get population metrics
generation.time(sp1_mpm$matA)
generation.time(sp2_mpm$matA)
life_expect_mean(sp1_mpm$matU)
life_expect_mean(sp2_mpm$matU)
lambda(sp1_mpm$matA)
lambda(sp2_mpm$matA)
net.reproductive.rate(sp1_mpm$matA)
net.reproductive.rate(sp2_mpm$matA)

# temp_model <- alter_k_mat_mult(st_indep_m, vr_type = "f",
#                                mat_rows = c(1,1,2,2), mat_cols = c(1,3,1,3),
#                                alt_vals = c(-0.0005, -0.000025, -0.000025, -0.0005))
temp_model <- def_sdcomp_model(n_stages, sp1_vrs, sp2_vrs, k_mat_stdep,
                               sp1_vr_sd, sp2_vr_sd)

calc_avg_int_str_mult(temp_model, "f", mat_rows = c(1,1,2,2), species = c(1,2,1,2))
calc_avg_int_str_mult(st_indep_m, "f", mat_rows = c(1,1,2,2), species = c(1,2,1,2))

x <- c()
y <- c()
env_cond <- c()

for(i in 1:50) {
  cat(sprintf("%d\n", i))
  
  set.seed(i)
  env_cond <- c(env_cond, rnorm(1))
  sp1_rnorms <- lapply(1:1000, function(t) {
    list(
      s = rnorm(2), # Small & large survival
      g = rnorm(1), # Maturation rate
      r = rnorm(1), # Retrogression rate
      f = rnorm(1) # Fecundity
    )
  })
  sp2_rnorms <- lapply(1:1000, function(t) {
    list(
      s = rnorm(2), # Small & large survival
      g = rnorm(1), # Maturation rate
      r = rnorm(1), # Retrogression rate
      f = rnorm(1) # Fecundity
    )
  })
  
  test_project <- sdcomp_project_vr(temp_model, sp1_0, sp2_0, timestep = 1000,
                                    sp1_rnorms, sp2_rnorms)
  stindep_project <- sdcomp_project_vr(st_indep_m, sp1_0, sp2_0, timestep = 1000,
                                       sp1_rnorms, sp2_rnorms)
  
  x <- c(x, var(test_project$sp1_s1[500:1000] + test_project$sp1_s2[500:1000]))
  y <- c(y, var(stindep_project$sp1_s1[500:1000] + stindep_project$sp1_s2[500:1000]))
}

test_project$vr_sp1_f1 / stindep_project$vr_sp1_f1

plot(x, y)

t.test(x, y, paired = TRUE)
mean(sqrt(x))
mean(sqrt(y))

# var.test(stindep_project$sp2_s1[500:5000] + stindep_project$sp2_s2[500:5000],
#          test_project$sp2_s1[500:5000] + test_project$sp2_s2[500:5000])

set.seed(1)
sp1_rnorms <- lapply(1:1000, function(t) {
  list(
    s = rnorm(2), # Small & large survival
    g = rnorm(1), # Maturation rate
    r = rnorm(1), # Retrogression rate
    f = rnorm(1) # Fecundity
  )
})
sp2_rnorms <- lapply(1:1000, function(t) {
  list(
    s = rnorm(2), # Small & large survival
    g = rnorm(1), # Maturation rate
    r = rnorm(1), # Retrogression rate
    f = rnorm(1) # Fecundity
  )
})

test_project <- sdcomp_project_vr(temp_model, sp1_0, sp2_0, timestep = 1000,
                                  sp1_rnorms, sp2_rnorms)
stindep_project <- sdcomp_project_vr(st_indep_m, sp1_0, sp2_0, timestep = 1000,
                                     sp1_rnorms, sp2_rnorms)

test_proj_comm <- long_comm_from_proj(test_project) %>%
  group_by(species, stage) %>%
  filter(!is.na(size))
# test_proj_vr <- extract_vr_from_proj(test_project)
sdin_proj_comm <- long_comm_from_proj(stindep_project) %>%
  group_by(species, stage) %>%
  filter(!is.na(size))

ggplot(sdin_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
  geom_line()
ggplot(test_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
  geom_line()

# ggplot(sdin_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
#   geom_line()

ggplot(stindep_project, aes(x = time, y = sp2_s1 + sp2_s2)) +
  geom_line()

ggplot(test_project, aes(x = time, y = sp2_s1 + sp2_s2)) +
  geom_line()

ggplot(stindep_project, aes(x = time, y = sp2_s1 + sp2_s2)) +
  geom_line()

ggplot(test_project, aes(x = time, y = sp2_s1 + sp2_s2)) +
  geom_line()




sdin_end <- sdcomp_project_end_avg(st_indep_m, sp1_0, sp2_0, 10000, 5000)
sdin_end$sp1_pop[1] / sdin_end$sp1_pop[2]
sdin_end$sp2_pop[1] / sdin_end$sp2_pop[2]
test_end <- sdcomp_project_end_avg(temp_model, sp1_0, sp2_0, 10000, 5000)
test_end$sp1_pop[1] / (sdin_end$sp1_pop[1] + sdin_end$sp1_pop[2])

mean(stindep_project[500:900,"sp1_s1"] + stindep_project[500:900,"sp1_s2"])
mean(stindep_project[500:900,"sp2_s1"] + stindep_project[500:900,"sp2_s2"])






# 
# arima(stindep_project$sp1_s1[400:900], c(9, 1, 0))
# arima(stindep_project$sp1_s1[400:900], c(9, 1, 0))
# AIC(
#   arima(stindep_project$sp1_s1[400:900], c(0, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(1, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(2, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(3, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(4, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(5, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(6, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(7, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(8, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(9, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(10, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(11, 1, 0)),
#   arima(stindep_project$sp1_s1[400:900], c(12, 1, 0))
# )
# 
# cpgram(resid(arima(stindep_project$sp1_s1, c(1, 0, 0))), main = "Cumulative periodogram of residuals")
# cpgram(resid(arima(stindep_project$sp2_s1, c(1, 0, 0))), main = "Cumulative periodogram of residuals")

# ggplot(test_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
#   geom_line()




sdin_proj_demo <- demo_from_proj(stindep_project, 2)
sdin_proj_demo

ggplot(sdin_proj_demo, aes(x = time, y = lambda, color = species)) +
  geom_line() +
  geom_vline(xintercept = 100) +
  geom_smooth()

test_time <- 96
t.test(sdin_proj_demo$lambda[test_time:(test_time+100)] - 1)


test_tseries <- (sdin_proj_comm %>%
  filter(species == "1", stage == "1") %>%
  arrange(time))$size

t.test(test_tseries[500:200], test_tseries[200:300],
       alternative = "two.sided")


ggplot(test_proj_vr %>% filter(vr_type == "s"), aes(x = vr_val)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  facet_wrap(~ vr_name)

stindep_proj_comp <- stindep_project
stindep_proj_comp$time <- stindep_proj_comp$time + 1
stindep_proj_comp <- full_join(stindep_project, stindep_proj_comp,
                               by = join_by(time), suffix = c("", "_prev"))

ggplot(stindep_proj_comp %>% filter(time > 1000), aes(x = vr_sp1_s1, y = sp1_s1 / sp1_s1_prev)) +
  geom_point()

test_lm <- lm(sp1_s1 ~ sp1_s1_prev * vr_sp1_s1, data = stindep_proj_comp)

summary(test_lm)

# test_proj_comm <- long_comm_from_proj(test_project)
# test_proj_vr <- long_vr_from_proj(test_project)
# 
# ggplot(test_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
#   geom_line()
# 
# ggplot(test_proj_vr %>% filter(vr_type == "s"), aes(x = vr_val)) +
#   geom_histogram(binwidth = 0.05, boundary = 0) +
#   facet_wrap(~ vr_name)
# 


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
                                 mat_rows = c(1, 1, 3, 3),
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
stage_labs <- c("Small", "Large")
names(stage_labs) <- c("s1", "s2")

ggplot(sd_models_slh_m_long, aes(x = st_dep, y = avg)) +
  geom_point() +
  theme_jun2() +
  labs(
    x = "Degree of stage-dependence in competition",
    y = "Long-term mean of density"
  ) +
  facet_grid(stage ~ species, scales = "free_y", labeller = labeller(species = species_labs, stage = stage_labs))

ggsave(here("StageDependentComp","figures","prelim","stoch_slowlh_varstdep_mean_zoom_mult.png"), width = 8, height = 5)

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
stage_labs <- c("Small", "Large")
names(stage_labs) <- c("s1", "s2")

# y limit setter
# ylim_setter <- tibble(
#   st_dep = c(0.2, 0.2, 0.2, 0.2),
#   sd = c(0, 10, 0, 5),
#   species = c("sp1", "sp1", "sp1", "sp1"),
#   stage = c("s1", "s1", "s2", "s2")
# )
ylim_setter <- tibble(
  st_dep = c(0.2, 0.2, 0.2, 0.2),
  sd = c(25, 53, 8, 25),
  species = c("sp1", "sp1", "sp1", "sp1"),
  stage = c("s1", "s1", "s2", "s2")
)

ggplot(sd_models_slh_sd_long, aes(x = st_dep, y = sd, color = species)) +
  geom_point(size = 0.7) +
  geom_point(data = ylim_setter, alpha = 0) +
  theme_jun1() +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = "Degree of stage-dependence",
    y = "Standard deviation of density"
  ) +
  guides(color = "none") +
  facet_grid(stage ~ species, scales = "free_y", labeller = labeller(species = species_labs, stage = stage_labs)) +
  theme(
    # panel.background = element_rect(fill='transparent'), #transparent panel bg
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent', color=NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color=NA), #transparent legend panel
  )

# Save plot
ggsave(here("StageDependentComp","figures","prelim","stoch_slowlh_varstdep_sd_zoom_mult.png"), width = 4, height = 3.8, scale = 0.8)

# ====== B5. Test Spearman correlation =======

for(sp in c("sp1", "sp2")) {
  for(st in c("s1", "s2")) {
    sd_models_slh_sd_sub <- sd_models_slh_sd_long %>%
      filter(species == sp & stage == st)
    
    print(paste(sp, st))
    sd_model_cor_test <- cor.test(sd_models_slh_sd_sub$st_dep, sd_models_slh_sd_sub$sd,
                                  method = "spearman")
    print(sd_model_cor_test)
    print(nrow(sd_models_slh_sd_sub))
  }
}

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
ggsave(here("StageDependentComp/figures/prelim/stoch_ann_elemvals_mult.png"), width = 8, height = 5)

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


# ======= Scrap test ========


set.seed(1)
stoch_rnorms <- matrix(
  rnorm(1000 * 2 * (2+1+1+1)), ncol = 2 * (2+1+1+1)
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
stoch_rnorms

stoch_rnorms_sp1 <- lapply(1:1000, function(t) {
  list(
    s = c(stoch_rnorms[t, c(1,2)]), # Small & large survival
    g = stoch_rnorms[t, 3], # Maturation rate
    r = stoch_rnorms[t, 4], # Retrogression rate
    f = stoch_rnorms[t, 5] # Fecundity
  )
})

stoch_rnorms_sp2 <- lapply(1:1000, function(t) {
  list(
    s = c(stoch_rnorms[t, c(6,7)]), # Small & large survival
    g = stoch_rnorms[t, 8], # Maturation rate
    r = stoch_rnorms[t, 9], # Retrogression rate
    f = stoch_rnorms[t, 10] # Fecundity
  )
})

sd_test_model <- alter_stdep(
  base_comm_models$sj$F2S2,
  0, "sj", pop_eq_structs$F2$sj, pop_eq_structs$S2$sj
)

sd_test_model

# test_project <- sdcomp_project_vr(base_comm_models$sj$F2S2,
#                                   n_eq * pop_eq_structs$F2$sj * 0.05,
#                                   n_eq * pop_eq_structs$S2$sj * 0.05,
#                                   timestep = 1000,
#                                   stoch_rnorms_sp1, stoch_rnorms_sp2)

test_project <- sdcomp_project_vr(sd_test_model,
                                  n_eq * pop_eq_structs$F2$sj * 0.05,
                                  n_eq * pop_eq_structs$F2$sj * 0.05,
                                  timestep = 1000,
                                  stoch_rnorms_sp1, stoch_rnorms_sp2)

test_proj_comm <- long_comm_from_proj(test_project) %>%
  group_by(species, stage) %>%
  filter(!is.na(size))

ggplot(test_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
  geom_line()

test_project$sp1 <- test_project$sp1_s1 + test_project$sp1_s2
test_project$sp2 <- test_project$sp2_s1 + test_project$sp2_s2
test_project$sp1_struct <- test_project$sp1_s1 / test_project$sp1
test_project$sp2_struct <- test_project$sp2_s1 / test_project$sp2

future_sp1 <- c(test_project$sp1[2:(length(test_project$sp1))], NA)
test_project$future_sp1 <- future_sp1
future_sp2 <- c(test_project$sp2[2:(length(test_project$sp2))], NA)
test_project$future_sp2 <- future_sp2

test_project <- bind_cols(test_project, as.data.frame(stoch_rnorms))

mean(test_project$sp1[250:900])
mean(test_project$sp2[250:900])

nrow(test_project)

test_project

test_lm <- lm(
  (future_sp2 / sp2) ~ sp1 + sp2 + sp1_struct + sp2_struct +
    stoch_sj_2 + stoch_sa_2 + stoch_g_2 + stoch_r_2 + stoch_f_2,
  test_project[250:900,]
)

summary(test_lm)
Anova(test_lm)


# Anova Table (Type II tests)
# 
# Response: (future_sp2/sp2)
# Sum Sq  Df    F value    Pr(>F)    
# sp1        0.00187   1    55.3562 3.245e-13 ***
#   sp2        0.05447   1  1615.6249 < 2.2e-16 ***
#   sp1_struct 0.00000   1     0.1404    0.7080    
# sp2_struct 0.16914   1  5016.4432 < 2.2e-16 ***
#   stoch_sj_2 0.25576   1  7585.4034 < 2.2e-16 ***
#   stoch_sa_2 0.01984   1   588.3338 < 2.2e-16 ***
#   stoch_g_2  0.00005   1     1.4111    0.2353    
# stoch_r_2  0.00008   1     2.5043    0.1140    
# stoch_f_2  0.92382   1 27399.0477 < 2.2e-16 ***
#   Residuals  0.02161 641                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# test_project <- sdcomp_project_vr(sd_test_model,
#                                   n_eq * pop_eq_structs$F2$sj, n_eq * pop_eq_structs$S2$sj,
#                                   timestep = 1000,
#                                   stoch_rnorms, stoch_rnorms)



test_project <- sdcomp_project_vr(base_comm_models$sj$F2F2,
                                  n_eq * pop_eq_structs$F2$sj * 0.05,
                                  n_eq * pop_eq_structs$F2$sj * 0.05,
                                  timestep = 1000,
                                  stoch_rnorms_sp1, stoch_rnorms_sp2)

test_project

test_proj_comm <- long_comm_from_proj(test_project) %>%
  group_by(species, stage) %>%
  filter(!is.na(size))

ggplot(test_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
  geom_line()

ggplot(test_project, aes(x = time, y = sp1_s1+sp1_s2))

sd_test_model <- alter_stdep(
  base_comm_models$sj$F2F2,
  0, "sj", pop_eq_structs$F2$sj, pop_eq_structs$F2$sj
)

sd_test_model

sd_project <- sdcomp_project_vr(sd_test_model,
                                n_eq * pop_eq_structs$F2$sj * 0.05,
                                n_eq * pop_eq_structs$F2$sj * 0.05,
                                timestep = 1000,
                                stoch_rnorms_sp1, stoch_rnorms_sp2)

sd_project

sd_proj_comm <- long_comm_from_proj(sd_project) %>%
  group_by(species, stage) %>%
  filter(!is.na(size))

ggplot(sd_proj_comm, aes(x = time, y = size, color = species, linetype = stage)) +
  geom_line()
