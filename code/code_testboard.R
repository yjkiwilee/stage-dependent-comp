source("sdcompr.R")

pacman::p_load("tidyverse", "Rage")

n_stages <- 2

sp1_vrs <- list(
  s = c(1, 0.5, 0.8),
  g = c(0.6),
  f = c(4)
)

sp2_vrs <- list(
  s = c(1, 0.4, 0.8),
  g = c(0.6),
  f = c(4)
)

k_mat <- list(
  s = matrix(c(
    0, 0, 0, 0,
    0.007, 0.001, 0.004, 0.0005,
    0.001, 0.007, 0.0005, 0.004,
    0, 0, 0, 0,
    0.004, 0.0005, 0.007, 0.001,
    0.0005, 0.004, 0.001, 0.007
  ), ncol = 4, byrow = TRUE),
  g = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  f = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE)
)

k_mat_sd <- list(
  s = matrix(c(
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
  ), ncol = 5, byrow = TRUE),
  g = matrix(c(
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
  ), ncol = 5, byrow = TRUE),
  f = matrix(c(
    0, 0, 0, 0, 0.01,
    0, 0, 0, 0, 0.01
  ), ncol = 5, byrow = TRUE)
)

sdcomp_m1 <- def_sdcomp_model(n_stages, sp1_vrs, sp2_vrs, k_mat)
traj <- sdcomp_project(sdcomp_m1, c(100, 100), c(100, 100), timestep = 25)
traj_long <- traj %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")

ggplot(traj_long, aes(x = time, y = density, color = species, linetype = stage)) +
  geom_line()

sdcomp_m2 <- alter_k_mat(sdcomp_m1, "s", 2, 3, 0.01, iter = 70, scal_fact_l = -100, scal_fact_r = 200)
sdcomp_m2 <- alter_k_mat(sdcomp_m2, "s", 5, 1, 0.01, iter = 70, scal_fact_l = -100, scal_fact_r = 200)

sdcomp_m3 <- alter_k_mat(sdcomp_m1, "s", 5, 3, 0.001, iter = 70, scal_fact_l = -10, scal_fact_r = 200)
sdcomp_m3 <- alter_k_mat(sdcomp_m3, "s", 2, 1, 0.001, iter = 70, scal_fact_l = -10, scal_fact_r = 200)

calc_avg_int_str(sdcomp_m1, "s", 2, 1)
calc_avg_int_str(sdcomp_m2, "s", 2, 1)
sdcomp_m2
sdcomp_m3

traj_2 <- sdcomp_project(sdcomp_m2, c(100, 100), c(100, 100), timestep = 25)
traj_2_long <- traj_2 %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")

ggplot(traj_2_long, aes(x = time, y = density, color = species, linetype = stage)) +
  geom_line()


# Stoch
k_mat <- list(
  s = matrix(c(
    0, 0, 0, 0,
    0.02, 0.02, 0.007, 0.007,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0.007, 0.007, 0.02, 0.02,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  g = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE),
  f = matrix(c(
    0, 0, 0, 0,
    0, 0, 0, 0
  ), ncol = 4, byrow = TRUE)
)

sdcomp_stoch <- def_sdcomp_model(n_stages, sp1_vrs, sp2_vrs, k_mat, k_mat_sd)

traj_st <- sdcomp_project(sdcomp_stoch, c(100, 100), c(100, 100), timestep = 100)
traj_long_st <- traj_st %>%
  pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density")
sdcomp_stoch
ggplot(traj_long_st, aes(x = time, y = density, color = species, linetype = stage)) +
  geom_line()

sdcomp_project_end_avg(sdcomp_stoch, c(100, 100), c(100, 100), timestep = 4000, avgtf = 1000)

lsd <- traj_long_st %>%
  filter(time > 250 & time < 1000) %>%
  filter(species == "sp1" & stage == "s1")
lsd <- lsd$density

stage_dep_vec <- seq(0, 0.007, 0.000125)

sdcomp_sts <- list() # List for storing stage-dependent competition models

for(stage_dep in stage_dep_vec) {
  print(stage_dep)
  print("A")
  sdcomp_st_temp <- alter_k_mat(sdcomp_stoch, "s", 2, 1, 0.02 - stage_dep, iter = 20, scal_fact_l = 0, scal_fact_r = 100)
  print("B")
  sdcomp_st_temp <- alter_k_mat(sdcomp_st_temp, "s", 2, 3, 0.007 - stage_dep, iter = 20, scal_fact_l = 0, scal_fact_r = 100)
  print("C")
  sdcomp_st_temp <- alter_k_mat(sdcomp_st_temp, "s", 5, 3, 0.02 - stage_dep, iter = 20, scal_fact_l = 0, scal_fact_r = 100)
  print("D")
  sdcomp_st_temp <- alter_k_mat(sdcomp_st_temp, "s", 5, 1, 0.007 - stage_dep, iter = 20, scal_fact_l = 0, scal_fact_r = 100)
  sdcomp_sts[[as.character(stage_dep)]] <- sdcomp_st_temp
}

stoch_means <- c()
stoch_iqrs <- c()

for(stage_dep in stage_dep_vec) {
  print(stage_dep)
  
  sdcomp_st_temp <- sdcomp_sts[[as.character(stage_dep)]]
  
  traj_st_temp <- sdcomp_project(sdcomp_st_temp, c(30, 50), c(30, 50), timestep = 10000)
  
  traj_summ_temp <- traj_st_temp %>%
    pivot_longer(!time, names_sep = "_", names_to = c("species", "stage"), values_to = "density") %>%
    filter(time > 1000 & time < 10000) %>%
    filter(species == "sp1" & stage == "s1") %>%
    summarise(
      mean = mean(density),
      iqr = IQR(density)
    )
  stoch_means <- append(stoch_means, traj_summ_temp$mean)
  stoch_iqrs <- append(stoch_iqrs, traj_summ_temp$iqr)
}

stoch_iqrs
stoch_means

stage_dep_data <- tibble(
  stage_dep = stage_dep_vec,
  mean = stoch_means,
  iqr = stoch_iqrs
)

ggplot(stage_dep_data, aes(x = stage_dep, y = iqr)) +
  geom_point()

cor.test(stage_dep_data$stage_dep, stage_dep_data$iqr, method = "pearson")
