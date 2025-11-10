################################################################################
#' Script for calculating the crowding factor around individuals of the focal species.
#' ALPHA is specified outside the script.
#'
#' Written by Young Jun Lee
#' Nov 2024

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {
  # Unload stats packages that cause collision
  pacman::p_unload("arm", "lme4", "MASS")
  # Load census data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))
  
  # Load species vital rates
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_sp_vr.R"))
  
  # Load species stage thresholds
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_stage_thres.R"))
}

# Load relevant functions
source(here("StageDependentComp","code","census_analysis","func","stage_func.R"))
source(here("StageDependentComp","code","census_analysis","func","spatial_analysis.R"))

# ===== Constants =====

# Year & Plot number for test
TEST_PLOT_YEAR <- 2016
TEST_PLOT_ID <- "46"

# Alpha value to use (SPECIFIED EXTERNALLY)
# ALPHA <- 0.2

# Decay function to use; a step function
DECAY_FUN <- function(x, a) {
  ifelse(x <= a, 1, 0)
}
# DECAY_FUN <- function(x, a) exp(1)^(-a*(x^2))

# ===== Plot the crowding factor kernel =====

# # Populate the function output in a DF
# decay_fun_df <- tibble(
#   x = seq(0, 100, 1)
# ) %>%
#   mutate(y = DECAY_FUN(x, ALPHA))

# # Plot function
# decay_fun_plt <- ggplot(decay_fun_df, aes(x = x, y = y)) +
#   geom_line() +
#   theme_jun1() +
#   labs(
#     title = "Individual crowding factor against distance",
#     x = "Distance (cm)",
#     y = "Individual crowding factor",
#     caption = paste0("Alpha = ", ALPHA)
#   )
# 
# decay_fun_plt
# 
# # Save plot
# ggsave(here("StageDependentComp",
#             "figures",
#             "census_analysis",
#             "crowding",
#             "crowding_factor_test",
#             paste0("decay_kernel_step_", ALPHA, ".png")),
#        decay_fun_plt,
#        width = FIG_W, height = FIG_H)

# # ===== Subset test plot & visualise plot =====
# 
# # Test plot
# test_plot_data <- census_data %>%
#   filter(Year == TEST_PLOT_YEAR & Plot == TEST_PLOT_ID)
# 
# # Simplify data for plotting
# test_plot_simple_data <- test_plot_data %>%
#   mutate(
#     Taxon = as.factor(ifelse(Taxon %in% c(SP1, SP2), Taxon, "Other"))
#   ) %>%
#   mutate(
#     Taxon = fct_relevel(Taxon, c(SP1, SP2, "Other"))
#   )
# 
# # Visualise plot
# 
# # Plot points
# test_plot_plt <- ggplot() +
#   geom_point(data = test_plot_simple_data,
#              aes(x = X..cm., y = Y..cm., color = Taxon),
#              size = 1) +
#   coord_cartesian(xlim = c(0, 200), ylim = c(0, 200), expand = FALSE) +
#   scale_color_manual(
#     values = c("red", "blue", "grey")
#   ) +
#   theme_jun1() +
#   theme(aspect.ratio = 1) +
#   labs(
#     title = paste(c("Individuals in plot", TEST_PLOT_ID, "in", TEST_PLOT_YEAR), collapse = " "),
#     x = "X (cm)",
#     y = "Y (cm)"
#   )
# 
# test_plot_plt
# 
# # Save plot
# ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
#             "crowding_factor_test", paste0("test_plot_", ALPHA, ".png")),
#        test_plot_plt, width = FIG_W, height = FIG_H)
# 
# # ===== Calculate crowding factor =====
# 
# # Insert crowding factor columns into dataset
# test_crowd_data <- calc_crowd_coeff(test_plot_data, SP1, SP2, sp1_thres_len, sp2_thres_len, ALPHA, DECAY_FUN)
# 
# # Pivot longer
# test_crowd_data_long <- test_crowd_data %>%
#   pivot_longer(c(Sp1.S.crowding, Sp1.L.crowding, Sp2.S.crowding, Sp2.L.crowding, Other.crowding),
#                names_to = "Crowding.type",
#                values_to = "Crowding.factor") %>%
#   mutate(
#     Crowding.type = factor(Crowding.type)
#   ) %>%
#   mutate(
#     Crowding.type = fct_relevel(Crowding.type, c("Sp1.S.crowding", "Sp1.L.crowding", "Sp2.S.crowding", "Sp2.L.crowding", "Other.crowding"))
#   )
# 
# # ===== Visualise distribution of crowding factor =====
# 
# # Plot boxplot of crowding factors
# crowd_distrib_plt <- ggplot(test_crowd_data_long, aes(x = Crowding.factor, y = Crowding.type)) +
#   scale_y_discrete(breaks = c(
#       "Other.crowding", "Sp1.L.crowding", "Sp1.S.crowding",
#       "Sp2.L.crowding", "Sp2.S.crowding"),
#     labels = c(
#       "Other",
#       paste("Large,", SP2_SHORT),
#       paste("Small,", SP2_SHORT),
#       paste("Large,", SP1_SHORT),
#       paste("Small,", SP1_SHORT)
#   )) +
#   labs(
#     x = "Crowding factor",
#     y = "Category",
#     title = paste0("Distribution of crowding factor\nin plot ", TEST_PLOT_ID, ", ", TEST_PLOT_YEAR),
#     caption = paste0("N = ", nrow(test_crowd_data), ", Alpha = ", ALPHA)
#   ) +
#   geom_boxplot() +
#   theme_jun1()
# 
# # Save plot
# ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
#             "crowding_factor_test",
#             paste0("crowd_fct_distrib_", ALPHA, ".png")),
#        crowd_distrib_plt, width = FIG_W, height = FIG_H)
# 
# 
# ## ===== Plot individuals of the focal species only =====
# 
# crowd_distrib_plt_focal_only <- ggplot(test_crowd_data_long %>% filter(Taxon %in% c(SP1, SP2)), aes(x = Crowding.factor, y = Crowding.type)) +
#   scale_y_discrete(breaks = c(
#     "Other.crowding", "Sp1.L.crowding", "Sp1.S.crowding",
#     "Sp2.L.crowding", "Sp2.S.crowding"),
#     labels = c(
#       "Other",
#       paste("Large,", SP2_SHORT),
#       paste("Small,", SP2_SHORT),
#       paste("Large,", SP1_SHORT),
#       paste("Small,", SP1_SHORT)
#     )) +
#   labs(
#     x = "Crowding factor",
#     y = "Category",
#     title = paste0("Distribution of crowding factor\nin plot ", TEST_PLOT_ID, ", ", TEST_PLOT_YEAR, ", focal species only"),
#     caption = paste0("N = ", nrow(test_crowd_data %>% filter(Taxon %in% c(SP1, SP2))), ", Alpha = ", ALPHA)
#   ) +
#   geom_boxplot() +
#   theme_jun1()
# 
# crowd_distrib_plt_focal_only
# 
# # Save plot
# ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
#             "crowding_factor_test",
#             paste0("crowd_fct_distrib_focal_only_", ALPHA, ".png")),
#        crowd_distrib_plt_focal_only, width = FIG_W, height = FIG_H)
# 
# 
# ## ===== Plot focal category only, for focal species =====
# 
# crowd_distrib_plt_focal_cat_only <- ggplot(test_crowd_data_long %>% filter(Taxon %in% c(SP1, SP2) & Crowding.type != "Other.crowding"),
#                                            aes(x = Crowding.factor, y = Crowding.type)) +
#   scale_y_discrete(breaks = c(
#     "Sp1.L.crowding", "Sp1.S.crowding",
#     "Sp2.L.crowding", "Sp2.S.crowding"),
#     labels = c(
#       paste("Large,", SP2_SHORT),
#       paste("Small,", SP2_SHORT),
#       paste("Large,", SP1_SHORT),
#       paste("Small,", SP1_SHORT)
#     )) +
#   labs(
#     x = "Crowding factor",
#     y = "Category",
#     title = paste0("Distribution of focal crowding factor\nin plot ", TEST_PLOT_ID, ", ", TEST_PLOT_YEAR, ", focal species only"),
#     caption = paste0("N = ", nrow(test_crowd_data %>% filter(Taxon %in% c(SP1, SP2))), ", Alpha = ", ALPHA)
#   ) +
#   geom_boxplot() +
#   theme_jun1()
# 
# crowd_distrib_plt_focal_cat_only
# 
# # Save plot
# ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
#             "crowding_factor_test",
#             paste0("crowd_fct_distrib_focal_cat_only_", ALPHA, ".png")), 
#        crowd_distrib_plt_focal_cat_only, width = FIG_W, height = FIG_H)
# 
# ## ==== Plot crowding factor on histogram
# 
# # Labeller for facets
# crowd_type_labs <- c(
#   `Sp1.S.crowding` = paste("Small,", SP1_SHORT),
#   `Sp1.L.crowding` = paste("Large,", SP1_SHORT),
#   `Sp2.S.crowding` = paste("Small,", SP2_SHORT),
#   `Sp2.L.crowding` = paste("Large,", SP2_SHORT),
#   `Other.crowding` = "Other"
# )
# 
# crowd_distrib_hist_plt <- ggplot(test_crowd_data_long %>% filter(Taxon %in% c(SP1, SP2)),
#                                  aes(x = Crowding.factor)) +
#   geom_histogram(bins = 15) +
#   facet_wrap(~ Crowding.type,
#              scale = "free",
#              labeller = as_labeller(crowd_type_labs),
#              ncol = 2) +
#   theme_jun1() +
#   labs(
#     x = "Crowding factor",
#     y = "Count",
#     title = "Distribution of crowding factor in plot ", TEST_PLOT_ID, ", ", TEST_PLOT_YEAR,
#     caption = paste0("N = ", nrow(test_crowd_data %>% filter(Taxon %in% c(SP1, SP2))), ", Alpha = ", ALPHA, ", Focal individuals only")
#   )
# 
# crowd_distrib_hist_plt
# 
# # Save plot
# ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
#             "crowding_factor_test",
#             paste0("crowd_fct_distrib_hist_", ALPHA, ".png")),
#        crowd_distrib_hist_plt, width = FIG_W, height = FIG_H + 1.5)
# 
# # Commented code for scatter plot
# # ggplot(test_crowd_data, aes(x = Sp1.S.crowding, y = Sp1.L.crowding)) +
# #   geom_point()
# # 
# # ggplot(test_crowd_data, aes(x = Sp1.S.crowding)) +
# #   geom_histogram(binwidth = 0.1, boundary = 0) +
# #   theme_jun1()
# 
# # Test for normality
# # shapiro.test(test_crowd_data$Other.crowding)
# 
# # ===== Simulate data & populate null distribution for the crowding factor =====
# 
# # Number of iterations
# SIM_PLOT_N_ITER <- 100
# 
# # DF to store the simulated data
# sim_plot_data <- NULL
# 
# cat("Simulating null distribution for crowding factor...\n")
# 
# # For each iteration:
# for(i in 1:SIM_PLOT_N_ITER) {
#   # Temporary DF to store the batch of simulated data
#   sim_data_batch <- test_plot_data %>%
#     mutate(
#       Batch.ID = i # Insert batch ID
#     )
#   
#   # Randomly populate the X and Y coordinates
#   sim_data_batch <- sim_data_batch %>%
#     mutate(
#       X..cm. = runif(nrow(.), min = 0, max = 200),
#       Y..cm. = runif(nrow(.), min = 0, max = 200)
#     )
#   
#   # Calculate crowding
#   sim_crowd_data <- calc_crowd_coeff(sim_data_batch, SP1, SP2, sp1_thres_len, sp2_thres_len, ALPHA, DECAY_FUN)
#   
#   # Bind to sim_plot_data
#   if(is.null(sim_plot_data)) {
#     sim_plot_data <- sim_crowd_data
#   } else {
#     sim_plot_data <- bind_rows(sim_plot_data, sim_crowd_data)
#   }
#   
#   # cat(paste0(i, " / ", SIM_PLOT_N_ITER, " iterations done\n"))
# }
# 
# # cat("Plotting simulated dataset\n")
# 
# # Pivot longer
# sim_plot_data_long <- sim_plot_data %>%
#   pivot_longer(c(Sp1.S.crowding, Sp1.L.crowding, Sp2.S.crowding, Sp2.L.crowding, Other.crowding),
#                names_to = "Crowding.type",
#                values_to = "Crowding.factor")
# 
# # Plot simulated data
# sim_crowd_distrib_hist_plt <- ggplot(sim_plot_data_long %>% filter(Taxon %in% c(SP1, SP2)), aes(x = Crowding.factor)) +
#   geom_histogram(bins = 15) +
#   facet_wrap(~ factor(Crowding.type, levels = c("Sp1.S.crowding", "Sp1.L.crowding",
#                                                 "Sp2.S.crowding", "Sp2.L.crowding",
#                                                 "Other.crowding")),
#              scale = "free",
#              labeller = as_labeller(crowd_type_labs),
#              ncol = 2) +
#   theme_jun1() +
#   labs(
#     x = "Crowding factor",
#     y = "Count",
#     title = paste0("Simulated null distribution of crowding factor in plot ", TEST_PLOT_ID, ", ", TEST_PLOT_YEAR),
#     caption = paste0("N = ", nrow(sim_plot_data %>% filter(Taxon %in% c(SP1, SP2))),
#                                   ", Alpha = ", ALPHA,
#                                   ", N. iterations = ", SIM_PLOT_N_ITER,
#                                   ", Focal individuals only")
#   )
# 
# sim_crowd_distrib_hist_plt
# 
# # Save plot
# ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
#             "crowding_factor_test",
#             paste0("crowd_fct_null_distrib_", ALPHA, ".png")),
#        sim_crowd_distrib_hist_plt, width = FIG_W, height = FIG_H + 1.5)

# ===== Calculate crowding factor for all individuals in all plots =====

# Obtain plots & years with at least one of the focal individual
focal_plts <- census_data %>%
  filter(Taxon %in% c(SP1, SP2) & Died.this.census.final == 0) %>%
  group_by(Year, Plot) %>%
  summarise(
    n = n()
  ) %>%
  filter(n > 0) %>%
  select(Year, Plot)

# DF to store the processed rows
crowd_fct_df <- NULL

cat("Starting to calculate crowding factors for plots...\n")

# Iterate through these plots and calculate the crowding factor for focal species
for(i in 1:nrow(focal_plts)) {
  # Extract focal plot info
  focal_yr <- focal_plts[[i, "Year"]]
  focal_plt_id <- focal_plts[[i, "Plot"]]
  
  # Obtain corresponding subset from the census data
  plot_sub_data <- census_data %>%
    filter(Year == focal_yr & Plot == focal_plt_id)
  
  # Insert crowding factor columns into dataset & filter for live indivs of focal species
  plot_sub_data <- calc_crowd_coeff(plot_sub_data, SP1, SP2, sp1_thres_len, sp2_thres_len, ALPHA, DECAY_FUN) %>%
    filter(Taxon %in% c(SP1, SP2))
  
  # Rbind with crowd_fct_df
  if(is.null(crowd_fct_df)) {
    crowd_fct_df <- plot_sub_data
  } else {
    crowd_fct_df <- bind_rows(crowd_fct_df, plot_sub_data)
  }
  
  # cat(paste0(i, " / ", nrow(focal_plts), " plots processed\n"))
}

crowd_fct_df <- crowd_fct_df %>%
  select( # Only select crowding factors and identifying variablees
    c(Year, Plot, Tag, Taxon,
      Sp1.crowding, Sp1.S.crowding, Sp1.L.crowding,
      Sp2.crowding, Sp2.S.crowding, Sp2.L.crowding,
      Other.crowding, Total.crowding)
  )

# Save resulting df
write_csv(crowd_fct_df,
          here("StageDependentComp", "result_data", "wp2",
               "crowding_factor_test",
               paste0("focal_spp_crowding_", ALPHA, ".csv")))

# ===== Merge df with crowding with df with vital rate information =====

# Merge species data
sp_data <- bind_rows(sp1_data, sp2_data)

# Add ID column to crowd_fct_df and sp_data
crowd_fct_df_temp <- crowd_fct_df %>%
  mutate(
    ID = paste(Year, Plot, Tag)
  ) %>%
  select( # Only select ID and crowding factors
    !c(Year, Plot, Tag, Taxon)
  )
sp_data_temp <- sp_data %>%
  mutate(
    ID = paste(Year, Plot, Tag)
  )

crowd_fct_df_temp
sp_data_temp

# Join two dfs
merged_focal_data <- full_join(sp_data_temp, crowd_fct_df_temp, by = join_by(ID == ID)) %>%
  select(!ID) # Remove ID column

# Save resulting df
write_csv(merged_focal_data,
          here("StageDependentComp", "result_data", "wp2",
               "crowding_factor_test",
               paste0("focal_spp_full_data_", ALPHA, ".csv")))
