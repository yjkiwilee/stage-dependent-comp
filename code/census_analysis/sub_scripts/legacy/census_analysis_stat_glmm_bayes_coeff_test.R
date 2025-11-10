################################################################################
#' Script for fitting a GLMM to the demographic & environmental data
#' Using a Bayesian approach
#' This is different from census_analysis_stat_glmm in the specific models;
#' refer to research log on 27022025
#' Written by Young Jun Lee
#' Feb 2024
#' https://cran.r-project.org/web//packages//glmm/vignettes/vignettes.pdf


# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {
  # Load stats packages
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_stat_packages.R"))
  # Load focal species vital rate data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_full_sp_data.R"))
  # Load statistical data; conditional upon census_analysis_stat_glmm.R!
  
  # Load annual environmental data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_ann_env_data.R"))
  
  # Load recruitment factor data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_sp_rec_fact.R"))
  
  # Load functions for stats
  source(here("StageDependentComp", "code", "census_analysis", "func",
              "stat_fitting.R"))
}

# ===== Constants =====

# Width of the boundary area to exclude individuals from (to reduce edge effect)
# In cma
BOUNDARY_WIDTH = 25

# Plot dimension, in cm
PLOT_DIM = 200

# ===== Import & clean data =====

# Filter focal species vital rate data to only include live individuals with
# Non-NA vital rates
spp_stat_data <- focal_spp_data %>%
  filter(
    Died.this.census.final == 0 &
    !is.na(Survival.future) &
    !is.na(Growth.future) &
    !is.na(X..Capitulescences) &
    !is.na(Stage) &
    X..cm. >= BOUNDARY_WIDTH & X..cm. <= PLOT_DIM - BOUNDARY_WIDTH &
    Y..cm. >= BOUNDARY_WIDTH & Y..cm. <= PLOT_DIM - BOUNDARY_WIDTH
  ) %>%
  left_join( # Insert recruitment factor
    rec_fact,
    by = join_by(Year == Year, Taxon == Taxon)
  ) %>%
  left_join( # Insert environmental factors
    env_ann_summ,
    by = join_by(Year == Year)
  ) %>%
  dplyr::select( # Only select variables important for stats
    c(
      # Basic variables
      Taxon, Stage,
      # Reponse variables
      Progression, Retrogression, Survival.future, Growth.future, X..Capitulescences,
      Recruitment.factor,
      # Explanatory variables (+ one covariate)
      ends_with(".crowding"),
      # Covariates
      GS.Precipitation.inch, GS.Avg.temperature.F, GS.Var.temperature.F,
      Length..cm., Length..cm..prev, X..cm., Y..cm.,
      # Random effect
      Plot, Year,
      # Is_recruit to set previous length as zero if recruit
      Is.recruit
    )
  ) %>%
  mutate( # Convert plot & year to factor
    Plot = as.factor(Plot),
    Year = as.factor(Year),
    Length..cm..prev = ifelse(Is.recruit == 1, 0, Length..cm..prev)
  ) %>%
  dplyr::select( # Remove Is.recruit
    -Is.recruit
  )

# Store stats dataset
write_csv(spp_stat_data,
          here("StageDependentComp", "result_data", "wp2", "spp_stat_summ.csv"))

# Subdivide dataset
sp_stage_data <- list(
  sp1_s = spp_stat_data %>%
    filter(Taxon == SP1 & Stage == "S"),
  sp1_l = spp_stat_data %>%
    filter(Taxon == SP1 & Stage == "L"),
  sp2_s = spp_stat_data %>%
    filter(Taxon == SP2 & Stage == "S"),
  sp2_l = spp_stat_data %>%
    filter(Taxon == SP2 & Stage == "L")
)

# Generate summary for each subdivision
sp_stage_data_mean <- lapply(sp_stage_data, function(sp_stage_dat) {
  summarise(sp_stage_dat, across(Growth.future:Length..cm., mean))
}) %>%
  bind_rows(.id = "Category") %>%
  pivot_longer(!Category, names_to = "Variable", values_to = "Mean")
sp_stage_data_sd <- lapply(sp_stage_data, function(sp_stage_dat) {
  summarise(sp_stage_dat, across(Growth.future:Length..cm., sd))
}) %>%
  bind_rows(.id = "Category") %>%
  pivot_longer(!Category, names_to = "Variable", values_to = "Std.dev")
sp_stage_data_summ <- full_join(sp_stage_data_mean, sp_stage_data_sd,
                                by = join_by(Category, Variable))

# Store statistical summary
write_csv(sp_stage_data_summ, here("StageDependentComp", "result_data", "wp2",
                                   "variable_stat_summ.csv"))

# ===== Check distribution of vital rates =====
# 
# plot(fitdist(sp_stage_data$sp1_l$X..Capitulescences, distr = "norm"))
# plot(fitdist(sp_stage_data$sp2_l$X..Capitulescences, distr = "gamma"))

# Testing transformation of growth

# ggplot(sp_stage_data_norm$sp1_s, aes(x = sqrt(Growth.future - min(Growth.future)))) +
#   geom_histogram()
# 
# ggplot(sp_stage_data_norm$sp1_l, aes(x = sqrt(max(Growth.future - min(Growth.future)) - (Growth.future - min(Growth.future))))) +
#   geom_histogram()
# 
# ggplot(sp_stage_data_norm$sp2_s, aes(x = sqrt(Growth.future - min(Growth.future)))) +
#   geom_histogram()
# 
# ggplot(sp_stage_data_norm$sp2_l, aes(x = sqrt(max(Growth.future - min(Growth.future)) - (Growth.future - min(Growth.future))))) +
#   geom_histogram()

# ===== Transform/normalise variables =====

# Normalise continuous explanatory variables
sp_stage_data_norm <- lapply(names(sp_stage_data), function(sp_stage_name) {
  sp_stage_dat <- sp_stage_data[[sp_stage_name]]
  
  temp_dat <- sp_stage_dat %>%
    group_by(Taxon, Stage) %>%
    mutate( # Normalise for mean and std dev
      across(Sp1.crowding:Y..cm., function(x) {
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      })
    )
  
  # Transform growth based on class
  # if(sp_stage_name %in% c("sp1_s", "sp2_s")) { # If small
  #   # Individuals tend to grow; therefore take square root transformation
  #   temp_dat <- temp_dat %>%
  #     mutate(
  #       Growth.future = sqrt(Growth.future - min(Growth.future))
  #     )
  # } else { # If large
  #   # Individuals tend to shrink; therefore reflect then take square root transformation
  #   temp_dat <- temp_dat %>%
  #     mutate(
  #       Growth.future = sqrt(max(Growth.future - min(Growth.future)) - (Growth.future - min(Growth.future)))
  #     )
  # }
  
  temp_dat
})
names(sp_stage_data_norm) <- c("sp1_s", "sp1_l", "sp2_s", "sp2_l")

# ===== Fit GLMMs of varying complexity & document results =====

# # Model types
# stat_model_types <- as.list(c(
#   "M1" = "",
#   "M2" = "Total.crowding",
#   "M3" = "Sp1.crowding + Sp2.crowding + Other.crowding",
#   "M4" = "Sp1.S.crowding + Sp1.L.crowding +
#           Sp2.S.crowding + Sp2.L.crowding +
#           Other.crowding"
# ))
# 
# # Model varieties within each type
# stat_model_templates <- as.list(c(
#   ".0" = "GS.Precipitation.inch + (1 | Plot)",
#   ".L" = "GS.Precipitation.inch + Length..cm. + (1 | Plot)",
#   ".LPREV" = "GS.Precipitation.inch + Length..cm..prev + (1 | Plot)",
#   ".YR" = "(1 | Plot) + (1 | Year)",
#   ".L.YR" = "Length..cm. + (1 | Plot) + (1 | Year)",
#   ".LPREV.YR" = "Length..cm..prev + (1 | Plot) + (1 | Year)"
# ))
# 
# # Full list of models to fit
# stat_model_forms <- lapply(names(stat_model_types), function (type_name) {
#   stat_model_temp <- lapply(stat_model_templates, function (model_template) {
#     
#     # Build model equation string
#     if(type_name == "M1") {
#       model_str_interm <- paste0("~ ", model_template)
#     } else {
#       model_str_interm <- paste0("~ ", stat_model_types[type_name],
#                                  " + ", model_template)
#     }
#     model_str_interm
#   })
#   
#   # Add type prefix to model name
#   names(stat_model_temp) <- paste0(type_name, names(stat_model_temp), sep = "")
#   
#   # Return list
#   stat_model_temp
# })
# 
# # Flatten & convert back to list
# stat_model_forms <- as.list(unlist(stat_model_forms))
# 
# # Remove M1.YR which only includes random effects
# 
# # Create a list of list of model formulae
# stat_model_full_list <- lapply(stat_vrs, function(vr) {
#   lapply(stat_model_forms, function(model_form) {
#     as.formula(paste(vr, model_form))
#   })
# })
# names(stat_model_full_list) <- stat_vrs
# 
# # ===== Fit negative binomial model with random effects =====
# 
# ## ===== Fit models to the data =====
# stat_fitted_nb <- list(
#   "SP1.S" = fit_vr_models(stat_model_full_list, "S", sp_stage_data_norm$sp1_s,
#                                  fec_family = "nb"),
#   "SP1.L" = fit_vr_models(stat_model_full_list, "L", sp_stage_data_norm$sp1_l,
#                                  fec_family = "nb"),
#   "SP2.S" = fit_vr_models(stat_model_full_list, "S", sp_stage_data_norm$sp2_s,
#                                  fec_family = "nb"),
#   "SP2.L" = fit_vr_models(stat_model_full_list, "L", sp_stage_data_norm$sp2_l,
#                                  fec_family = "nb")
# )



## ===== Extract summary values from the fitted models =====

# stat_summ_df_nb <- get_stat_summ_df(stat_fitted_nb) %>%
#   group_by(Taxon, Stage, Vital.rate) %>%
#   arrange(Taxon, Stage, Vital.rate, AIC) # Arrange by AIC for each vital rate
# 
# # Store summary as a csv
# write_csv(stat_summ_df_nb,
#           here("StageDependentComp", "result_data", "wp2",
#                "stat_summary_nb.csv"))
# 
# # Get model with lowest AIC for each stage & vital rate
# stat_min_AIC_nb <- stat_summ_df_nb %>%
#   group_by(Taxon, Stage, Vital.rate) %>%
#   filter(
#     AIC <= min(AIC) + 2 # Select models with less than 2 units difference from minimum AIC
#   )
# 
# # Store summary as a csv
# write_csv(stat_min_AIC_nb,
#           here("StageDependentComp", "result_data", "wp2",
#                "stat_min_AIC_nb.csv"))

# ===== Inspect residuals =====

# plot(residuals(stat_fitted_nb$SP2.L$Survival.future$M1))
# print(which(residuals(stat_fitted_nb$SP2.L$Survival.future$M1.L) < -2))
# 
# vif(stat_fitted_nb$SP2.L$X..Capitulescences$M1.L)
# cor(sp_stage_data$sp2_l %>%
#       dplyr::select(GS.Precipitation.inch, Length..cm.))
# 
# summary(stat_fitted_nb$SP1.S$Progression$M3)

# ===== Calculate proportion of most predictive models across vital rates =====
# NB: Only for stat_fitted_nb

# Determine the simplest model out of those in stat_min_AIC_nb

# stat_bestm_nb <- stat_min_AIC_nb %>%
#   filter(
#     Model.name == min(Model.name)
#   )

# ===== Bayesian model fitting =====

## ===== Test scripts ======

# Investigate default priors

default_prior(X..Capitulescences ~ Sp1.S.crowding + Sp1.L.crowding +
                Sp2.S.crowding + Sp2.L.crowding +
                Other.crowding + Length..cm. + (1 | Plot) + (1 | Year),
              data = sp_stage_data_norm$sp1_l,
              family = negbinomial())

# Set priors

test_model_b <- brm(
  X..Capitulescences ~ Sp1.S.crowding + Sp1.L.crowding +
    Sp2.S.crowding + Sp2.L.crowding +
    Other.crowding + Length..cm. + (1 | Plot) + (1 | Year),
  family = negbinomial(),
  data = sp_stage_data_norm$sp1_l,
  chains = 2, # n of chains
  iter = 5000, # n of iterations, including burnin
  warmup = 1000, # burnin
  thin = 1
)

test_model_b

plot(test_model_b)
# 
# test_model_b_2 <- brm(
#   X..Capitulescences ~ GS.Precipitation.inch + Length..cm. + (1 | Plot),
#   family = negbinomial(),
#   data = sp_stage_data_norm$sp2_l,
#   chains = 2, # nb of chains
#   iter = 5000, # nb of iterations, including burnin
#   warmup = 1000, # burnin
#   thin = 1
# )
# 
# summary(test_model_b_2)
# plot(test_model_b_2)

# ===== Test plots of vital rate against density =====
# 
# ggplot(sp_stage_data$sp1_l, aes(x = Other.crowding, y = X..Capitulescences)) +
#   geom_point() +
#   theme_jun1()
# 
# ggplot(spp_stat_data, aes(x = log(Sp1.S.crowding))) +
#   geom_histogram()

