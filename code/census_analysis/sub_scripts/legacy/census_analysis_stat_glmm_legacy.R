################################################################################
# Script for fitting a GLMM to the demographic & environmental data
#
# Written by Young Jun Lee
# Nov 2024
# https://cran.r-project.org/web//packages//glmm/vignettes/vignettes.pdf


# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {
  # Load stats packages
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_stat_packages.R"))
  # Load focal species vital rate data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_full_sp_data.R"))
  
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
# In cm
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
      Year, Taxon, Stage,
      # Reponse variables
      Progression, Retrogression, Survival.future, Growth.future, X..Capitulescences,
      Recruitment.factor,
      # Explanatory variables (+ one covariate)
      ends_with(".crowding"),
      # Covariates
      GS.Precipitation.inch, GS.Avg.temperature.F, GS.Var.temperature.F,
      Length..cm., X..cm., Y..cm.,
      # Random effect
      Plot,
      # Extra info
      Year
    )
  ) %>%
  mutate( # Convert plot to factor
    Plot = as.factor(Plot)
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

ggplot(sp_stage_data_norm$sp1_s, aes(x = sqrt(Growth.future - min(Growth.future)))) +
  geom_histogram()

ggplot(sp_stage_data_norm$sp1_l, aes(x = sqrt(max(Growth.future - min(Growth.future)) - (Growth.future - min(Growth.future))))) +
  geom_histogram()

ggplot(sp_stage_data_norm$sp2_s, aes(x = sqrt(Growth.future - min(Growth.future)))) +
  geom_histogram()

ggplot(sp_stage_data_norm$sp2_l, aes(x = sqrt(max(Growth.future - min(Growth.future)) - (Growth.future - min(Growth.future))))) +
  geom_histogram()

# ===== Transform/normalise variables =====

# Normalise continuous explanatory variables
sp_stage_data_norm <- lapply(names(sp_stage_data), function(sp_stage_name) {
  sp_stage_dat <- sp_stage_data[[sp_stage_name]]
  
  temp_dat <- sp_stage_dat %>%
    group_by(Taxon, Stage) %>%
    mutate( # Normalise for mean and std dev
      across(Sp1.crowding:Y..cm., function(x) {
        (x - mean(x)) / sd(x)
      })
    )
  
  # Transform growth based on class
  if(sp_stage_name %in% c("sp1_s", "sp2_s")) { # If small
    # Individuals tend to grow; therefore take square root transformation
    temp_dat <- temp_dat %>%
      mutate(
        Growth.future = sqrt(Growth.future - min(Growth.future))
      )
  } else { # If large
    # Individuals tend to shrink; therefore reflect then take square root transformation
    temp_dat <- temp_dat %>%
      mutate(
        Growth.future = sqrt(max(Growth.future - min(Growth.future)) - (Growth.future - min(Growth.future)))
      )
  }
})
names(sp_stage_data_norm) <- c("sp1_s", "sp1_l", "sp2_s", "sp2_l")

# ===== Check distribution across plots =====

# Summarise distribution across plots
spp_plot_distrib <- spp_stat_data %>%
  group_by(Taxon, Stage, Plot) %>%
  summarise(
    n = n(),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  group_by(Taxon, Stage) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  mutate( # Add dummy order variable for plotting
    N.Order = 1:n()
  )

spp_plot_distrib

# Plot
plt_plot_distrib <- ggplot(spp_plot_distrib, aes(x = N.Order, y = n)) +
  geom_col() +
  facet_wrap(~ Taxon + Stage, scales = "free_y", drop = TRUE) +
  scale_x_discrete(labels = spp_plot_distrib$Plot) +
  theme_jun1() +
  labs(
    x = "Plot",
    y = "Count",
    title = "Distribution of individuals across plots"
  )

plt_plot_distrib

# Save plot
ggsave(here("StageDependentComp", "figures", "census_analysis", "focal_species",
            "spp_plot_distrib.png"), plt_plot_distrib,
       width = FIG_W, height = FIG_H)

# ===== Check spatial autocorrelation =====

# Vital rates of interest
stat_vrs <- c(
  "Progression",
  # "Retrogression", # Excluded due to insufficient instances of retrogression
  "Growth.future",
  "Survival.future",
  "X..Capitulescences"
)

test_spatial_corr <- NULL
test_data_sub <- sp_stage_data$sp1_l
test_data_plots <- unique(test_data_sub$Plot)
for(plot_id in test_data_plots) {
  # Subset data
  plot_data_temp <- test_data_sub %>%
    filter(Plot == plot_id)
  
  # Iterate through year
  test_data_years <- unique(plot_data_temp$Year)
  for(yr in test_data_years) {
    plot_data_yr <- plot_data_temp %>%
      filter(Year == yr)
    
    # Skip if there's only one individual in the plot
    if(nrow(plot_data_yr) <= 1) { next }
    
    # Iterate through all individuals
    for(i in 1:nrow(plot_data_yr)) {
      # Skip if it's the last indiv
      if(i == nrow(plot_data_yr)) {next}
      for(j in (i+1):nrow(plot_data_yr)) {
        # Euclidian distance between individuals
        indiv_dist <- sqrt((plot_data_yr[i,]$X..cm. - plot_data_yr[j,]$X..cm.)^2 +
                             (plot_data_yr[i,]$Y..cm. - plot_data_yr[j,]$Y..cm.)^2)
        
        diff_vr <- lapply(stat_vrs, function(vr) {
          abs(plot_data_yr[[i,vr]] - plot_data_yr[[j,vr]])
        })
        names(diff_vr) <- stat_vrs
        
        test_spat_row <- tibble(
          Year = yr,
          Plot = plot_id,
          Distance = indiv_dist
        ) %>%
          bind_cols(as.data.frame(diff_vr))
        
        if(is.null(test_spatial_corr)) {
          test_spatial_corr <- test_spat_row
        } else {
          test_spatial_corr <- bind_rows(test_spatial_corr, test_spat_row)
        }
      }
    }
  }
}

# Pivot longer
test_spatial_corr <- test_spatial_corr %>%
  pivot_longer(4:ncol(test_spatial_corr), names_to = "Vital.rate", values_to = "Vital.rate.diff")

# Plot with smoothing lines
plt_spat_autocorr <- ggplot(test_spatial_corr, aes(x = Distance, y = Vital.rate.diff)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ Vital.rate, scales = "free_y") +
  theme_jun1() +
  labs(
    x = "Euclidian distance (cm)",
    y = "Difference in vital rate",
    title = "Pairwise difference in vital rates"
  )

plt_spat_autocorr

# Save plot
ggsave(here("StageDependentComp", "figures", "census_analysis", "stat", "spat_autocorr.png"),
       width = FIG_W, height = FIG_H)

# ===== Fit GLMMs of varying complexity & document results =====

# Model varieties
stat_model_forms <- as.list(c(
  "M1" = "~ GS.Precipitation.inch + (1 | Plot)",
  # "M1.L" = "~ GS.Precipitation.inch + Length..cm. + (1 | Plot)",
  "M2" = "~ Total.crowding + GS.Precipitation.inch + (1 | Plot)",
  # "M2.L" = "~ Total.crowding + GS.Precipitation.inch + Length..cm. + (1 | Plot)",
  "M3" = "~ Sp1.crowding + Sp2.crowding + Other.crowding +
          GS.Precipitation.inch + (1 | Plot)",
  # "M3.L" = "~ Sp1.crowding + Sp2.crowding + Other.crowding +
  #         GS.Precipitation.inch + Length..cm. + (1 | Plot)",
  "M4" = "~ Sp1.S.crowding + Sp1.L.crowding +
          Sp2.S.crowding + Sp2.L.crowding +
          Other.crowding +
          GS.Precipitation.inch + (1 | Plot)"
  # "M4.L" = "~ Sp1.S.crowding + Sp1.L.crowding +
  #         Sp2.S.crowding + Sp2.L.crowding +
  #         Other.crowding +
  #         GS.Precipitation.inch + Length..cm. + (1 | Plot)"
))

# Create a list of list of model formulae
stat_model_full_list <- lapply(stat_vrs, function(vr) {
  lapply(stat_model_forms, function(model_form) {
    as.formula(paste(vr, model_form))
  })
})
names(stat_model_full_list) <- stat_vrs

## ===== Fit models to the data =====
stat_fitted <- list(
  "SP1.S" = fit_vr_models(stat_model_full_list, "S", sp_stage_data_norm$sp1_s),
  "SP1.L" = fit_vr_models(stat_model_full_list, "L", sp_stage_data_norm$sp1_l),
  "SP2.S" = fit_vr_models(stat_model_full_list, "S", sp_stage_data_norm$sp2_s),
  "SP2.L" = fit_vr_models(stat_model_full_list, "L", sp_stage_data_norm$sp2_l)
)

## ===== Extract summary values from the fitted models =====

stat_summ_df <- get_stat_summ_df(stat_fitted) %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  arrange(Taxon, Stage, Vital.rate, AIC) # Arrange by AIC for each vital rate

# Store summary as a csv
write_csv(stat_summ_df,
          here("StageDependentComp", "result_data", "wp2",
               "stat_summary.csv"))

# Get model with lowest AIC for each stage & vital rate
stat_min_AIC <- stat_summ_df %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  filter(
    AIC <= min(AIC) + 2 # Select models with less than 2 units difference from minimum AIC
  )

# Store summary as a csv
write_csv(stat_min_AIC,
          here("StageDependentComp", "result_data", "wp2",
               "stat_min_AIC.csv"))

# fixef()
# VarCorr()
# confint.merMod()
# 
# ggplot(sp_stage_data$sp1_l, aes(x = Length..cm., y = Sp1.S.crowding)) +
#   geom_point()

# ===== Test fitting without random effect of plot =====

stat_model_forms_worand <- as.list(c(
  "M1" = "~ GS.Precipitation.inch",
  # "M1.L" = "~ GS.Precipitation.inch + Length..cm.",
  "M2" = "~ Total.crowding + GS.Precipitation.inch",
  # "M2.L" = "~ Total.crowding + GS.Precipitation.inch + Length..cm.",
  "M3" = "~ Sp1.crowding + Sp2.crowding + Other.crowding +
          GS.Precipitation.inch",
  # "M3.L" = "~ Sp1.crowding + Sp2.crowding + Other.crowding +
  #         GS.Precipitation.inch + Length..cm.",
  "M4" = "~ Sp1.S.crowding + Sp1.L.crowding +
          Sp2.S.crowding + Sp2.L.crowding +
          Other.crowding +
          GS.Precipitation.inch"
  # "M4.L" = "~ Sp1.S.crowding + Sp1.L.crowding +
  #         Sp2.S.crowding + Sp2.L.crowding +
  #         Other.crowding +
  #         GS.Precipitation.inch + Length..cm."
))

# Create a list of list of model formulae
stat_model_full_list_worand <- lapply(stat_vrs, function(vr) {
  lapply(stat_model_forms_worand, function(model_form) {
    as.formula(paste(vr, model_form))
  })
})
names(stat_model_full_list_worand) <- stat_vrs

## ===== Fit models to the data =====
stat_fitted_worand <- list(
  "SP1.S" = fit_vr_models_worand(stat_model_full_list_worand, "S", sp_stage_data_norm$sp1_s),
  "SP1.L" = fit_vr_models_worand(stat_model_full_list_worand, "L", sp_stage_data_norm$sp1_l),
  "SP2.S" = fit_vr_models_worand(stat_model_full_list_worand, "S", sp_stage_data_norm$sp2_s),
  "SP2.L" = fit_vr_models_worand(stat_model_full_list_worand, "L", sp_stage_data_norm$sp2_l)
)

## ===== Extract summary values from the fitted models =====

stat_summ_df_worand <- get_stat_summ_df_worand(stat_fitted_worand) %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  arrange(Taxon, Stage, Vital.rate, AIC) # Arrange by AIC for each vital rate

# Store summary as a csv
write_csv(stat_summ_df_worand,
          here("StageDependentComp", "result_data", "wp2",
               "stat_summary_worand.csv"))

# Get model with lowest AIC for each stage & vital rate
stat_min_AIC_worand <- stat_summ_df_worand %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  filter(
    AIC <= min(AIC) + 2 # Select models with less than 2 units difference from minimum AIC
  )

# Store summary as a csv
write_csv(stat_min_AIC_worand,
          here("StageDependentComp", "result_data", "wp2",
               "stat_min_AIC_worand.csv"))

# ===== Fit negative binomial model without random effect =====

## ===== Fit models to the data =====
stat_fitted_worand_nb <- list(
  "SP1.S" = fit_vr_models_worand(stat_model_full_list_worand, "S", sp_stage_data_norm$sp1_s,
                                 fec_family = "nb"),
  "SP1.L" = fit_vr_models_worand(stat_model_full_list_worand, "L", sp_stage_data_norm$sp1_l,
                                 fec_family = "nb"),
  "SP2.S" = fit_vr_models_worand(stat_model_full_list_worand, "S", sp_stage_data_norm$sp2_s,
                                 fec_family = "nb"),
  "SP2.L" = fit_vr_models_worand(stat_model_full_list_worand, "L", sp_stage_data_norm$sp2_l,
                                 fec_family = "nb")
)

## ===== Extract summary values from the fitted models =====

stat_summ_df_worand_nb <- get_stat_summ_df_worand(stat_fitted_worand_nb) %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  arrange(Taxon, Stage, Vital.rate, AIC) # Arrange by AIC for each vital rate

# Store summary as a csv
write_csv(stat_summ_df_worand_nb,
          here("StageDependentComp", "result_data", "wp2",
               "stat_summary_worand_nb.csv"))

# Get model with lowest AIC for each stage & vital rate
stat_min_AIC_worand_nb <- stat_summ_df_worand_nb %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  filter(
    AIC <= min(AIC) + 2 # Select models with less than 2 units difference from minimum AIC
  )

# Store summary as a csv
write_csv(stat_min_AIC_worand_nb,
          here("StageDependentComp", "result_data", "wp2",
               "stat_min_AIC_worand_nb.csv"))

# ===== Fit negative binomial model with random effects =====

## ===== Fit models to the data =====
stat_fitted_nb <- list(
  "SP1.S" = fit_vr_models(stat_model_full_list, "S", sp_stage_data_norm$sp1_s,
                                 fec_family = "nb"),
  "SP1.L" = fit_vr_models(stat_model_full_list, "L", sp_stage_data_norm$sp1_l,
                                 fec_family = "nb"),
  "SP2.S" = fit_vr_models(stat_model_full_list, "S", sp_stage_data_norm$sp2_s,
                                 fec_family = "nb"),
  "SP2.L" = fit_vr_models(stat_model_full_list, "L", sp_stage_data_norm$sp2_l,
                                 fec_family = "nb")
)



## ===== Extract summary values from the fitted models =====

stat_summ_df_nb <- get_stat_summ_df(stat_fitted_nb) %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  arrange(Taxon, Stage, Vital.rate, AIC) # Arrange by AIC for each vital rate

# Store summary as a csv
write_csv(stat_summ_df_nb,
          here("StageDependentComp", "result_data", "wp2",
               "stat_summary_nb.csv"))

# Get model with lowest AIC for each stage & vital rate
stat_min_AIC_nb <- stat_summ_df_nb %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  filter(
    AIC <= min(AIC) + 2 # Select models with less than 2 units difference from minimum AIC
  )

# Store summary as a csv
write_csv(stat_min_AIC_nb,
          here("StageDependentComp", "result_data", "wp2",
               "stat_min_AIC_nb.csv"))

# ===== Inspect residuals =====

plot(residuals(stat_fitted_nb$SP2.L$Survival.future$M1))
print(which(residuals(stat_fitted_nb$SP2.L$Survival.future$M1.L) < -2))

vif(stat_fitted_nb$SP2.L$X..Capitulescences$M1.L)
cor(sp_stage_data$sp2_l %>%
      dplyr::select(GS.Precipitation.inch, Length..cm.))

summary(stat_fitted_nb$SP1.S$Progression$M3)

# ===== Calculate proportion of most predictive models across vital rates =====
# NB: Only for stat_fitted_nb

# Determine the simplest model out of those in stat_min_AIC_nb

stat_bestm_nb <- stat_min_AIC_nb %>%
  filter(
    Model.name == min(Model.name)
  )

# ===== Bayesian model fitting =====

# test_model_b <- brm(
#   X..Capitulescences ~ Sp1.S.crowding + Sp1.L.crowding +
#     Sp2.S.crowding + Sp2.L.crowding +
#     Other.crowding +
#     GS.Precipitation.inch + Length..cm. + (1 | Plot),
#   family = negbinomial(),
#   data = sp_stage_data_norm$sp1_l,
#   chains = 2, # nb of chains
#   iter = 5000, # nb of iterations, including burnin
#   warmup = 1000, # burnin
#   thin = 1
# )
# 
# plot(test_model_b)
# 
# stat_fitted_b_nb <- list(
#   "SP1.S" = fit_vr_models_b(stat_model_full_list, "S", sp_stage_data_norm$sp1_s,
#                           fec_family = "nb"),
#   "SP1.L" = fit_vr_models_b(stat_model_full_list, "L", sp_stage_data_norm$sp1_l,
#                           fec_family = "nb"),
#   "SP2.S" = fit_vr_models_b(stat_model_full_list, "S", sp_stage_data_norm$sp2_s,
#                           fec_family = "nb"),
#   "SP2.L" = fit_vr_models_b(stat_model_full_list, "L", sp_stage_data_norm$sp2_l,
#                           fec_family = "nb")
# )
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

