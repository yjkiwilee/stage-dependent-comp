################################################################################
#' Script for fitting a GLMM to the demographic & environmental data
#' This is different from census_analysis_stat_glmm in the specific models;
#' refer to research log on 27022025
#' CRUCIALLY, THIS SCRIPT USES AN EXTERNAL VALUE OF ALPHA
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
  
  # Load annual environmental data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_ann_env_data.R"))
  
  # Load recruitment factor data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_sp_rec_fact.R"))
  
  # Load stat data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_load_stat_data_coeff.R"))
  
  # Load functions for stats
  source(here("StageDependentComp", "code", "census_analysis", "func",
              "stat_fitting.R"))
}

# ===== Fit GLMMs of varying complexity & document results =====

# Model types
stat_model_types <- as.list(c(
  "M1" = "",
  "M2" = "Total.crowding",
  "M3" = "Sp1.crowding + Sp2.crowding + Other.crowding",
  "M4" = "Sp1.S.crowding + Sp1.L.crowding +
          Sp2.S.crowding + Sp2.L.crowding +
          Other.crowding"
))

# Model varieties within each type
stat_model_templates <- as.list(c(
  # ".0" = "GS.Precipitation.inch + (1 | Plot)",
  # ".L" = "GS.Precipitation.inch + Length..cm. + (1 | Plot)",
  # ".LPREV" = "GS.Precipitation.inch + Length..cm..prev + (1 | Plot)",
  # ".YR" = "(1 | Plot) + (1 | Year)",
  ".L.YR" = "Length..cm. + (1 | Plot) + (1 | Year)"
  # ".LPREV.YR" = "Length..cm..prev + (1 | Plot) + (1 | Year)"
))

# Full list of models to fit
stat_model_forms <- lapply(names(stat_model_types), function (type_name) {
  stat_model_temp <- lapply(stat_model_templates, function (model_template) {
    
    # Build model equation string
    if(type_name == "M1") {
      model_str_interm <- paste0("~ ", model_template)
    } else {
      model_str_interm <- paste0("~ ", stat_model_types[type_name],
                                 " + ", model_template)
    }
    model_str_interm
  })
  
  # Add type prefix to model name
  names(stat_model_temp) <- paste0(type_name, names(stat_model_temp), sep = "")
  
  # Return list
  stat_model_temp
})

# Flatten & convert back to list
stat_model_forms <- as.list(unlist(stat_model_forms))

# Create a list of list of model formulae
stat_model_full_list <- lapply(stat_vrs, function(vr) {
  lapply(stat_model_forms, function(model_form) {
    as.formula(paste(vr, model_form))
  })
})
names(stat_model_full_list) <- stat_vrs

# ====== Fit models to data =======

stat_fitted <- list(
  "SP1.S" = fit_vr_lmm(stat_model_full_list, "S", sp_stage_data_norm$sp1_s,
                          maxit = 7e4),
  "SP1.L" = fit_vr_lmm(stat_model_full_list, "L", sp_stage_data_norm$sp1_l,
                          maxit = 7e4),
  "SP2.S" = fit_vr_lmm(stat_model_full_list, "S", sp_stage_data_norm$sp2_s,
                          maxit = 7e4),
  "SP2.L" = fit_vr_lmm(stat_model_full_list, "L", sp_stage_data_norm$sp2_l,
                          maxit = 7e4)
)

## ===== Generate plot of model prediction vs data for each model ======

predictions_df <- tibble()

for(model_class in names(stat_fitted)) {
  for(model_vr in names(stat_fitted[[model_class]])) {
    for(model_name in names(stat_fitted[[model_class]][[model_vr]])) {
      model_temp <- stat_fitted[[model_class]][[model_vr]][[model_name]]
      
      model_species <- if(strsplit(model_class, "\\,")[[1]][1] == "SP1") { SP1 } else { SP2 }
      model_stage <- strsplit(model_class, "\\.")[[1]][2]
      
      # Get predictions
      stat_predictions <- model_temp@frame
      stat_predictions$Predict <- predict(model_temp)
      
      stat_predictions$Vital.rate.value = stat_predictions[[model_vr]]
      stat_predictions[[model_vr]] <- NULL
      stat_predictions <- stat_predictions %>%
        mutate(
          Class = model_class,
          Vital.rate = model_vr,
          Model.name = model_name
        )
      
      predictions_df <- bind_rows(predictions_df, stat_predictions)
    }
  }
}

# Plot as facets
ggplot(predictions_df, aes(x = Vital.rate.value, y = Predict, color = Model.name)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(Class ~ Vital.rate + Model.name, scales = "free",
             dir = "v", nrow = length(stat_model_forms)) +
  scale_color_brewer(palette = "Set2") +
  theme_jun1()

ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
            "crowding_factor_test", paste0("data_predict_lmm_", ALPHA, ".png")), width = 22, height = 11)


## ===== Extract summary values from the fitted models =====

stat_summ_full_df <- get_stat_summ_df(stat_fitted) %>%
  group_by(Taxon, Stage, Vital.rate) %>%
  arrange(Taxon, Stage, Vital.rate, AIC) # Arrange by AIC for each vital rate

### ====== Rescale model coefficients to account for normalisation ======

# Map between class names in variable summary and taxon / stage
class_taxon <- c(
  `sp1_s` = "SP1",
  `sp1_l` = "SP1",
  `sp2_s` = "SP2",
  `sp2_l` = "SP2"
)
class_stage <- c(
  `sp1_s` = "S",
  `sp1_l` = "L",
  `sp2_s` = "S",
  `sp2_l` = "L"
)

# Add column of associated mean and standard deviations
sp_stage_data_summ_proc <- sp_stage_data_summ %>%
  dplyr::mutate(
    Taxon = class_taxon[as.factor(Category)],
    Stage = class_stage[as.factor(Category)]
  ) %>%
  dplyr::select(!Category)

# Rescale all coefficients
stat_summ_rescaled_df <- left_join(stat_summ_full_df, sp_stage_data_summ_proc,
          by = join_by(Taxon, Stage, Variable),
          suffix = c("", ".y")) %>%
  dplyr::mutate(
    Coefficient.r = Coefficient / Std.dev,
    Lower.95CI.r = Lower.95CI / Std.dev,
    Upper.95CI.r = Upper.95CI / Std.dev
  )

# Calculate transformation to be applied to intercepts
stat_summ_temp_df <- stat_summ_rescaled_df %>%
  dplyr::filter(
    Variable != "(Intercept)"
  ) %>%
  dplyr::group_by(Taxon, Stage, Vital.rate) %>%
  dplyr::summarise(
    Coefficient.trans = sum(-(Coefficient * Mean / Std.dev)),
    Lower.95CI.trans = sum(-(Lower.95CI * Mean / Std.dev)),
    Upper.95CI.trans = sum(-(Upper.95CI * Mean / Std.dev))
  ) %>%
  dplyr::mutate(Variable = "(Intercept)")

# Calculate intercepts
stat_summ_rescaled_df <- dplyr::left_join(stat_summ_rescaled_df, stat_summ_temp_df,
                                   by = join_by(Taxon, Stage, Variable),
                                   suffix = c("", ".y")) %>%
  dplyr::mutate(
    Coefficient.r = ifelse(Variable == "(Intercept)",
                     Coefficient + Coefficient.trans,
                     Coefficient.r),
    Lower.95CI.r = ifelse(Variable == "(Intercept)",
                     Lower.95CI + Lower.95CI.trans,
                     Lower.95CI.r),
    Upper.95CI.r = ifelse(Variable == "(Intercept)",
                     Upper.95CI + Upper.95CI.trans,
                     Upper.95CI.r)
  ) %>%
  dplyr::select(!c(Coefficient.trans, Lower.95CI.trans, Upper.95CI.trans))

# Store summary as a csv
write_csv(stat_summ_rescaled_df,
          here("StageDependentComp", "result_data", "wp2",
               "crowding_factor_test",
               paste0("stat_summary_full_lmm_", ALPHA, ".csv")))

