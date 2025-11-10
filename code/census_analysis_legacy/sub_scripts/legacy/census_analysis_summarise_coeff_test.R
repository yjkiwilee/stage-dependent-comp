################################################################################
#' Script for summarising the crowding function coefficient test results
#' 
#' Written by Young Jun Lee
#' Feb 2024
#' 

# ======= Define model type ======

# ========= Load the output files from the coefficient tests ===========

# Get the list of coefficient values for which output file is available
output_file_names <- list.files(
  path = here("StageDependentComp", "result_data", "wp2", "crowding_factor_test"),
  full.names = FALSE
)

# Filter stat_summary only
output_file_names <- output_file_names[
  grepl(paste0("^stat_summary_full_", MODEL_TYPE, "_[0-9]+(\\.[0-9]+)?.csv"),
        output_file_names)
]

# Extract the coefficients
crowd_coeffs <- as.numeric(str_match(output_file_names, "[0-9]+(\\.[0-9]+)?")[,1])
crowd_coeffs

# ======= Extract summary values for each coefficient value ======

# DF to store the summary values
model_summ_df <- tibble()

for(crowd_coeff in crowd_coeffs) {
  # Read the relevant summary file
  model_fit_summ <- read_csv(
    here("StageDependentComp", "result_data", "wp2", "crowding_factor_test",
         paste0("stat_summary_full_", MODEL_TYPE, "_", crowd_coeff, ".csv"))
  ) %>%
    mutate( # Add column containing crowding function coefficient
      Crowding.func.coeff = crowd_coeff
    )
  
  model_summ_df <- bind_rows(model_summ_df, model_fit_summ)
}

# Merge Taxon & Stage into Class
model_summ_df <- model_summ_df %>%
  mutate(
    Class = paste0(Taxon, ".", Stage)
  ) %>%
  dplyr::select(!c(Taxon, Stage)) %>%
  group_by(Class, Vital.rate)

# ====== Plot exploratory plots ========

# AIC ~ crowding coefficient
ggplot(model_summ_df, aes(x = Crowding.func.coeff, y = AIC, color = Model.name)) +
  geom_point() +
  geom_point(data = filter(model_summ_df, Has.converged != TRUE), color = "black") +
  facet_wrap(Class ~ Vital.rate, scales = "free_y") +
  scale_color_brewer(palette = "Set2") +
  theme_jun1() +
  labs(
    x = "Alpha"
  )

ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
            "crowding_factor_test", paste0("AIC_coeff_summ_", MODEL_TYPE, ".png")),
       width = 12, height = 8)

# Random variation ~ crowding coefficient
ggplot(filter(model_summ_df, Random.variance < 50),
       aes(x = Crowding.func.coeff, y = Random.variance, color = Model.name)) +
  geom_point() +
  geom_point(data = filter(model_summ_df, Random.variance < 50 & Has.converged != TRUE),
             color = "black") +
  facet_wrap(Class ~ Vital.rate, scales = "free_y") +
  scale_color_brewer(palette = "Set2") +
  theme_jun1() +
  labs(
    x = "Alpha"
  )

ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
            "crowding_factor_test", paste0("randvar_coeff_summ_", MODEL_TYPE, ".png")),
       width = 12, height = 8)

# # AIC ~ crowding coefficient with restricted x axis
# model_summ_df_sub <- model_summ_df %>%
#   filter(Crowding.func.coeff <= 1)
# 
# ggplot(model_summ_df_sub, aes(x = Crowding.func.coeff, y = AIC, color = Model.name)) +
#   geom_point() +
#   geom_point(data = filter(model_summ_df_sub, Has.converged != TRUE), color = "black") +
#   facet_wrap(Class ~ Vital.rate, scales = "free_y") +
#   theme_jun1()
# 
# ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
#             "crowding_factor_test", paste0("AIC_coeff_summ_zoom_", MODEL_TYPE, ".png")),
#        width = 12, height = 8)

# Coefficient 95%CI ~ crowding coefficient
ggplot(filter(model_summ_df, Model.name == "M4.L.YR"),
       aes(x = Crowding.func.coeff, y = Coefficient, ymin = Lower.95CI,
           ymax = Upper.95CI, color = Has.converged)) +
  geom_point() +
  geom_errorbar() +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) + # Add hline at zero
  facet_wrap(Class ~ Vital.rate + Variable, scales = "free_y") +
  theme_jun1()

ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
            "crowding_factor_test", paste0("coeff_CI_summ_", MODEL_TYPE, ".png")), width = 50, height = 50, scale = 0.5)

# Coefficient 95%CI ~ crowding coefficient with removal of outliers
ggplot(filter(model_summ_df,
              Model.name == "M4.L.YR" &
              (Coefficient >= (median(Coefficient) - 1.5 * IQR(Coefficient)) &
                Coefficient <= (median(Coefficient) + 1.5 * IQR(Coefficient)))
              ),
       aes(x = Crowding.func.coeff, y = Coefficient, ymin = Lower.95CI, ymax = Upper.95CI, color = Has.converged)) +
  geom_point() +
  geom_errorbar() +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) + # Add hline at zero
  facet_wrap(Class ~ Vital.rate + Variable, scales = "free_y") +
  theme_jun1()

ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
            "crowding_factor_test", paste0("coeff_CI_summ_nooutliers_", MODEL_TYPE, ".png")),
       width = 50, height = 50, scale = 0.5)

# MSE ~ crowding coefficient for M2

ggplot(filter(model_summ_df, Model.name == "M2.L.YR"),
       aes(x = Crowding.func.coeff, y = MSE, color = Has.converged)) +
  geom_point() +
  facet_wrap(Class ~ Vital.rate, scales = "free_y") +
  theme_jun1()

ggsave(here("StageDependentComp", "figures", "census_analysis", "crowding",
            "crowding_factor_test", paste0("MSE_coeff_summ_", MODEL_TYPE, ".png")),
       width = 12, height = 8)
