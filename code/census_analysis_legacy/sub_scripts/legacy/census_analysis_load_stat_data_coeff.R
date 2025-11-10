################################################################################
#' Script for loading statistical data produced by
#' census_analysis_stat_glmm.R
#' Written by Young Jun Lee
#' Feb 2024


# Read stats dataset
spp_stat_data <- read_csv(here("StageDependentComp", "result_data", "wp2",
                               "crowding_factor_test",
                               paste0("stat_nonnorm_data_", ALPHA, ".csv")))

# Read statistical summary
sp_stage_data_summ <- read_csv(here("StageDependentComp", "result_data", "wp2",
                                    "crowding_factor_test",
                                    paste0("variable_stat_summ_", ALPHA, ".csv")))

# Read processed data
sp_stage_data_df <- read_csv(here("StageDependentComp", "result_data", "wp2",
                                       "crowding_factor_test",
                                       paste0("stat_nonnorm_data_", ALPHA, ".csv"))) %>%
  mutate( # Convert plot & year to factor
    Plot = as.factor(Plot),
    Year = as.factor(Year)
  )
sp_stage_data_norm_df <- read_csv(here("StageDependentComp", "result_data", "wp2",
                                       "crowding_factor_test",
                                       paste0("stat_norm_data_", ALPHA, ".csv"))) %>%
  mutate( # Convert plot & year to factor
     Plot = as.factor(Plot),
     Year = as.factor(Year)
   )
sp_stage_data_normtrans_df <- read_csv(here("StageDependentComp", "result_data", "wp2",
                                            "crowding_factor_test",
                                            paste0("stat_normtrans_data_", ALPHA, ".csv"))) %>%
  mutate( # Convert plot & year to factor
    Plot = as.factor(Plot),
    Year = as.factor(Year)
  )

# Subdivide dataset
sp_stage_data <- list(
  sp1_s = sp_stage_data_df %>%
    filter(Taxon == SP1 & Stage == "S"),
  sp1_l = sp_stage_data_df %>%
    filter(Taxon == SP1 & Stage == "L"),
  sp2_s = sp_stage_data_df %>%
    filter(Taxon == SP2 & Stage == "S"),
  sp2_l = sp_stage_data_df %>%
    filter(Taxon == SP2 & Stage == "L")
)
sp_stage_data_norm <- list(
  sp1_s = sp_stage_data_norm_df %>%
    filter(Taxon == SP1 & Stage == "S"),
  sp1_l = sp_stage_data_norm_df %>%
    filter(Taxon == SP1 & Stage == "L"),
  sp2_s = sp_stage_data_norm_df %>%
    filter(Taxon == SP2 & Stage == "S"),
  sp2_l = sp_stage_data_norm_df %>%
    filter(Taxon == SP2 & Stage == "L")
)
sp_stage_data_normtrans <- list(
  sp1_s = sp_stage_data_normtrans_df %>%
    filter(Taxon == SP1 & Stage == "S"),
  sp1_l = sp_stage_data_normtrans_df %>%
    filter(Taxon == SP1 & Stage == "L"),
  sp2_s = sp_stage_data_normtrans_df %>%
    filter(Taxon == SP2 & Stage == "S"),
  sp2_l = sp_stage_data_normtrans_df %>%
    filter(Taxon == SP2 & Stage == "L")
)
