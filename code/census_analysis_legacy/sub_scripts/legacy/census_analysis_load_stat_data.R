################################################################################
#' Script for loading statistical data produced by
#' census_analysis_stat_glmm.R
#' Written by Young Jun Lee
#' Feb 2024


# Read stats dataset
spp_stat_data <- read_csv(here("StageDependentComp", "result_data", "wp2",
               "spp_stat_summ.csv"))

# Read statistical summary
sp_stage_data_summ <- read_csv(here("StageDependentComp", "result_data", "wp2",
                                   "variable_stat_summ.csv"))

# Save transformed data
sp_stage_data_norm_df <- read_csv(here("StageDependentComp", "result_data", "wp2",
               "stat_norm_data.csv")) %>%
  mutate( # Convert plot & year to factor
     Plot = as.factor(Plot),
     Year = as.factor(Year)
   )
sp_stage_data_normtrans_df <- read_csv(here("StageDependentComp", "result_data", "wp2",
                                       "stat_normtrans_data.csv")) %>%
  mutate( # Convert plot & year to factor
    Plot = as.factor(Plot),
    Year = as.factor(Year)
  )

# Subdivide dataset
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
