################################################################################
# Script for calculating matrices for the focal species
#
# Written by Young Jun Lee
# Oct 2024

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {
  
  # Load species vital rate data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_sp_vr.R"))
  
  # Load annual recruitment factor
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_sp_rec_fact.R"))
  
}

# Load relevant functions
source(here("StageDependentComp","code","census_analysis","func","stage_func.R"))
source(here("StageDependentComp","code","census_analysis","func","mat_func.R"))

# ===== Calculate annual vital rates =====

# NB: Here, I am using the number of new tags for a species to directly be the number of new recruits.



# Get annual vital rates
sp1_ann_vrs <- summ_ann_vrs(sp1_data, sp1_rec_fact)
sp2_ann_vrs <- summ_ann_vrs(sp2_data, sp2_rec_fact)

# print(sp1_ann_vrs)
# print(sp2_ann_vrs)

# Save annual vital rates as csv
write_csv(sp1_ann_vrs, here("StageDependentComp", "result_data", "wp2", "sp1_vr_summ.csv"))
write_csv(sp2_ann_vrs, here("StageDependentComp", "result_data", "wp2", "sp2_vr_summ.csv"))

# ===== Calculate annual matrices based on vital rates =====

# Calculate annual matrices
sp1_ann_mats <- calc_ann_mat(sp1_ann_vrs)
sp2_ann_mats <- calc_ann_mat(sp2_ann_vrs)

# Dataframe to store the annual matrices
ann_mat_df <- tibble(
  Taxon = character(),
  Year = numeric(),
  S.S = numeric(),
  S.L = numeric(),
  L.S = numeric(),
  L.L = numeric(),
  Lambda = numeric(),
  Generation.time = numeric()
)

# Unpack annual matrices into dataframe & merge
sp1_ann_mat_df <- summ_ann_mat_stats(sp1_ann_mats)
sp1_ann_mat_df$Taxon <- rep(SP1, nrow(sp1_ann_mat_df))
sp2_ann_mat_df <- summ_ann_mat_stats(sp2_ann_mats)
sp2_ann_mat_df$Taxon <- rep(SP2, nrow(sp2_ann_mat_df))
ann_mat_df <- ann_mat_df %>%
  bind_rows(sp1_ann_mat_df) %>%
  bind_rows(sp2_ann_mat_df)

# Save merged dataframe
write_csv(ann_mat_df, here("StageDependentComp","result_data","wp2","sp_annual_matrices.csv"))

# ====== Calculate pooled matrices across years =======

# Calculate pooled matrices
sp1_pooled_mat <- calc_pooled_mat(sp1_ann_vrs)
sp2_pooled_mat <- calc_pooled_mat(sp2_ann_vrs)

# Calculate statistics from pooled matrices & add taxon info
sp1_pooled_mat_stats <- calc_mat_stats(sp1_pooled_mat) %>%
  mutate(Taxon = SP1) %>%
  select(Taxon, everything())
sp2_pooled_mat_stats <- calc_mat_stats(sp2_pooled_mat) %>%
  mutate(Taxon = SP2) %>%
  select(Taxon, everything())

# Bind dfs and save to csv
pooled_mat_stats <- bind_rows(sp1_pooled_mat_stats, sp2_pooled_mat_stats)
write_csv(pooled_mat_stats, here("StageDependentComp","result_data","wp2","sp_pooled_matrices.csv"))

