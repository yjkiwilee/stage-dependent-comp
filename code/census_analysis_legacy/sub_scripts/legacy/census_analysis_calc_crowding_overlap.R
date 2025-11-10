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
  plot_sub_data <- calc_crowd_overlap(plot_sub_data, SP1, SP2, sp1_thres_len, sp2_thres_len) %>%
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
