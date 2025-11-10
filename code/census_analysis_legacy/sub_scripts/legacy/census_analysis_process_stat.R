################################################################################
#' Script for processing the data from calc_crowding for statistical analysis
#' 
#' Written by Young Jun Lee
#' Feb 2024

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
}

# ===== Import & clean data =====

# Filter focal species vital rate data to only include live individuals with
# Non-NA vital rates
spp_stat_data <- focal_spp_data %>%
  filter(
    Died.this.census.final == 0 &
      !is.na(Survival.future) &
      !is.na(Growth.future) &
      !is.na(X..Capitulescences) &
      !is.na(Stage)
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
          here("StageDependentComp", "result_data", "wp2",
               "spp_stat_summ.csv"))

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
  summarise(sp_stage_dat, across(Progression:Y..cm., mean(.x, na.rm = TRUE)))
}) %>%
  bind_rows(.id = "Category") %>%
  pivot_longer(!Category, names_to = "Variable", values_to = "Mean")
sp_stage_data_sd <- lapply(sp_stage_data, function(sp_stage_dat) {
  summarise(sp_stage_dat, across(Progression:Y..cm., sd(.x, na.rm = TRUE)))
}) %>%
  bind_rows(.id = "Category") %>%
  pivot_longer(!Category, names_to = "Variable", values_to = "Std.dev")
sp_stage_data_summ <- full_join(sp_stage_data_mean, sp_stage_data_sd,
                                by = join_by(Category, Variable))

# Store statistical summary
write_csv(sp_stage_data_summ, here("StageDependentComp", "result_data", "wp2",
                                   "variable_stat_summ.csv"))

# ===== Transform/normalise variables =====

# Save completely untransformed data
write_csv(reduce(sp_stage_data, full_join),
          here("StageDependentComp", "result_data", "wp2",
               "stat_nonnorm_data.csv"))

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
  
  temp_dat
})

names(sp_stage_data_norm) <- c("sp1_s", "sp1_l", "sp2_s", "sp2_l")

# Save transformed data
write_csv(reduce(sp_stage_data_norm, full_join),
          here("StageDependentComp", "result_data", "wp2",
               "stat_norm_data.csv"))

# Normalise & transform continuous explanatory variables
sp_stage_data_normtrans <- lapply(names(sp_stage_data), function(sp_stage_name) {
  sp_stage_dat <- sp_stage_data[[sp_stage_name]]
  
  temp_dat <- sp_stage_dat %>%
    group_by(Taxon, Stage) %>%
    mutate( # Normalise for mean and std dev
      across(Sp1.crowding:Y..cm., function(x) {
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
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
        Growth.future = sqrt(max(Growth.future - min(Growth.future)) -
                               (Growth.future - min(Growth.future)))
      )
  }
  
  temp_dat
})

names(sp_stage_data_normtrans) <- c("sp1_s", "sp1_l", "sp2_s", "sp2_l")

# Save transformed data
write_csv(reduce(sp_stage_data_normtrans, full_join),
          here("StageDependentComp", "result_data", "wp2",
                                   "stat_normtrans_data.csv"))

