################################################################################
# Script for pre-processing the census data
#
# Written by Young Jun Lee
# Oct 2024

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {

  if(system.file(package="here") == "") {
    install.packages("here")
    library("here")
  }
  
  # Load config parameters
  source(here("StageDependentComp","code","census_analysis","census_analysis_config.R"))
  
  # Load necessary packages
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_packages.R"))
  
  # Load census data as "census_data"
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_data.R"))
  
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_ann_env_data.R"))
}

# ===== Process vital rate data =====

# DF used for determining future survival and length in next year
census_data_future <- census_data %>%
  mutate(
    Year = Year - 1
  ) %>%
  rename(
    Survival.future = Survival,
    Length..cm..next = Length..cm.
  ) %>%
  dplyr::select(Year, Plot, Tag, Survival.future, Length..cm..next)

# DF used for determining the age of individuals
census_data_age <- census_data %>%
  group_by(Tag) %>%
  summarise( # Get the year each individual was first observed
    Init.year = .data[["Year"]][.data[["Year"]] == min(.data[["Year"]])],
    Init.was.recruit = .data[["Is.recruit"]][.data[["Year"]] == min(.data[["Year"]])]
  ) %>%
  filter(!is.na(Init.was.recruit)) %>% # Remove individuals that aren't recruits
  dplyr::select(Tag, Init.year)

# Merge future survival and length
census_data_temp <- census_data %>%
  left_join(census_data_future, by = join_by(Year, Plot, Tag))

# Determine age of each individual
census_data_temp <- census_data_temp %>%
  left_join(census_data_age, by = join_by(Tag)) %>%
  mutate(
    Age = Year - Init.year
  ) %>%
  dplyr::select(!Init.year)

# Add Is.flowering & Growth.future
census_data_temp <- census_data_temp %>%
  mutate(
    Is.flowering = ifelse(X..Capitulescences > 0, 1, 0),
    Growth.future = Length..cm..next - Length..cm.
  )

# ======= Add stage boundaries & determine life stage of each individual ======

# Calculate stage boundaries for each taxon
census_stage_bounds <- census_data %>%
  group_by(Taxon) %>%
  filter(Is.seedling == 1) %>%
  summarise(
    Stage.boundary.cm. = max(Length..cm., na.rm = TRUE)
  )
# census_stage_bounds <- census_data %>%
#   group_by(Taxon) %>%
#   filter(X..Capitulescences > 0) %>%
#   summarise(
#     Stage.boundary.cm. = min(Length..cm., na.rm = TRUE)
#   )

# Merge stage boundaries into census data and insert stage column
census_data_temp <- left_join(census_data_temp, census_stage_bounds,
                            by = join_by(Taxon)) %>%
  mutate( # Insert column for life stage
    Stage = ifelse(Length..cm. <= Stage.boundary.cm. & 
                     Is.flowering == 0, "S", "L"),
    Stage.future = ifelse(Length..cm..next <= Stage.boundary.cm. &
                            Is.flowering == 0, "S", "L")
  ) %>%
  mutate(
    Progression = ifelse(Stage == "S" & Stage.future == "L", 1, 0),
    Retrogression = ifelse(Stage == "L" & Stage.future == "S", 1, 0)
  )

census_data <- census_data_temp

# ===== Get summary of each taxon ======

census_data_summ_tax <- census_data %>%
  group_by(Year, Taxon, Stage) %>%
  summarise(
    N.live = sum(.data$Survival == 1, na.rm = TRUE),
    N.dead = sum(.data$Survival == 0, na.rm = TRUE),
    N.NA = sum(is.na(.data$Survival))
  ) %>%
  arrange(Taxon)

# ====== Add environmental data =====

census_data <- census_data %>%
  left_join(
    env_ann_summ,
    by = join_by(Year)
  )

# ======= Save processed data =======

# Save updated census data
write_csv(census_data, here("StageDependentComp",
                            "result_data",
                            "data_cleanup",
                            "demography_2014-2022_processed.csv"))

# Save stage boundaries of each taxon
write_csv(census_stage_bounds, here("StageDependentComp",
                                    "result_data",
                                    "wp2",
                                    "taxa_stage_bounds.csv"))

# Save summary of taxa
write_csv(census_data_summ_tax, here("StageDependentComp",
                                     "result_data",
                                     "wp2",
                                     "taxa_summary.csv"))

# ===== Clear unnecessary data ========
rm(
  census_data,
  census_data_age,
  census_data_future,
  census_data_summ_tax,
  census_data_temp,
  census_stage_bounds
)

