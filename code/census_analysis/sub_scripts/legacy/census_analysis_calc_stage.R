################################################################################
# Script for determining the stage boundary of focal species & growth / retrogression
# and storing the result
#
# Written by Young Jun Lee
# Oct 2024

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {

  # Load modified census data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))
  
  # Choose focal species & subset data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_set_sp.R"))

}

# test_data <- census_data %>%
#   filter(Taxon == "Senecio crassulus")
# 
# seedling_data <- test_data %>%
#   filter(Is.seedling == 1) %>%
#   select(Year, Tag) %>%
#   rename(Seedling.year = Year)
# 
# test_data <- test_data %>%
#   left_join(seedling_data, by = join_by(Tag)) %>%
#   mutate(
#     Age = Year - Seedling.year
#   )
# 
# ggplot(test_data, aes(x = Age, y = Length..cm.)) +
#   geom_point()

# seedling_dat <- sp2_data %>%
#   filter(Is.seedling == 1 & Survival == 1) 
# 
# ggplot(seedling_dat, aes(x = Length..cm.)) +
#   geom_boxplot()

# ===== Initial setup =====

# Load function needed for stage boundary determination
source(here("StageDependentComp", "code", "census_analysis", "func", "stage_func.R"))

# Target proportion of reproductive individuals in small class
SMALL_REPROD_PROP <- 0

# ===== Determine stage boundary =====

# Calculate stage boundaries
# sp1_thres_len <- calc_thres_len(SMALL_REPROD_PROP, sp1_data %>% filter(Died.this.census.final == 0))
# sp2_thres_len <- calc_thres_len(SMALL_REPROD_PROP, sp2_data %>% filter(Died.this.census.final == 0))
sp1_thres_len <- median(filter(sp1_data, Died.this.census.final == 0)$Length..cm.)
sp2_thres_len <- median(filter(sp1_data, Died.this.census.final == 0)$Length..cm.)

# Count total number of individuals & individuals in the small stage
sp1_n <- nrow(sp1_data %>% filter(Died.this.census.final == 0))
sp2_n <- nrow(sp2_data %>% filter(Died.this.census.final == 0))
sp1_nsmall <- nrow(sp1_data %>% filter(Length..cm. <= sp1_thres_len & Died.this.census.final == 0))
sp2_nsmall <- nrow(sp2_data %>% filter(Length..cm. <= sp2_thres_len & Died.this.census.final == 0))

# Store info in a df
thres_len_df <- tibble(
  Taxon = c(SP1, SP2),
  Stage.boundary..cm = c(sp1_thres_len, sp2_thres_len),
  N.total = c(sp1_n, sp2_n),
  N.small = c(sp1_nsmall, sp2_nsmall)
)

# Calculate additional variables
thres_len_df <- thres_len_df %>%
  mutate(
    N.large = N.total - N.small,
    Proportion.small = N.small / N.total,
    Proportion.large = (N.total - N.small) / N.total
  )

# Store result in csv
write_csv(thres_len_df, here("StageDependentComp", "result_data", "wp2", "sp_stage_struct.csv"))

# Insert stage classification into the dataset
sp1_data$Stage <- ifelse(sp1_data$Length..cm. <= sp1_thres_len, "S", "L")
sp2_data$Stage <- ifelse(sp2_data$Length..cm. <= sp2_thres_len, "S", "L")

sp1_data$Stage.future <- ifelse(sp1_data$Length..cm..next <= sp1_thres_len, "S", "L")
sp2_data$Stage.future <- ifelse(sp2_data$Length..cm..next <= sp2_thres_len, "S", "L")

# ===== Calculate stage-dependent vital rates =====

## ===== Growth =====

sp1_g <- calc_stage_prog(sp1_data, sp1_thres_len)
sp2_g <- calc_stage_prog(sp2_data, sp2_thres_len)
# Incorporate into the data table
sp1_data$Progression <- sp1_g
sp2_data$Progression <- sp2_g

## ===== Retrogression =====

sp1_r <- calc_stage_retrog(sp1_data, sp1_thres_len)
sp2_r <- calc_stage_retrog(sp2_data, sp2_thres_len)
# Incorporate into the data table
sp1_data$Retrogression <- sp1_r
sp2_data$Retrogression <- sp2_r

# Validate that these vital rates are NA for last survey time point
sp1_data %>%
  group_by(Stage, Stage.future, Progression, Retrogression) %>%
  count()
sp2_data %>%
  group_by(Stage, Stage.future, Progression, Retrogression) %>%
  count()

# Save data to csv files
write_csv(sp1_data, here("StageDependentComp","result_data","wp2","sp1_vr.csv"))
write_csv(sp2_data, here("StageDependentComp","result_data","wp2","sp2_vr.csv"))
