################################################################################
# Script for calculating the annual 'recruitment factor' & storing the results
#
# Written by Young Jun Lee
# Oct 2024

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {
  
  # Load modified census data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))
  
  # Load species data with vital rates
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_sp_vr.R"))
  
}

# Load function needed for annual recruitment factor calculation
source(here("StageDependentComp", "code", "census_analysis", "func", "stage_func.R"))

# Calculate recruitment factor
sp1_rec_fact <- calc_rec_fact(sp1_data) %>%
  mutate(Taxon = SP1) %>%
  relocate(Taxon)
sp2_rec_fact <- calc_rec_fact(sp2_data) %>%
  mutate(Taxon = SP2) %>%
  relocate(Taxon)

# Merge data
rec_fact <- bind_rows(sp1_rec_fact, sp2_rec_fact)

# Store recruitment factor data
write_csv(rec_fact, here("StageDependentComp","result_data","wp2","spp_rec_fact.csv"))
