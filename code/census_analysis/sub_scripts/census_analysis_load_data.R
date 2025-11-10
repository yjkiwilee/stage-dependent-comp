################################################################################
# Script for loading the raw census data
#
# Written by Young Jun Lee
# Oct 2024

# Load census data
census_data <-  readr::read_csv(here("StageDependentComp", "data", "demography_2014-2022.csv"))