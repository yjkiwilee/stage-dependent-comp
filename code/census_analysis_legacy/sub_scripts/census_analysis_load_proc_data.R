################################################################################
# Script for loading the pre-processed census data
#
# Written by Young Jun Lee
# Oct 2024

# Load census data
census_data <- read_csv(here("StageDependentComp",
                             "result_data",
                             "data_cleanup",
                             "demography_2014-2022_processed.csv"))
