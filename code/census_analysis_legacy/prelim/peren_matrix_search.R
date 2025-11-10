# ===== Extract stage-structured matrices on perennial plant species from COMPADRE =====

# Load packages
pacman::p_load("tidyverse", "Rcompadre")

cdb <- cdb_fetch("compadre")

glimpse(cdb@data)

count(cdb@data, CensusType)
