################################################################################
# Script for setting the focal species & subsetting census data accordingly
#
# Written by Young Jun Lee
# Oct 2024

# Two focal species
# SP1 <- "Heterotheca villosa"
# SP2 <- "Senecio crassulus"
SP1 <- "Lupinus argenteus"
SP2 <- "Eriogonum umbellatum"

# Shorthand for species
SP1_SHORT <- "L.arg"
SP2_SHORT <- "E.umb"

# Subset data
sub_data <- census_data %>%
  filter(
    (Taxon %in% c(SP1, SP2))
  )

# Extract individual species only
sp1_data <- sub_data %>%
  filter(Taxon == SP1)
sp2_data <- sub_data %>%
  filter(Taxon == SP2)