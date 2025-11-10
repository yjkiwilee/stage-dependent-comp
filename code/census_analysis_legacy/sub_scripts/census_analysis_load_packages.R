################################################################################
# Script for loading the packages needed for census data analysis as part of WP2.
#
# Written by Young Jun Lee
# Oct 2024

# Check if pacman is installed, if not, install pacman
if(system.file(package="pacman") == "") {
  install.packages("pacman")
}
# Install & load the rest of the packages
pacman::p_load("popbio", "tidyverse", "ggforce", "here", "patchwork", "parallel",
               "RColorBrewer")
