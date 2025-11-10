################################################################################
# Script for loading stats packages
#
# Written by Young Jun Lee
# Nov 2024

# Check if pacman is installed, if not, install pacman
if(system.file(package="pacman") == "") {
  install.packages("pacman")
}
# Install & load the rest of the packages
# pacman::p_load("tidyverse", "here", "arm", "lme4", "fitdistrplus", "lmtest",
#                "brms", "car", "nlme", "posterior", "tidybayes", "rstan", "rstanarm",
#                "bayesplot")
pacman::p_load("tidyverse", "here", "lme4", "car", "lmtest", "boot")
