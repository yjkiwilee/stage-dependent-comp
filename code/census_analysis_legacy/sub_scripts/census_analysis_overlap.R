################################################################################
#' Script for analysing stage-dependence in biotic interactions
#'
#' Written by Young Jun Lee
#' Apr 2025

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {
  
  if(system.file(package="here") == "") {
    install.packages("here")
    library("here")
  }
  
  library("here")
  
  # Load necessary packages
  source(here("StageDependentComp",
              "code",
              "census_analysis",
              "sub_scripts",
              "census_analysis_load_packages.R"))
  
  source(here("StageDependentComp",
              "code",
              "census_analysis",
              "sub_scripts",
              "census_analysis_load_stat_packages.R"))
  
  # Load processed census data as "census_data"
  source(here("StageDependentComp",
              "code",
              "census_analysis",
              "sub_scripts",
              "census_analysis_load_proc_data.R"))
}

# Load functions used for spatial analysis
source(here("StageDependentComp","code","census_analysis","func","spatial_analysis.R"))

# ====== Initial parameters ======

# Focal species to test stage-dependence in
FOCAL_SPP <- c(
  "Lupinus argenteus",
  "Ivesia gordonii",
  "Eriogonum umbellatum"
  # "Elymus lanceolatus"
)

# Combinations to test stage-dependence in
# FOCAL_SPP_PAIRS <- expand.grid(f_sp = FOCAL_SPP, n_sp = FOCAL_SPP) %>%
#   arrange(f_sp)

FOCAL_SPP_PAIRS <- data.frame(
  f_sp = FOCAL_SPP[3],
  n_sp = FOCAL_SPP[3]
)

FOCAL_SPP_PAIRS

# Vital rates to test stage-dependence in
FOCAL_VRS <- c(
  "Growth.future", # Growth from t to t+1, conditional on survival
  "Survival.future", # Survival from t to t+1
  "Is.flowering", # Flowering status at t, conditional on survival
  "X..Capitulescences" # Number of flowers at t, conditional on flowering
)

# Plot dimension in cm
PLOT_DIM <- 200
# The width of the boundary area; 0 means no individual is excluded from analysis
BOUNDARY_WIDTH <- 0

# Number of cores to use
N_CPUS <- 16

# Number of bootstrap replicates to calculate
BOOTSTRAP_SIZE <- 29999

# Random seed to use for bootstrapping
# This seed is applied for each call of the boot() function
RANDOM_SEED <- 1

# ======= Analyse each combination ======

for(rowi in 1:nrow(FOCAL_SPP_PAIRS)) {
  for(vri in 1:length(FOCAL_VRS)) {
    # Focal species (species receiving the density-dependent effect)
    SP_FOCAL <- FOCAL_SPP_PAIRS$f_sp[rowi]
    # Neighbour species (species exerting the density-dependent effect)
    SP_NEI <- FOCAL_SPP_PAIRS$n_sp[rowi]
    
    SP_FOCAL_ABB <- tolower(substr(SP_FOCAL, 1, 3))
    SP_NEI_ABB <- tolower(substr(SP_NEI, 1, 3))
    
    FOCAL_VR <- FOCAL_VRS[vri]
    
    # Check if file already exists
    spp_fname <- here(
      "StageDependentComp",
      "result_data",
      "wp2",
      "bootstrap",
      sprintf(
        "%s_%s_%s_bootstrap.csv",
        SP_FOCAL_ABB,
        SP_NEI_ABB,
        FOCAL_VR
      )
    )
    
    if(file.exists(spp_fname)) {
      next
    } else {
      # Analyse interaction and get statistical summary
      # source(here("StageDependentComp",
      #             "code",
      #             "census_analysis",
      #             "sub_scripts",
      #             "census_analysis_overlap_single.R"),
      #        local = TRUE)
      source(here("StageDependentComp",
                  "code",
                  "census_analysis",
                  "sub_scripts",
                  "census_analysis_overlap_single_svr.R"),
             local = TRUE)
    }
  }
}

cat(sprintf("%s\n", format(Sys.time(), "%c")))
cat("Done!\n")

