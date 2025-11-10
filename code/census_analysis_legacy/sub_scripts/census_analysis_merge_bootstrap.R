# Code for merging across bootstrap replicates


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
FOCAL_SPP_PAIRS <- expand.grid(f_sp = FOCAL_SPP, n_sp = FOCAL_SPP) %>%
  arrange(f_sp)

FOCAL_SPP_PAIRS

# Vital rates to test stage-dependence in
FOCAL_VRS <- c(
  "Growth.future", # Growth from t to t+1, conditional on survival
  "Survival.future", # Survival from t to t+1
  "Is.flowering", # Flowering status at t, conditional on survival
  "X..Capitulescences" # Number of flowers at t, conditional on flowering
)

# ====== Merge bootstrap datasets ======

res_df <- NULL

for(rowi in 1:nrow(FOCAL_SPP_PAIRS)) {
  for(vri in 1:length(FOCAL_VRS)) {
    # Focal species (species receiving the density-dependent effect)
    SP_FOCAL <- FOCAL_SPP_PAIRS$f_sp[rowi]
    # Neighbour species (species exerting the density-dependent effect)
    SP_NEI <- FOCAL_SPP_PAIRS$n_sp[rowi]
    
    SP_FOCAL_ABB <- tolower(substr(SP_FOCAL, 1, 3))
    SP_NEI_ABB <- tolower(substr(SP_NEI, 1, 3))
    
    FOCAL_VR <- FOCAL_VRS[vri]
    
    # Load file
    spp_pair_df <- read_csv(here(
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
    ))
    
    # Merge
    if(is.null(res_df)) {
      res_df <- spp_pair_df
    } else {
      res_df <- bind_rows(res_df, spp_pair_df)
    }
  }
}

res_df[res_df$Has.converged == 0,5:9] <- NA

write_csv(
  res_df,
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "all_bootstrap.csv"
  )
)



