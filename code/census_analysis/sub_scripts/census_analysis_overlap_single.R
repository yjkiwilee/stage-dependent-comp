################################################################################
#' Script for analysing stage-dependence in biotic interactions
#' for a single species pair (heterospecific or conspecific)
#' To be called within census_analysis_overlap.R
#'
#' Written by Young Jun Lee
#' Apr 2025

cat(sprintf("%s\n", format(Sys.time(), "%c")))
cat(sprintf("Testing effect of %s on %s\n", SP_NEI, SP_FOCAL))

# ====== Extract & process focal data ======

# Subset census data to only include focal species
focal_data <- census_data %>%
  filter(
    Taxon == SP_FOCAL &
      !is.na(Length..cm.) &
      !is.na(X..cm.) &
      !is.na(Y..cm.) &
      X..cm. >= BOUNDARY_WIDTH & X..cm. <= PLOT_DIM - BOUNDARY_WIDTH &
      Y..cm. >= BOUNDARY_WIDTH & Y..cm. <= PLOT_DIM - BOUNDARY_WIDTH
  )

# ====== Extract & process focal & other neighbour data =======

# Focal neighbour data, excluding dormant or dead individuals (i.e. Length = 0)
nei_data <- census_data %>%
  filter(
    Taxon == SP_NEI &
      !is.na(Length..cm.) &
      Length..cm. > 0 &
      !is.na(X..cm.) &
      !is.na(Y..cm.)
  )

# Subset into small & large individuals
nei_st <- list(
  S = nei_data %>% filter(Stage == "S"),
  L = nei_data %>% filter(Stage == "L")
)

# Other neighbour data, excluding dormant or dead individuals (i.e. Length = 0)
other_data <- census_data %>%
  filter(
    Taxon != SP_NEI &
      !is.na(Length..cm.) &
      Length..cm. > 0 &
      !is.na(X..cm.) &
      !is.na(Y..cm.)
  )

# ====== Calculate local cover density =====

focal_data_stat <- NULL

# Check if file already exists
dens_dat_fname <- here(
  "StageDependentComp",
  "result_data",
  "wp2",
  "density",
  sprintf(
    "%s_%s_density.csv",
    SP_FOCAL_ABB,
    SP_NEI_ABB
  )
)

test_df <- read_csv(here(
  "StageDependentComp",
  "result_data",
  "wp2",
  "density",
  "lup_eri_density.csv"
))

test_df %>% group_by(Stage) %>%
  summarise(
    n_zero_small = sum(Overlap.small == 0, na.rm = TRUE),
    n_nonzero_small = sum(Overlap.small != 0, na.rm = TRUE),
    n_zero_large = sum(Overlap.large == 0, na.rm = TRUE),
    n_nonzero_large = sum(Overlap.large != 0, na.rm = TRUE)
  )

# If exists, load data
if(file.exists(dens_dat_fname)) {
  cat(sprintf("%s\n", format(Sys.time(), "%c")))
  cat("Loading density data\n")
  
  focal_data_stat <- read_csv(dens_dat_fname)
} else { # Calculate if doesn't exist
  cat(sprintf("%s\n", format(Sys.time(), "%c")))
  cat("Calculating focal neighbour densities\n")
  
  # Get pairwise local cover densities
  nei_pairs <- lapply(nei_st, function(nei_dat) {
    calc_dens_overlap_pairs(focal_data, nei_dat)
  })
  
  cat(sprintf("%s\n", format(Sys.time(), "%c")))
  cat("Calculating non-focal neighbour densities\n")
  
  other_pairs <- calc_dens_overlap_pairs(focal_data, other_data)
  
  # Add total local cover density columns to focal data
  focal_data_stat <- focal_data %>%
    mutate(
      Overlap.small = calc_overlap_tot_frompair(nei_pairs$S, focal_data),
      Overlap.large = calc_overlap_tot_frompair(nei_pairs$L, focal_data),
      Overlap.other = calc_overlap_tot_frompair(other_pairs, focal_data)
    ) %>%
    mutate(
      Overlap.focal = Overlap.small + Overlap.large,
      Overlap = Overlap.small + Overlap.large + Overlap.other
    )
  
  # z-transform length & overlap variables
  norm_vars <- c("Length..cm.",
                 "Overlap.focal", "Overlap.small", "Overlap.large",
                 "Overlap.other", "Overlap")
  
  for(norm_var in norm_vars) {
    focal_data_stat[[paste0(norm_var, ".norm")]] <-
      (focal_data_stat[[norm_var]] - mean(focal_data_stat[[norm_var]])) /
      sd(focal_data_stat[[norm_var]])
  }
  
  # Save data
  write_csv(focal_data_stat, dens_dat_fname)
}

# ====== Bootstrap GLMM for each vital rate =====

# DF to store boostrap data
glmm_boot_df <- tibble(
  Focal.taxon = character(),
  Neighbour.taxon = character(),
  Focal.vital.rate = character(),
  Focal.stage = character(),
  Intercept = numeric(),
  Length.cm.slope = numeric(),
  Overlap.small.slope = numeric(),
  Overlap.large.slope = numeric(),
  Overlap.other.slope = numeric(),
  Error.occurred = numeric(),
  Has.converged = numeric(),
  Is.singular = numeric(),
  Sample.size.unique = numeric(),
  Sample.size = numeric(),
  Global.sample.size = numeric()
)

# For each vital rate
for(FOCAL_VR in FOCAL_VRS) {
  cat(sprintf("\t%s\n", format(Sys.time(), "%c")))
  cat(sprintf("\tTesting effect on %s\n", FOCAL_VR))
  
  # Subset focal data depending on the type of vital rate being analysed
  focal_data_vr <- focal_data_stat %>%
    filter(
      !is.na(.data[[FOCAL_VR]])
    )
  
  if(FOCAL_VR %in% c("Growth")) {
    # Growth is conditional upon survival
    focal_data_vr <- focal_data_vr %>%
      filter(
        Survival == 1
      )
  } else if(FOCAL_VR == "Is.flowering") {
    # Flowering status is conditional upon survival and non-seedling status
    focal_data_vr <- focal_data_vr %>%
      filter(
          Survival == 1 &
          Stage == "L"
      )
  } else if(FOCAL_VR == "X..Capitulescences") {
    # Fecundity is conditional upon survival and flowering status
    focal_data_vr <- focal_data_vr %>%
      filter(
        Is.flowering == 1 &
          Survival == 1 &
          Stage == "L"
      )
  }
  
  # Model equation
  model_eq <- paste0(FOCAL_VR, " ~ Length..cm..norm +
    Overlap.small.norm +
    Overlap.large.norm +
    Overlap.other.norm + 
    (1 | Year) + (1 | Plot)")
  # Convert to R equation
  model_eq <- as.formula(model_eq)
  
  # Partition focal data into different life stages
  focal_data_st <- list(
    S = focal_data_vr %>% filter(Stage == "S"),
    L = focal_data_vr %>% filter(Stage == "L")
  )
  
  # For each life stage
  for(FOCAL_STAGE in names(focal_data_st)) {
    # Skip flowering probability and fecundity for small individuals
    if(FOCAL_STAGE == "S" && FOCAL_VR %in% c("Is.flowering", "X..Capitulescences")) {
      next
    }
    
    cat(sprintf("\t\t%s\n", format(Sys.time(), "%c")))
    cat(sprintf("\t\tTesting effect on stage %s\n", FOCAL_STAGE))
    
    focal_data_sub <- focal_data_st[[FOCAL_STAGE]]
    
    # z-transform length & overlap variables
    norm_vars <- c("Length..cm.",
                   "Overlap.focal", "Overlap.small", "Overlap.large",
                   "Overlap.other", "Overlap")
    
    for(norm_var in norm_vars) {
      focal_data_sub[[paste0(norm_var, ".norm")]] <-
        (focal_data_sub[[norm_var]] - mean(focal_data_sub[[norm_var]])) /
        sd(focal_data_sub[[norm_var]])
    }
    
    # Bootstrap
    boot_res <- boot(
      focal_data_sub,
      function(df, idx) {
        sub_df <- df[idx,]
        
        # Fit model
        stat_model <- tryCatch(
          {
            if(FOCAL_VR %in% c("Growth.future", "Growth", "X..Capitulescences")) {
              lmer(
                model_eq,
                data = sub_df,
                control = lmerControl(
                  check.rankX = "stop.deficient"
                )
              )
            } else if(FOCAL_VR %in% c("Survival.future", "Survival", "Is.flowering")) {
              glmer(
                model_eq,
                family = "binomial",
                data = sub_df,
                control = glmerControl(
                  optimizer = "bobyqa",
                  check.rankX = "stop.deficient"
                )
              )
            }
          },
          error = function(e) {
            cat(paste0(e, "\n"))
            return(NULL)
          }
        )
        
        # If error occurred, skip
        if(is.null(stat_model)) {
          return(c(
            0, 0, 0, 0, 0,
            1, 0, 0, # 1 indicates error occurred
            length(unique(idx)), length(idx)
          ))
        } else {
          # Get fixed effect slopes
          mod_fixefs <- fixef(stat_model)
          
          return(c(
            mod_fixefs,
            0, if(is.null(unlist(stat_model@optinfo$conv$lme4))) { 1 } else { 0 },
            isSingular(stat_model),
            length(unique(idx)), length(idx)
          ))
        }
      },
      BOOTSTRAP_SIZE,
      parallel = "multicore",
      ncpus = N_CPUS
    )
    
    # Populate df
    glmm_boot_part <- tibble(
      Focal.taxon = SP_FOCAL,
      Neighbour.taxon = SP_NEI,
      Focal.stage = FOCAL_STAGE,
      Focal.vital.rate = FOCAL_VR,
      Intercept = as.numeric(boot_res$t[,1]),
      Length.cm.slope = as.numeric(boot_res$t[,2]) /
        sd(focal_data_sub$Length..cm., na.rm = TRUE),
      Overlap.small.slope = as.numeric(boot_res$t[,3]) /
        sd(focal_data_sub$Overlap.small, na.rm = TRUE),
      Overlap.large.slope = as.numeric(boot_res$t[,4]) /
        sd(focal_data_sub$Overlap.large, na.rm = TRUE),
      Overlap.other.slope = as.numeric(boot_res$t[,5]) /
        sd(focal_data_sub$Overlap.other, na.rm = TRUE),
      Error.occurred = as.numeric(boot_res$t[,6]),
      Has.converged = as.numeric(boot_res$t[,7]),
      Is.singular = as.numeric(boot_res$t[,8]),
      Sample.size.unique = as.numeric(boot_res$t[,9]),
      Sample.size = as.numeric(boot_res$t[,10]),
      Global.sample.size = as.numeric(nrow(focal_data_sub))
    )
    
    # Merge with main df
    glmm_boot_df <- bind_rows(glmm_boot_df, glmm_boot_part)
    
    # Save intermediate df
    write_csv(
      glmm_boot_df,
      here(
        "StageDependentComp",
        "result_data",
        "wp2",
        "bootstrap",
        sprintf(
          "%s_%s_bootstrap_interm.csv",
          SP_FOCAL_ABB,
          SP_NEI_ABB
        )
      )
    )
  }
}

# Save df to file
write_csv(
  glmm_boot_df,
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    sprintf(
      "%s_%s_bootstrap.csv",
      SP_FOCAL_ABB,
      SP_NEI_ABB
    )
  )
)

