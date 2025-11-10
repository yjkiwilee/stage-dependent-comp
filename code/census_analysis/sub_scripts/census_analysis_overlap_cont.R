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
  Min.length..cm. = numeric(),
  Max.length..cm. = numeric(),
  Global.min.length..cm. = numeric(),
  Global.max.length..cm. = numeric(),
  Intercept = numeric(),
  Length.cm.slope = numeric(),
  Overlap.small.slope = numeric(),
  Overlap.small.length.interact = numeric(),
  Overlap.large.slope = numeric(),
  Overlap.large.length.interact = numeric(),
  Overlap.other.slope = numeric(),
  Overlap.other.length.interact = numeric(),
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
  focal_data_vr <- focal_data_stat
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
        Is.seedling == 0 &
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
    Overlap.small.norm + Length..cm..norm:Overlap.small.norm +
    Overlap.large.norm + Length..cm..norm:Overlap.large.norm +
    Overlap.other.norm + Length..cm..norm:Overlap.other.norm +
                     (1 | Year) + (1 | Plot)")
  # Convert to R equation
  model_eq <- as.formula(model_eq)
  
  # z-transform length & overlap variables
  norm_vars <- c("Length..cm.",
                 "Overlap.focal", "Overlap.small", "Overlap.large",
                 "Overlap.other", "Overlap")

  for(norm_var in norm_vars) {
    focal_data_vr[[paste0(norm_var, ".norm")]] <-
      (focal_data_vr[[norm_var]] - mean(focal_data_vr[[norm_var]])) /
      sd(focal_data_vr[[norm_var]])
  }

  # Bootstrap
  boot_res <- boot(
    focal_data_vr,
    function(df, idx) {
      sub_df <- df[idx,]

      did_throw_error <- 0

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
          min(sub_df[["Length..cm."]], na.rm = TRUE),
          max(sub_df[["Length..cm."]], na.rm = TRUE),
          0, 0, 0, 0, 0, 0, 0, 0,
          1, 0, 0, # 1 indicates error occurred
          length(unique(idx)), length(idx)
        ))
      } else {
        # Get fixed effect slopes
        mod_fixefs <- fixef(stat_model)

        return(c(
          min(sub_df[["Length..cm."]], na.rm = TRUE),
          max(sub_df[["Length..cm."]], na.rm = TRUE),
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
    Focal.vital.rate = FOCAL_VR,
    Min.length..cm. = as.numeric(boot_res$t[,1]),
    Max.length..cm. = as.numeric(boot_res$t[,2]),
    Global.min.length..cm. = as.numeric(min(focal_data_vr[["Length..cm."]], na.rm = TRUE)),
    Global.max.length..cm. = as.numeric(max(focal_data_vr[["Length..cm."]], na.rm = TRUE)),
    Intercept = as.numeric(boot_res$t[,3]),
    Length.cm.slope = as.numeric(boot_res$t[,4]) /
      sd(focal_data_vr$Length..cm., na.rm = TRUE),
    Overlap.small.slope = as.numeric(boot_res$t[,5]) /
      sd(focal_data_vr$Overlap.small, na.rm = TRUE),
    Overlap.small.length.interact = as.numeric(boot_res$t[,6]) / 
      (sd(focal_data_vr$Overlap.small, na.rm = TRUE) *
         sd(focal_data_vr$Length..cm., na.rm = TRUE)),
    Overlap.large.slope = as.numeric(boot_res$t[,7]) /
      sd(focal_data_vr$Overlap.large, na.rm = TRUE),
    Overlap.large.length.interact = as.numeric(boot_res$t[,8]) / 
      (sd(focal_data_vr$Overlap.large, na.rm = TRUE) *
         sd(focal_data_vr$Length..cm., na.rm = TRUE)),
    Overlap.other.slope = as.numeric(boot_res$t[,9]) /
      sd(focal_data_vr$Overlap.other, na.rm = TRUE),
    Overlap.other.length.interact = as.numeric(boot_res$t[,10]) / 
      (sd(focal_data_vr$Overlap.other, na.rm = TRUE) *
         sd(focal_data_vr$Length..cm., na.rm = TRUE)),
    Error.occurred = as.numeric(boot_res$t[,11]),
    Has.converged = as.numeric(boot_res$t[,12]),
    Is.singular = as.numeric(boot_res$t[,13]),
    Sample.size.unique = as.numeric(boot_res$t[,14]),
    Sample.size = as.numeric(boot_res$t[,15]),
    Global.sample.size = as.numeric(nrow(focal_data_vr))
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

glmm_boot_df <- read_csv(
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "lup_lup_bootstrap.csv"
  )
)

glmm_boot_plot <- glmm_boot_df %>%
  group_by(Focal.taxon, Neighbour.taxon, Focal.vital.rate) %>%
  mutate(
    Iteration = row_number()
  ) %>%
  group_by(Iteration, .add = TRUE) %>%
  mutate(
    Length..cm.start = min(focal_data_stat$Length..cm.),
    Length..cm.end = max(focal_data_stat$Length..cm.)
  ) %>%
  mutate(
    Overlap.small.slope.start = Overlap.small.slope +
      Overlap.small.length.interact * Length..cm.start,
    Overlap.small.slope.end = Overlap.small.slope +
      Overlap.small.length.interact * Length..cm.end
  ) %>%
  select(Focal.taxon, Neighbour.taxon, Focal.vital.rate, Iteration,
         Length..cm.start, Length..cm.end,
         Overlap.small.slope.start,
         Overlap.small.slope.end,
         Error.occurred) %>%
  pivot_longer(
    c(Overlap.small.slope.start, Overlap.small.slope.end),
    names_to = NULL,
    values_to = "Overlap.small.slope"
  ) %>%
  pivot_longer(
    c(Length..cm.start, Length..cm.end),
    names_to = NULL,
    values_to = "Length..cm."
  )

glmm_boot_plot

ggplot(
  glmm_boot_plot %>% filter(
    Focal.vital.rate == "X..Capitulescences" &
      Error.occurred != 1
  ),
  aes(
    x = Length..cm.,
    y = Overlap.small.slope,
    group = interaction(
      Focal.taxon,
      Neighbour.taxon,
      Focal.vital.rate,
      Iteration
    )
  )
) +
  geom_line()


