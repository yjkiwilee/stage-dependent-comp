# Code test board

source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))

source(here("StageDependentComp", "code", "census_analysis", "func", "spatial_analysis.R"))

# ====== Extract & process focal data ======

FOCAL_SP <- "Eriogonum umbellatum"
FOCAL_VR <- "Survival.future"

# Subset census data to only include focal species
focal_data <- census_data %>%
  mutate(
    Focal.VR = .data[[FOCAL_VR]]
  ) %>%
  filter(
    Taxon == FOCAL_SP &
    !is.na(Length..cm.) &
    Length..cm. > 0 &
    !is.na(X..cm.) &
    !is.na(Y..cm.) &
    !is.na(Focal.VR)
  ) %>%
  transform(
    Length..cm.norm = scale(Length..cm.),
    Precip.norm = scale(GS.Precipitation.inch)
  )

# ====== Extract & process focal & other neighbour data =======

# Focal neighbour data, excluding dormant or dead individuals (i.e. Length = 0)
nei_data <- census_data %>%
  filter(
    !is.na(Length..cm.) &
    Length..cm. > 0 &
    !is.na(X..cm.) &
    !is.na(Y..cm.)
  )

focal_data %>%
  arrange(as.numeric(Year), as.numeric(Plot))

f_data_grouped <- focal_data %>%
  group_by(Plot, Year)

overlap_data <- f_data_grouped %>%
  group_modify(
    function(f_sub_df, key_df) {
      cat(paste0(
        key_df$Plot[1], " ",
        key_df$Year[1], "\n"
      ))
      
      n_sub_df <- nei_data %>%
        filter(
          Plot == key_df$Plot[1] &
            Year == key_df$Year[1]
        ) %>%
        select(
          X..cm., Y..cm., Length..cm., Taxon, Tag
        ) %>%
        rename(
          X..cm..nei = X..cm.,
          Y..cm..nei = Y..cm.,
          Length..cm..nei = Length..cm.,
          Taxon..nei = Taxon,
          Tag..nei = Tag
        )
      
      pair_data <- tidyr::crossing(
        f_sub_df,
        n_sub_df
      ) %>%
        filter( # Filter self
          Tag != Tag..nei
        )
      
      pair_data <- pair_data %>%
        mutate(
          Distance..cm. = sqrt(
            (X..cm. - X..cm..nei)^2 + (Y..cm. - Y..cm..nei)^2
          )
        )
      
      if(nrow(pair_data) > 0) {
        pair_data <- pair_data %>%
          mutate( # Calculate 'density'; note that focal and neighbour are swapped
            Prop.overlap = prop_overlap(
              X..cm..nei, Y..cm..nei, Length..cm..nei,
              X..cm., Y..cm., Length..cm.,
              0, 200, 0, 200
            ) /
              area_circle_in_bound(X..cm., Y..cm., Length..cm.,
                                   0, 200, 0, 200) * # Divide by area
              100^2 # Convert to density per square metre
          )
      } else {
        pair_data$Prop.overlap = numeric()
      }
      
      return(pair_data)
    }
  ) %>%
  mutate(
    Is.conspecific = Taxon == Taxon..nei
  ) %>%
  transform(
    Distance..cm.norm = scale(Distance..cm.)
  ) %>%
  ungroup() %>%
  mutate( # Calculate relative neighbour length
    Length..cm..nei.rel = Length..cm..nei / Length..cm.
  )

overlap_data

# ===== MCMC using Metropolis–Hastings algorithm ======

# Initial values of parameters
PARAMS_0 <- c(
  sd = 0,
  l = 0,
  sd_l = 0
)

# Matrix to store parameter trajectory
params_traj <- matrix(PARAMS_0, nrow = 1)
colnames(params_traj) <- names(PARAMS_0)

# Vectors to store log-likelihood
loglik_traj <- c(NA)

## ===== Evaluate initial likelihood =====

# Collapse overlap data
overlap_summ <- overlap_data %>%
  group_by(Plot, Year, Tag) %>%
  summarise(
    Total.overlap.coeff.c = sum(
      (
        Prop.overlap *
          exp(
            Length..cm. * PARAMS_0["l"] +
            Length..cm..nei * PARAMS_0["sd"] +
            Length..cm. * Length..cm..nei * PARAMS_0["sd_l"]
          )
      )[
        .data$Is.conspecific == TRUE
      ],
      na.rm = TRUE),
    Total.overlap.coeff.h = sum(
      (
        Prop.overlap *
          exp(
            Length..cm. * PARAMS_0["l"] +
              Length..cm..nei * PARAMS_0["sd"] +
              Length..cm. * Length..cm..nei * PARAMS_0["sd_l"]
          )
      )[
        .data$Is.conspecific == FALSE
      ],
      na.rm = TRUE),
    .groups = "drop"
  )

focal_data_overlap <- focal_data %>%
  left_join(
    overlap_summ,
    by = join_by(Plot, Year, Tag)
  ) %>%
  transform(
    Total.overlap.coeff.c.norm = scale(Total.overlap.coeff.c),
    Total.overlap.coeff.h.norm = scale(Total.overlap.coeff.h),
    Length..cm.norm = scale(Length..cm.)
  )

# Fit model to collapsed data
model_form <- Focal.VR ~
  Total.overlap.coeff.c.norm + Total.overlap.coeff.c.norm:Precip.norm +
  Total.overlap.coeff.h.norm + Total.overlap.coeff.h.norm:Precip.norm +
  Length..cm.norm + Length..cm.norm:Precip.norm +
  Precip.norm +
  (1 | Plot)

overlap_fit <- NA

overlap_fit <- tryCatch(
  {
    if(FOCAL_VR %in% c("X..Capitulescences", "Growth.future")) {
      overlap_fit <- lmer(
        model_form,
        data = focal_data_overlap,
        control = lmerControl(
          check.rankX = "stop.deficient"
        )
      )
    } else {
      overlap_fit <- glmer(
        model_form,
        data = focal_data_overlap,
        family = "binomial",
        control = glmerControl(
          optimizer = "bobyqa",
          check.rankX = "stop.deficient"
        )
      )
    }
  },
  error = function(e) {
    # Error if model fails to converge
    cat("error! ")
    return(NA)
  }
)

summary(overlap_fit)
logLik(overlap_fit)
AIC(overlap_fit)
# confint(overlap_fit)

# Vector to store log-likelihood
loglik_traj <- c(
  logLik(overlap_fit)[1]
)

# Vector to store acceptance
accept_traj <- c(
  NA
)

## ===== Run MCMC iterations =====

# Standard deviation for the Gaussian proposal function
# PARAMS_SD <- c(
#   sd = 0.02,
#   l = 0.04,
#   sd_l = 0.004,
#   kernel = 0.000001
# )
# Larg
# PARAMS_SD <- c(
#   sd = 0.02,
#   l = 0.05,
#   sd_l = 0.002,
#   kernel = 0.00001
# )
# Eumb Fertility
# PARAMS_SD <- c(
#   sd = 0.03,
#   l = 0.03,
#   sd_l = 0.00005
# )
PARAMS_SD <- c(
  sd = 0.5,
  l = 0.5,
  sd_l = 0.05
)

# Number of repetitions
MCMC_N_REP <- 50000

# Produce plot & save results every N iterations
PLOT_PERIOD <- 20

# Random seed
set.seed(1)

# Iterate for MCMC with Gaussian random walk
# https://www.r-bloggers.com/2020/06/a-slightly-more-advanced-mcmc-example/
for(mcmc_i in 1:MCMC_N_REP) {
  cat(paste0(mcmc_i, ": "))
  
  params_new <- params_traj[mcmc_i,] + rnorm(ncol(params_traj), sd = PARAMS_SD)
  
  # Collapse overlap data
  overlap_summ <- overlap_data %>%
    group_by(Plot, Year, Tag) %>%
    summarise(
      Total.overlap.coeff.c = sum(
        (
          Prop.overlap *
            exp(
              Length..cm. * params_new["l"] +
                Length..cm..nei * params_new["sd"] +
                Length..cm. * Length..cm..nei * params_new["sd_l"]
            )
        )[
          .data$Is.conspecific == TRUE
        ],
        na.rm = TRUE),
      Total.overlap.coeff.h = sum(
        (
          Prop.overlap *
            exp(
              Length..cm. * params_new["l"] +
                Length..cm..nei * params_new["sd"] +
                Length..cm. * Length..cm..nei * params_new["sd_l"]
            )
        )[
          .data$Is.conspecific == FALSE
        ],
        na.rm = TRUE),
      .groups = "drop"
    )
  
  focal_data_overlap <- focal_data %>%
    left_join(
      overlap_summ,
      by = join_by(Plot, Year, Tag)
    ) %>%
    transform(
      Total.overlap.coeff.c.norm = scale(Total.overlap.coeff.c),
      Total.overlap.coeff.h.norm = scale(Total.overlap.coeff.h),
      Length..cm.norm = scale(Length..cm.)
    )
  
  # Fit model to collapsed data
  loglik_new <- NA
  log_accept_ratio <- -Inf
  
  overlap_fit <- NULL
  
  overlap_fit <- tryCatch(
    {
      if(FOCAL_VR %in% c("X..Capitulescences", "Growth.future")) {
        overlap_fit <- lmer(
          model_form,
          data = focal_data_overlap,
          control = lmerControl(
            check.rankX = "stop.deficient"
          )
        )
      } else {
        overlap_fit <- glmer(
          model_form,
          data = focal_data_overlap,
          family = "binomial",
          control = glmerControl(
            optimizer = "bobyqa",
            check.rankX = "stop.deficient"
          )
        )
      }
    },
    error = function(e) {
      # Error if model fails to converge
      cat("error! ")
      return(NA)
    }
  )
  
  # Evaluate new log likelihood and the acceptance ratio
  if(!is.null(overlap_fit) & !is.na(loglik_traj[mcmc_i])) {
    # Evaluate new log likelihood and the acceptance ratio
    loglik_new <- logLik(overlap_fit)[1]
    log_accept_ratio <- loglik_new - loglik_traj[mcmc_i]
  }
  
  # Accept or reject
  if(log(runif(1)) <= log_accept_ratio) {
    cat(paste0("accepted!\n"))
    # Accept
    params_traj <- rbind(params_traj, params_new)
    loglik_traj <- c(loglik_traj, loglik_new)
    accept_traj <- c(accept_traj, 1)
  } else {
    cat(paste0("rejected...\n"))
    # Reject
    params_traj <- rbind(params_traj, params_traj[mcmc_i,])
    loglik_traj <- c(loglik_traj, loglik_traj[mcmc_i])
    accept_traj <- c(accept_traj, 0)
  }
  
  # Plot if needed
  if(mcmc_i %% PLOT_PERIOD == 0) {
    # Make DF for plotting
    plot_df <- tibble(
      iter = 1:(mcmc_i + 1),
      ll = loglik_traj,
      accept = accept_traj
    ) %>%
      bind_cols(params_traj)
    
    # Save DF
    write_csv(
      plot_df,
      here(
        "StageDependentComp",
        "result_data",
        "wp2",
        "mcmc",
        paste0(
          "ipm_mcmc.csv"
        )
      )
    )
    
    # Log likelihood
    ll_plot <- ggplot(plot_df, aes(x = iter, y = ll)) +
      geom_line() +
      theme_jun1()
    
    ggsave(
      here(
        "StageDependentComp",
        "figures",
        "census_analysis",
        "mcmc",
        paste0(
          "ipm_ll.png"
        )
      ),
      width = 6,
      height = 4,
      ll_plot
    )
    
    # Acceptance rate
    accept_plot <- ggplot(plot_df, aes(x = iter, y = accept)) +
      geom_point() +
      geom_smooth() +
      theme_jun1()
    
    ggsave(
      here(
        "StageDependentComp",
        "figures",
        "census_analysis",
        "mcmc",
        paste0(
          "ipm_accept.png"
        )
      ),
      width = 6,
      height = 4,
      accept_plot
    )
    
    # Parameters
    plot_df_long <- plot_df %>%
      pivot_longer(
        4:ncol(plot_df),
        names_to = "param",
        values_to = "param_val"
      )
    
    params_plot <- ggplot(plot_df_long, aes(x = iter, y = param_val, color = param)) +
      geom_line() +
      facet_wrap(~ param, ncol = 2, scales = "free_y") +
      theme_jun1()
    
    ggsave(
      here(
        "StageDependentComp",
        "figures",
        "census_analysis",
        "mcmc",
        paste0(
          "ipm_params.png"
        )
      ),
      width = 6,
      height = 9,
      params_plot
    )
  }
}
