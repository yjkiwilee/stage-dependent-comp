# Code test board

source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))

source(here("StageDependentComp", "code", "census_analysis", "func", "spatial_analysis.R"))

# ====== Extract & process focal data ======

# Subset census data to only include focal species
focal_data <- census_data %>%
  filter(
    !is.na(Length..cm.) &
    !is.na(X..cm.) &
    !is.na(Y..cm.) &
    !is.na(X..Capitulescences)
  )

# ====== Extract & process focal & other neighbour data =======

# Focal neighbour data, excluding dormant or dead individuals (i.e. Length = 0)
nei_data <- census_data %>%
  filter(
    !is.na(Length..cm.) &
    !is.na(X..cm.) &
    !is.na(Y..cm.) &
    !is.na(X..Capitulescences)
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
  )

overlap_data

# ===== MCMC using Metropolis–Hastings algorithm ======

# Initial values of parameters
ALPHA_0 <- 0

# Vectors to store parameter trajectory
alpha_traj <- c(ALPHA_0)

# Vectors to store log-likelihood
loglik_traj <- c(NA)

## ===== Evaluate initial likelihood =====

# Collapse overlap data
overlap_summ <- overlap_data %>%
  ungroup() %>%
  mutate( # Normalise neighbour length
    Length..cm..nei.norm = scale(Length..cm..nei)
  ) %>%
  group_by(Plot, Year, Tag) %>%
  summarise(
    Total.overlap.coeff.c = sum(
      (Prop.overlap *
        exp(Length..cm..nei.norm * ALPHA_0))[
          .data$Is.conspecific == TRUE
        ],
      na.rm = TRUE),
    Total.overlap.coeff.h = sum(
      (Prop.overlap *
         exp(Length..cm..nei.norm * ALPHA_0))[
           .data$Is.conspecific == FALSE
         ],
      na.rm = TRUE)
  ) %>%
  filter(
    abs(Total.overlap.coeff.c) != Inf &
      !is.na(Total.overlap.coeff.c) &
      abs(Total.overlap.coeff.h) != Inf &
      !is.na(Total.overlap.coeff.h)
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
# overlap_fit <- lmer(
#   X..Capitulescences ~ Total.overlap.coeff.norm + Length..cm.norm +
#     (1 | Plot) + (1 | Year) + (1 | Taxon),
#   data = focal_data_overlap,
#   control = lmerControl(
#     check.rankX = "stop.deficient"
#   )
# )
overlap_fit <- lmer(
  X..Capitulescences ~ Total.overlap.coeff.c.norm +
    Total.overlap.coeff.h.norm +
    Length..cm.norm +
    GS.Precipitation.inch +
    Length..cm.norm : Total.overlap.coeff.c.norm +
    GS.Precipitation.inch : Total.overlap.coeff.c.norm +
    Length..cm.norm : Total.overlap.coeff.h.norm +
    GS.Precipitation.inch : Total.overlap.coeff.h.norm +
    (1 | Plot) + (1 | Taxon),
  data = focal_data_overlap,
  control = lmerControl(
    check.rankX = "stop.deficient"
  )
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
ALPHA_SD <- 0.03

# Number of repetitions
MCMC_N_REP <- 50000

# Produce plot & save results every N iterations
PLOT_PERIOD <- 100

# Random seed
set.seed(1)

# Iterate for MCMC with Gaussian random walk
# https://www.r-bloggers.com/2020/06/a-slightly-more-advanced-mcmc-example/
for(mcmc_i in 1:MCMC_N_REP) {
  cat(paste0(mcmc_i, ": "))
  
  alpha_current <- alpha_traj[mcmc_i]
  
  alpha_new <- alpha_current + rnorm(1, sd = ALPHA_SD)
  
  cat(paste0(alpha_new, ": "))
  
  # Collapse overlap data
  overlap_summ <- overlap_data %>%
    ungroup() %>%
    mutate( # Normalise neighbour length
      Length..cm..nei.norm = scale(Length..cm..nei)
    ) %>%
    group_by(Plot, Year, Tag) %>%
    summarise(
      Total.overlap.coeff.c = sum(
        (Prop.overlap *
           exp(Length..cm..nei.norm * alpha_new))[
             .data$Is.conspecific == TRUE
           ],
        na.rm = TRUE),
      Total.overlap.coeff.h = sum(
        (Prop.overlap *
           exp(Length..cm..nei.norm * alpha_new))[
             .data$Is.conspecific == FALSE
           ],
        na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      abs(Total.overlap.coeff.c) != Inf &
        !is.na(Total.overlap.coeff.c) &
        abs(Total.overlap.coeff.h) != Inf &
        !is.na(Total.overlap.coeff.h)
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
  # overlap_fit <- lmer(
  #   X..Capitulescences ~ Total.overlap.coeff.norm + Length..cm.norm +
  #     (1 | Plot) + (1 | Year) + (1 | Taxon),
  #   data = focal_data_overlap,
  #   control = lmerControl(
  #     check.rankX = "stop.deficient"
  #   )
  # )
  # overlap_fit <- lmer(
  #   X..Capitulescences ~ Total.overlap.coeff.c.norm +
  #     Total.overlap.coeff.h.norm +
  #     Length..cm.norm +
  #     GS.Precipitation.inch +
  #     (1 | Plot),
  #   data = focal_data_overlap,
  #   control = lmerControl(
  #     check.rankX = "stop.deficient"
  #   )
  # )
  overlap_fit <- lmer(
    X..Capitulescences ~ Total.overlap.coeff.c.norm +
      Total.overlap.coeff.h.norm +
      Length..cm.norm +
      GS.Precipitation.inch +
      Length..cm.norm : Total.overlap.coeff.c.norm +
      GS.Precipitation.inch : Total.overlap.coeff.c.norm +
      Length..cm.norm : Total.overlap.coeff.h.norm +
      GS.Precipitation.inch : Total.overlap.coeff.h.norm +
      (1 | Plot) + (1 | Taxon),
    data = focal_data_overlap,
    control = lmerControl(
      check.rankX = "stop.deficient"
    )
  )
  
  # Evaluate new log likelihood and the acceptance ratio
  loglik_new <- logLik(overlap_fit)[1]
  log_accept_ratio <- loglik_new - loglik_traj[mcmc_i]
  
  # Accept or reject
  if(log(runif(1)) <= log_accept_ratio) {
    cat(paste0("accepted!\n"))
    # Accept
    alpha_traj <- c(alpha_traj, alpha_new)
    loglik_traj <- c(loglik_traj, loglik_new)
    accept_traj <- c(accept_traj, 1)
  } else {
    cat(paste0("rejected...\n"))
    # Reject
    alpha_traj <- c(alpha_traj, alpha_traj[mcmc_i])
    loglik_traj <- c(loglik_traj, loglik_traj[mcmc_i])
    accept_traj <- c(accept_traj, 0)
  }
  
  # Plot if needed
  if(mcmc_i %% PLOT_PERIOD == 0) {
    # Make DF for plotting
    plot_df <- tibble(
      iter = 1:(mcmc_i + 1),
      ll = loglik_traj,
      accept = accept_traj,
      alpha = alpha_traj
    )
    
    # Save DF
    write_csv(
      plot_df,
      here(
        "StageDependentComp",
        "result_data",
        "wp2",
        "mcmc",
        paste0(
          "single_alpha_mcmc.csv"
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
          "single_alpha_ll.png"
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
          "single_alpha_accept.png"
        )
      ),
      width = 6,
      height = 4,
      accept_plot
    )
    
    # Alpha parameter
    alpha_plot <- ggplot(plot_df, aes(x = iter, y = alpha)) +
      geom_line() +
      theme_jun1()
    
    ggsave(
      here(
        "StageDependentComp",
        "figures",
        "census_analysis",
        "mcmc",
        paste0(
          "single_alpha_alpha.png"
        )
      ),
      width = 6,
      height = 4,
      alpha_plot
    )
  }
}
