# Code test board

# Install here library
if(system.file(package="here") == "") {
  install.packages("here")
}
library("here")

# Load config parameters
source(here("StageDependentComp","code","census_analysis","census_analysis_config.R"))

# Load packages
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_packages.R"))

# Setup ggplot
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_setup_ggplot.R"))

source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))

source(here("StageDependentComp", "code", "census_analysis", "func", "spatial_analysis.R"))

FOCAL_SPP <- c(
  "Lupinus argenteus",
  "Ivesia gordonii",
  "Eriogonum umbellatum",
  "Heterotheca villosa"
)

N_SPP <- length(FOCAL_SPP)

# FOCAL_VR <- "Survival.future"
# CHAIN_ID <- 2
# SD_SD <- 0.01
# SD_SPP <- 0.0001

print(
  paste(
    FOCAL_VR,
    CHAIN_ID,
    SD_SD,
    SD_SPP
  )
)

# ====== Extract & process focal data ======

# Subset census data to only include focal species
focal_data <- census_data %>%
  filter(
    Taxon %in% FOCAL_SPP &
      !is.na(Length..cm.) &
      !is.na(X..cm.) &
      !is.na(Y..cm.) &
      !is.na(.data[[FOCAL_VR]]) &
      Length..cm. > 0
  ) %>%
  transform(
    Length..cm.norm = scale(Length..cm.),
    Precip.norm = scale(GS.Precipitation.inch)
  ) %>%
  mutate(
    Focal.VR = .data[[FOCAL_VR]]
  )

# ====== Extract & process focal & other neighbour data =======

# Focal neighbour data, excluding dormant or dead individuals (i.e. Length = 0)
nei_data <- focal_data

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
        filter( # Filter self pairings
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
  group_by(Taxon..nei) %>%
  transform( # Normalise neighbour length within each neighbour species
    Length..cm..nei.norm = scale(Length..cm..nei)
  ) %>%
  ungroup() 

# ===== MCMC using Metropolis–Hastings algorithm ======

## ===== Model parameters ======

mp <- list(
  # Stage-dependence coefficient for each species pair
  # NB: When neighbour size is being normalised, it's scaled
  # for each neighbour species.
  sd = matrix(
    rep(0, N_SPP^2),
    nrow = 1
  ),
  
  # Species pair-specific interaction coefficient
  spp = matrix(
    rep(0, N_SPP^2),
    nrow = 1
  )
)

# Mapping matrix
sp_map_mat <- matrix(
  1:(N_SPP^2),
  ncol = N_SPP
)
colnames(sp_map_mat) <- FOCAL_SPP
rownames(sp_map_mat) <- FOCAL_SPP

# Reverse mapping vectors
foc_sp_map <- rep(FOCAL_SPP, times = N_SPP)
nei_sp_map <- rep(FOCAL_SPP, each = N_SPP)

## ===== Evaluate initial likelihood =====

# Collapse overlap data
overlap_summ <- overlap_data %>%
  filter(
    !is.na(Prop.overlap) &
      !is.infinite(Prop.overlap) &
      !is.nan(Prop.overlap)
  ) %>%
  group_by(Plot, Year, Tag) %>%
  summarise(
    Overlap.coeff = sum(
      Prop.overlap *
        mp$spp[1, sp_map_mat[
          matrix(c(Taxon, Taxon..nei), ncol = 2)
        ]] *
        exp(
          Length..cm..nei.norm *
            mp$sd[1, sp_map_mat[
              matrix(c(Taxon, Taxon..nei), ncol = 2)
            ]]
        ),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

focal_data_overlap <- focal_data %>%
  left_join(
    overlap_summ,
    by = join_by(Plot, Year, Tag)
  ) 

focal_data_overlap

# Fit model to collapsed data
overlap_fit <- NULL

# model_form <- as.formula(paste(
#   FOCAL_VR,
#   " ~ Length..cm.norm + Precip.norm +",
#   "Length..cm.norm : Precip.norm +",
#   "(1 | Plot) + (1 | Taxon)"
# ))

model_form <- Focal.VR ~
  offset(Overlap.coeff) +
  Length..cm.norm + Precip.norm +
  Length..cm.norm:Precip.norm +
  (1 | Plot) + (1 | Taxon)

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

summary(overlap_fit)
logLik(overlap_fit)
AIC(overlap_fit)
# confint(overlap_fit)

# Vector to store log-likelihood
loglik_hist <- c(
  logLik(overlap_fit)[1]
)

# Vector to store acceptance history
accept_hist <- c(
  NA
)

## ===== Run MCMC iterations =====

# Number of repetitions
MCMC_N_REP <- 50000
# Plot the parameters every...
PLOT_PERIOD <- 500

# Random seed
set.seed(CHAIN_ID)

# Iterate for MCMC with Gaussian random walk
# https://www.r-bloggers.com/2020/06/a-slightly-more-advanced-mcmc-example/
for(mcmc_i in 1:MCMC_N_REP) {
  cat(paste0(mcmc_i, ": "))
  
  sd_curr <- mp$sd[mcmc_i,]
  spp_curr <- mp$spp[mcmc_i,]
  
  sd_new <- sd_curr + rnorm(N_SPP^2, sd = SD_SD)
  spp_new <- spp_curr + rnorm(N_SPP^2, sd = SD_SPP)
  
  # Collapse overlap data
  overlap_summ <- overlap_data %>%
    filter(
      !is.na(Prop.overlap) &
        !is.infinite(Prop.overlap) &
        !is.nan(Prop.overlap)
    ) %>%
    group_by(Plot, Year, Tag) %>%
    summarise(
      Overlap.coeff = sum(
        Prop.overlap *
          spp_new[sp_map_mat[
            matrix(c(Taxon, Taxon..nei), ncol = 2)
          ]] *
          exp(
            Length..cm..nei.norm *
              sd_new[sp_map_mat[
                matrix(c(Taxon, Taxon..nei), ncol = 2)
              ]]
          ),
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  focal_data_overlap <- focal_data %>%
    left_join(
      overlap_summ,
      by = join_by(Plot, Year, Tag)
    )
  
  # Fit model to collapsed data
  overlap_fit <- NULL
  
  # model_form <- as.formula(paste(
  #   FOCAL_VR,
  #   " ~ Length..cm.norm + Precip.norm +",
  #   "Length..cm.norm : Precip.norm +",
  #   "(1 | Plot) + (1 | Taxon)"
  # ))
  
  model_form <- Focal.VR ~
    offset(Overlap.coeff) +
    Length..cm.norm + Precip.norm +
    Length..cm.norm:Precip.norm +
    (1 | Plot) + (1 | Taxon)
  
  loglik_new <- NA
  log_accept_ratio <- -1
  
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
  
  if(!is.na(overlap_fit) & !is.na(loglik_hist[mcmc_i])) {
    # Evaluate new log likelihood and the acceptance ratio
    loglik_new <- logLik(overlap_fit)[1]
    log_accept_ratio <- loglik_new - loglik_hist[mcmc_i]
  }
  
  # Accept or reject
  if(log(runif(1)) <= log_accept_ratio) {
    cat(paste0("accepted!\n"))
    # Accept
    mp$sd <- rbind(mp$sd, sd_new)
    mp$spp <- rbind(mp$spp, spp_new)
    loglik_hist <- c(loglik_hist, loglik_new)
    accept_hist <- c(accept_hist, 1)
  } else {
    cat(paste0("rejected...\n"))
    # Reject
    mp$sd <- rbind(mp$sd, mp$sd[mcmc_i,])
    mp$spp <- rbind(mp$spp, mp$spp[mcmc_i,])
    loglik_hist <- c(loglik_hist, loglik_hist[mcmc_i])
    accept_hist <- c(accept_hist, 0)
  }
  
  # Plot if needed
  if(mcmc_i %% PLOT_PERIOD == 0) {
    # Make DF for plotting
    sd_df <- mp$sd
    colnames(sd_df) <- paste0("sd_", as.character(1:(N_SPP^2)))
    sd_df <- as_tibble(sd_df) 
    spp_df <- mp$spp
    colnames(spp_df) <- paste0("spp_", as.character(1:(N_SPP^2)))
    spp_df <- as_tibble(spp_df)
    plot_df <- tibble(
      iter = 1:(mcmc_i + 1),
      ll = loglik_hist,
      accept = accept_hist
    ) %>%
      bind_cols(sd_df, spp_df) %>%
      pivot_longer(
        starts_with("sd_") | starts_with("spp_"),
        names_to = c(".value", "sp_comb"),
        names_sep = "_"
      ) %>%
      mutate(sp_comb = as.numeric(sp_comb)) %>%
      mutate(
        foc_sp = foc_sp_map[sp_comb],
        nei_sp = nei_sp_map[sp_comb]
      ) %>%
      select(!sp_comb)
    
    # Save DF
    write_csv(
      plot_df,
      here(
        "StageDependentComp",
        "result_data",
        "wp2",
        "mcmc",
        paste0(
          FOCAL_VR,
          "_",
          CHAIN_ID,
          ".csv"
        )
      )
    )
    
    # SD
    sd_plot <- ggplot(plot_df, aes(x = iter, y = sd, color = foc_sp, linetype = nei_sp)) +
        geom_line() +
        theme_jun1()
    
    ggsave(
      here(
        "StageDependentComp",
        "figures",
        "census_analysis",
        "mcmc",
        paste0(
          FOCAL_VR,
          "_",
          CHAIN_ID,
          "_sd.png"
        )
      ),
      width = 6,
      height = 4,
      sd_plot
    )
    
    # SPP
    spp_plot <- ggplot(plot_df, aes(x = iter, y = spp, color = foc_sp, linetype = nei_sp)) +
        geom_line() +
        theme_jun1()
    
    ggsave(
      here(
        "StageDependentComp",
        "figures",
        "census_analysis",
        "mcmc",
        paste0(
          FOCAL_VR,
          "_",
          CHAIN_ID,
          "_spp.png"
        )
      ),
      width = 6,
      height = 4,
      spp_plot
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
          FOCAL_VR,
          "_",
          CHAIN_ID,
          "_ll.png"
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
          FOCAL_VR,
          "_",
          CHAIN_ID,
          "_accept.png"
        )
      ),
      width = 6,
      height = 4,
      accept_plot
    )
  }
}

