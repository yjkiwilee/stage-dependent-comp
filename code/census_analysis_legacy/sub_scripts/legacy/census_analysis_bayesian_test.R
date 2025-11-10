################################################################################
#' Script for testing the use of a Bayesian model to determine stage-dependence
#' Written by Young Jun Lee
#' Feb 2024


# Load packages
pacman::p_load("here")
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
            "census_analysis_load_packages.R"))
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
            "census_analysis_load_stat_packages.R"))

# Load spatial analysis functions
source(here("StageDependentComp", "code", "census_analysis", "func",
            "spatial_analysis.R"))

# Load processed census data
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
            "census_analysis_load_proc_data.R"))

# Get focal species
source(here("StageDependentComp", "code", "census_analysis", "sub_scripts",
            "census_analysis_set_sp.R"))

# ===== Subset data ======

BOUNDARY_WIDTH <- 50
PLOT_DIM <- 200

dat_f <- sp1_data %>%
  filter(
    !is.na(Growth.future) &
      X..cm. >= BOUNDARY_WIDTH & X..cm. <= PLOT_DIM - BOUNDARY_WIDTH &
      Y..cm. >= BOUNDARY_WIDTH & Y..cm. <= PLOT_DIM - BOUNDARY_WIDTH
  )
dat_c <- sp1_data
dat_h <- sp2_data
dat_o <- census_data %>% filter(!(Taxon %in% c(SP1, SP2)))

# Generate mapping for year & plot
yr_map <- tibble(Year = unique(dat_f$Year))
yr_map$Id <- seq(1, nrow(yr_map))
plt_map <- tibble(Plot = unique(dat_f$Plot))
plt_map$Id <- seq(1, nrow(plt_map))

# Replace years & plots with incremental value
dat_f$Year.id <- as.factor(yr_map$Id[match(dat_f$Year, yr_map$Year)])
dat_f$Plot.id <- as.factor(plt_map$Id[match(dat_f$Plot, plt_map$Plot)])

# ===== Run ANOVA with year & plot to obtain random effect priors =====

dat_rand_m <- lm(Growth.future ~ Year.id + Plot.id, data = dat_f)

year_effs <- dat_rand_m$coefficients[
  grepl("^Year.id", names(dat_rand_m$coefficients))
  ]
names(year_effs) <- sub("^Year.id([0-9]+)", "\\1", names(year_effs))
year_effs["1"] <- 0
year_effs <- year_effs[match(as.character(seq(1, length(year_effs))), names(year_effs))]

plot_effs <- dat_rand_m$coefficients[
  grepl("^Plot.id", names(dat_rand_m$coefficients))
]
names(plot_effs) <- sub("^Plot.id([0-9]+)", "\\1", names(plot_effs))
plot_effs["1"] <- 0
plot_effs <- plot_effs[match(as.character(seq(1, length(plot_effs))), names(plot_effs))]

rand_intercept <- dat_rand_m$coefficients["(Intercept)"]

# ===== Construct distance matrices =====

dist_c <- calc_dist_mat_across(dat_f, dat_c)
dist_h <- calc_dist_mat_across(dat_f, dat_h)
dist_o <- calc_dist_mat_across(dat_f, dat_o)

# Replace Inf with large distances
dist_c[dist_c == Inf] <- 1e4
dist_h[dist_h == Inf] <- 1e4
dist_o[dist_o == Inf] <- 1e4

# ===== Collect data into list =====

# Gather data
bayesian_data <- list(
  n = nrow(dat_f),
  n_c = nrow(dat_c),
  n_h = nrow(dat_h),
  n_o = nrow(dat_o),
  n_tot = nrow(dat_c) + nrow(dat_h) + nrow(dat_o),
  n_plots = nrow(plt_map),
  n_years = nrow(yr_map),
  
  l_f = dat_f$Length..cm.,
  l_all = c(dat_c$Length..cm., dat_h$Length..cm., dat_o$Length..cm.),
  
  dist = cbind(dist_c, dist_h, dist_o),
  
  vr = dat_f$Growth.future,
  
  plot = as.numeric(dat_f$Plot.id),
  year = as.numeric(dat_f$Year.id),
  
  exp_mean_d_t = 30,
  
  rand_intercept = rand_intercept,
  plot_effs = plot_effs,
  year_effs = year_effs,
  
  vr_sd = sd(dat_f$Growth.future)
)

# ===== Fit Bayesian model ======

b_fit <- stan(
  file = here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_bmodel_g_2.stan"),
  data = bayesian_data,
  # chains = 4,
  warmup = 80,
  iter = 100,
  core = 4,
  control = list(
    adapt_delta = 0.99
  )
)

b_fit

plot(b_fit)




