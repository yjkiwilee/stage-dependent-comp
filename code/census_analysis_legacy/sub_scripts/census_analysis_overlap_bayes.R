# Code test board

source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))

source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_set_sp.R"))

source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_stat_packages.R"))

source(here("StageDependentComp", "code", "census_analysis", "func", "spatial_analysis.R"))

# Create census dataset with unique ID
census_data_id <- census_data %>%
  mutate(
    ID = paste(Plot, Year, Tag, sep = "/")
  )

focal_data <- census_data_id %>%
  filter(
    Taxon == SP1 &
      Died.this.census.final == 0 &
      !is.na(Growth.future) &
      Length..cm..next != 0
  )

nei_data <- census_data_id %>%
  filter(
    Taxon == SP1
  )

other_data <- census_data_id %>%
  filter(
    !(ID %in% nei_data$ID)
  )


# Get pairwise overlap & total overlap effects

focal_pairs <- calc_overlap_pairs(focal_data, nei_data)
other_pairs <- calc_overlap_pairs(focal_data, other_data)

focal_data$Prop.overlap.focal <- calc_overlap_tot_frompair(focal_pairs, focal_data)
focal_data$Prop.overlap.other <- calc_overlap_tot_frompair(other_pairs, focal_data)

focal_om <- calc_overlap_mat(focal_pairs, focal_data$ID, nei_data$ID)
other_om <- calc_overlap_mat(other_pairs, focal_data$ID, other_data$ID)

# ====== Bayesian test ======

focal_data_bayes <- focal_data

# Code-ify year and plot

# Generate mapping for year & plot
yr_map <- tibble(Year = unique(focal_data_bayes$Year))
yr_map$Id <- seq(1, nrow(yr_map))
plt_map <- tibble(Plot = unique(focal_data_bayes$Plot))
plt_map$Id <- seq(1, nrow(plt_map))

# Replace years & plots with incremental value
focal_data_bayes$Year.id <- as.factor(yr_map$Id[match(focal_data_bayes$Year, yr_map$Year)])
focal_data_bayes$Plot.id <- as.factor(plt_map$Id[match(focal_data_bayes$Plot, plt_map$Plot)])

# Model to extract effect of year & plot
yp_model <- lm(Growth.future ~ Year.id + Plot.id, data = focal_data_bayes)

summary(yp_model)

year_effs <- yp_model$coefficients[
  grepl("^Year.id", names(yp_model$coefficients))
]
names(year_effs) <- sub("^Year.id([0-9]+)", "\\1", names(year_effs))
year_effs["1"] <- 0
year_effs <- year_effs[match(as.character(seq(1, length(year_effs))), names(year_effs))]

plot_effs <- yp_model$coefficients[
  grepl("^Plot.id", names(yp_model$coefficients))
]
names(plot_effs) <- sub("^Plot.id([0-9]+)", "\\1", names(plot_effs))
plot_effs["1"] <- 0
plot_effs <- plot_effs[match(as.character(seq(1, length(plot_effs))), names(plot_effs))]

rand_intercept <- yp_model$coefficients["(Intercept)"]

# Gather data
bayesian_data <- list(
  n = nrow(focal_data_bayes),
  n_n = nrow(nei_data),
  n_plots = nrow(plt_map),
  n_years = nrow(yr_map),
  
  l_f = focal_data_bayes$Length..cm.,
  l_n = nei_data$Length..cm.,
  
  overlap = focal_om,
  overlap_other = focal_data_bayes$Prop.overlap.other,
  
  vr = focal_data_bayes$Growth.future,
  
  plot = as.numeric(focal_data_bayes$Plot.id),
  year = as.numeric(focal_data_bayes$Year.id),
  
  rand_intercept = rand_intercept,
  plot_effs = plot_effs,
  year_effs = year_effs,
  
  vr_sd = sd(focal_data_bayes$Growth.future)
)

# ===== Fit Bayesian model ======

b_fit <- stan(
  file = here("StageDependentComp", "code", "census_analysis", "sub_scripts",
              "census_analysis_bmodel_g_overlap.stan"),
  data = bayesian_data,
  chains = 2,
  warmup = 50,
  iter = 60,
  core = 2
)

b_fit

plot(b_fit)# Gather data
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
