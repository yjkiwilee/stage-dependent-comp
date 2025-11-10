################################################################################
# Script for processing the environmental data
#
# Written by Young Jun Lee
# Nov 2024

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {
  # Load packages
  pacman::p_load("tidyverse", "here", "stringr")
}

# Load relevant functions
source(here("StageDependentComp", "code", "census_analysis", "func", "env_func.R"))

# ===== Constants =====

ENV_START_YR <- 2013
ENV_END_YR <- 2023

# ===== Load environmental data & clean-up ======

## ===== Daily average temperature =====

# Load data
env_daily_temp <- read_csv(here("StageDependentComp", "data",
                                "schofield_pass_daily_avg_temp_f.csv"))

# Data clean-up
env_daily_temp <- clean_env_data(env_daily_temp, "temp.F", ENV_START_YR, ENV_END_YR)

# Replace temperature outliers with NA
env_daily_temp <- env_daily_temp %>%
  mutate(
    temp.F = ifelse(temp.F > 200 | temp.F < -50, NA, temp.F)
  )

# Plot monthly data
plt_env_temp <- ggplot(env_daily_temp, aes(x = as.factor(month), y = temp.F)) +
  geom_boxplot() +
  theme_jun1() +
  labs(
    x = "Month",
    y = "Temperature (F)",
    title = "Daily average temperature at Schofield Pass by month",
    caption = paste0("Data from ", ENV_START_YR, " to ", ENV_END_YR)
  )

plt_env_temp

# Save plot
ggsave(here("StageDependentComp", "figures", "census_analysis", "env_data",
            "daily_avg_temp_plot.png"),
       width = FIG_W, height = FIG_H)

## ===== Snow depth =====

# Load data
env_snow_depth <- read_csv(here("StageDependentComp", "data",
                                "schofield_pass_snow_depth_inch.csv"))

# Data clean-up
env_snow_depth <- clean_env_data(env_snow_depth, "snow.depth.inch", ENV_START_YR, ENV_END_YR)

# Plot monthly data
plt_env_snow <- ggplot(env_snow_depth, aes(x = as.factor(month), y = snow.depth.inch)) +
  geom_boxplot() +
  theme_jun1() +
  labs(
    x = "Month",
    y = "Snow depth (inches)",
    title = "Daily snow depth at Schofield Pass by month",
    caption = paste0("Data from ", ENV_START_YR, " to ", ENV_END_YR)
  )

plt_env_snow

# Save plot
ggsave(here("StageDependentComp", "figures", "census_analysis", "env_data",
            "snow_depth_plot.png"),
       width = FIG_W, height = FIG_H)

## ===== Cumulative precipitation =====

# Load data
env_cum_precip <- read_csv(here("StageDependentComp", "data",
                                "schofield_pass_precip_accum_inch.csv"))

# Data clean-up
env_cum_precip <- clean_env_data(env_cum_precip, "cumul.precip.inch", ENV_START_YR, ENV_END_YR)

# Empty DF to store processed cum precip data
proc_cum_precip <- NULL

# Force cumulative precipitation to start at 0 on 1st Jan
# NB: THIS SETS CUMULATIVE PRECIP ON 30TH SEP = CUMULATIVE PRECIP ON 1ST OCT
for(yr in unique(env_cum_precip$year)) {
  # Get subset of data for year
  yr_subset <- env_cum_precip %>% filter(year == yr)
  # Get cumulative precipitation on 1st Jan
  jan_cum_precip <- (env_cum_precip %>% filter(year == yr & day.in.year == 0))[[1, "cumul.precip.inch"]]
  # Get cumulative precipitation on 30th Sep
  sep_cum_precip <- (env_cum_precip %>% filter(year == yr & month == 9 & day == 30))[[1, "cumul.precip.inch"]]
  # Conditional mutate
  yr_subset <- yr_subset %>%
    mutate(
      cumul.precip.inch = ifelse(month <= 9, # If date is before 30th Sep
                                 cumul.precip.inch - jan_cum_precip, # Subtract 1st Jan precip
                                 cumul.precip.inch + sep_cum_precip - jan_cum_precip)
                                # Otherwise, add end-of-september precipitation before
                                # subtracting january precipitation
    )
  # Merge with result DF
  if(is.null(proc_cum_precip)) {
    proc_cum_precip <- yr_subset
  } else {
    proc_cum_precip <- bind_rows(proc_cum_precip, yr_subset)
  }
}

# Set new DF as env_cum_precip
env_cum_precip <- proc_cum_precip
# Remove temporary variable
rm(proc_cum_precip)

# Plot data
plt_env_cum_precip <- ggplot(env_cum_precip,
                             aes(x = day.in.year, y = cumul.precip.inch, group = year)) +
  geom_line(alpha = 0.7) +
  theme_jun1() +
  labs(
    x = "Day in year",
    y = "Cumulative precipitation (inches)",
    title = "Cumulative precipitation over year",
    caption = paste0("Data from ", ENV_START_YR, " to ", ENV_END_YR)
  )

plt_env_cum_precip

# Save plot
ggsave(here("StageDependentComp", "figures", "census_analysis", "env_data",
            "cum_precip_plot.png"),
       width = FIG_W, height = FIG_H)

# ===== Store cleaned data =====

# Merge environmental data into one
env_data <- env_daily_temp %>%
  full_join(env_snow_depth) %>%
  full_join(env_cum_precip)

env_data

# Write to CSV
write_csv(env_data,
          here("StageDependentComp", "result_data", "data_cleanup", "env_data_cleaned.csv"))

# ===== Calculate biologically relevant summary statistics for each year =====

# DF to store annual environmental summary
env_ann_summ <- NULL

# Iterate through years in the data
for(yr in unique(env_data$year)) {
  # Subset data
  yr_data <- env_data %>%
    filter(year == yr)
  
  # Get growth season summary
  ann_summ <- calc_ann_env_summ(yr_data)
  
  # Merge with DF
  if(is.null(env_ann_summ)) {
    env_ann_summ <- ann_summ
  } else {
    env_ann_summ <- bind_rows(env_ann_summ, ann_summ)
  }
}

# Store growth season summaries
write_csv(env_ann_summ,
          here("StageDependentComp", "result_data", "wp2", "annual_env_data.csv"))

