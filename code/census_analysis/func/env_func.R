################################################################################
# Functions for processing the environmental data
#
# Written by Young Jun Lee
# Nov 2024

# Data clean-up function
clean_env_data <- function(env_data, env_var_name, start_year, end_year) {
  env_data %>%
    select(!c(`10%`, `30%`, `70%`, `90%`,
              Min, Max, starts_with("Median"))) %>% # Remove summary statistics
    pivot_longer(!date, names_to = "year", values_to = env_var_name) %>% # Pivot longer
    filter(year <= end_year & year >= start_year) %>% # Filter for year
    mutate( # Split month & day
      year = as.numeric(year),
      month = as.numeric(str_split_fixed(date, "-", 2)[,1]),
      day = as.numeric(str_split_fixed(date, "-", 2)[,2])
    ) %>%
    mutate( # Calculate days since reference date & day in year
      days.since.ref = interval(ISOdate(start_year, 1, 1, hour = 0),
                                ISOdate(year, month, day, hour = 0)) %/% days(1),
      day.in.year = interval(ISOdate(year, 1, 1, hour = 0),
                             ISOdate(year, month, day, hour = 0)) %/% days(1)
    ) %>%
    select(!date) %>% # Remove date column
    relocate(all_of(env_var_name), .after = last_col()) %>% # Move temperature to last column
    relocate(days.since.ref) %>% # Move days since ref to start
    arrange(days.since.ref) %>% # Sort by days since ref
    filter(!is.na(days.since.ref)) # Remove days that don't exist (i.e. 29th Feb on non-leap years)
}

# Function for calculating annual environment summary
# NB: The input data must only be subsetted to contain a year
calc_ann_env_summ <- function(ann_env_data) {
  # Determine the growth season
  # Following Inouye et al. (2000) (https://doi.org/10.1073/pnas.97.4.1630)
  
  # Get year
  yr <- ann_env_data$year[[1]]
  
  # Determine the first day in the year without snow
  gr_start <- ann_env_data %>%
    filter(snow.depth.inch == 0) %>%
    filter(day.in.year == min(day.in.year))
  
  # Determine the first day of severe freeze in fall
  gr_end <- ann_env_data %>%
    filter(day.in.year > gr_start$day.in.year[[1]] & # After first day without snow
             temp.F < 25) %>% # Average temperature below 25 F
    filter(day.in.year == min(day.in.year))
  
  # Determine the difference in precipitation
  precip_gr_season <- gr_end$cumul.precip.inch[[1]] - 
    gr_start$cumul.precip.inch[[1]]
  
  # Additional variables
  # Average temperature over growth season
  avg_temp_gr_season <- mean((ann_env_data %>%
    filter(day.in.year >= gr_start$day.in.year[[1]] &
             day.in.year <= gr_end$day.in.year[[1]]))$temp.F, na.rm = TRUE)
  # Variance in temperature over growth season
  var_temp_gr_season <- var((ann_env_data %>%
    filter(day.in.year >= gr_start$day.in.year[[1]] &
             day.in.year <= gr_end$day.in.year[[1]]))$temp.F, na.rm = TRUE)
  
  # Build result dataframe
  res_df <- tibble(
    Year = yr,
    GS.Start.d.i.y = gr_start$day.in.year[[1]],
    GS.Start.month = gr_start$month[[1]],
    GS.Start.day = gr_start$day[[1]],
    GS.End.d.i.y = gr_end$day.in.year[[1]],
    GS.End.month = gr_end$month[[1]],
    GS.End.day = gr_end$day[[1]],
    GS.Window.day = gr_end$day.in.year[[1]] - gr_start$day.in.year[[1]],
    GS.Precipitation.inch = precip_gr_season,
    GS.Avg.temperature.F = avg_temp_gr_season,
    GS.Var.temperature.F = var_temp_gr_season
  )
  
  # Return result DF
  res_df
}



