# Functions needed to delineate stage boundary & calculate vital rates based on these stage boundaries

# ====== Function for finding the threshold length given a desired proportion of reproductives in the juvenile stage =====
calc_thres_len <- function(target_prop, sp_data, len_step = 0.1) {
  # Get max length in the dataset
  max_len <- max(sp_data$Length..cm.)
  
  # Stage boundary
  thres_len <- 0
  
  # Iterate until the proportion of reproductives is greater than the target prop
  while(thres_len <= max_len) {
    # Temporary threshold length
    temp_thres_len <- thres_len + len_step
    
    # Filter for plants whose length is smaller than or equal to the threshold
    temp_sub <- sp_data %>%
      filter(Length..cm. <= temp_thres_len)
    
    # Calculate proportion of reproductives
    prop_reprod <- nrow(filter(temp_sub, X..Capitulescences > 0)) / nrow(temp_sub)
    
    # If proportion of reproductives is greater than the threshold, exit loop
    if(prop_reprod > target_prop) { break }
    
    # Otherwise increase threshold length
    thres_len <- temp_thres_len
  }
  
  # Return the calculated stage boundary
  thres_len
}

# ===== Function for calculating progression given a subset of the data & the juvenile/adult boundary =====
calc_stage_prog <- function(sub_data, thres_len) {
  # Calculate progression
  growth_vec <- ifelse(is.na(sub_data$Length..cm..next) | sub_data$Survival.future == 0, NA, sub_data$Length..cm. <= thres_len & sub_data$Length..cm..next > thres_len)
  # Convert TRUE & FALSE to 1 and 0
  growth_vec <- 1 * growth_vec
  # Return progression vector
  growth_vec
}

# ===== Function for calculating retrogression given a subset of the data & the juvenile/adult boundary =====
calc_stage_retrog <- function(sub_data, thres_len) {
  # Calculate retrogression
  growth_vec <- ifelse(is.na(sub_data$Length..cm..next) | sub_data$Survival.future == 0, NA, sub_data$Length..cm. > thres_len & sub_data$Length..cm..next <= thres_len)
  # Convert TRUE & FALSE to 1 and 0
  growth_vec <- 1 * growth_vec
  # Return growth vector
  growth_vec
}

# ===== Function for calculating the recruitment & recruitment factor for each year ======
# NB: THIS FUNCTION ASSUMES THAT THE DATA HAS BEEN SUBSETTED SO IT ONLY CONTAINS ONE SPECIES
calc_rec_fact <- function(sp_data) {
  # Tibble to store the recruitment factor
  rec_fact_df <- tibble(
    Year = numeric(),
    Total.capitulescences = numeric(),
    N.new.tags.next.year = numeric(),
    Recruitment.factor = numeric()
  )
  
  # Get years
  yrs <- sort(unique(sp_data$Year))
  
  # Iterate through year and count new unique individuals
  for(i in 1:(length(yrs) - 1)) { # Finish loop before the end since we don't have data for recruitment in 2023
    # Get current year
    yr <- yrs[i]
    # Count the total capitulescences across the population for that year
    cap_tot <- sum((sp_data %>% filter(Year == yr))$X..Capitulescences)
    # Count recruits
    new_indivs <- nrow(sp_data %>% filter(Year == (yr + 1) & Is.recruit == 1))
    # Calculate recruitment factor
    rec_fact <- new_indivs / cap_tot
    # Create new row in df
    rec_fact_df <- rbind(rec_fact_df, tibble(
      Year = yr,
      Total.capitulescences = cap_tot,
      N.new.tags.next.year = new_indivs,
      Recruitment.factor = rec_fact
    ))
  }
  
  # Return the df
  rec_fact_df
}

# ====== Function for obtaining a dataframe summarising the overall annual vital rates of a given species =======
summ_ann_vrs <- function(vr_data) {
  # Tibble for storing the annual vital rates
  vr_df <- tibble(
    Year = numeric(),
    Small.survival = numeric(),
    Large.survival = numeric(),
    Progression = numeric(),
    Retrogression = numeric(),
    Fertility = numeric(),
    N.small = numeric(),
    N.large = numeric()
  )
  
  # Get the years represented in the dataset
  yrs <- sort(unique(vr_data$Year))
  # Remove last element from the represented years
  yrs <- yrs[1:(length(yrs)-1)]
  
  # Iterate through the years to calculate annual vital rates
  for(yr in yrs) {
    # Subset data based on year
    vr_data_sub <- vr_data %>%
      filter(Year == yr)
    rec_fact_sub <- rec_fact %>%
      filter(Year == yr)
    
    # Small & large survival
    ss <- (vr_data_sub %>%
      filter(Stage == "S" & Died.this.census.final == 0 & !is.na(Survival.future)) %>%
      summarise(
        prop_surv = sum(Survival.future) / n()
      ))[[1, 1]]
    sl <- (vr_data_sub %>%
      filter(Stage == "L" & Died.this.census.final == 0 & !is.na(Survival.future)) %>%
      summarise(
        prop_surv = sum(Survival.future) / n()
      ))[[1, 1]]
    
    # Progression
    g <- (vr_data_sub %>%
      filter(Stage == "S" & !is.na(Progression)) %>%
      summarise(
        prop_prog = sum(Progression) / n()
      ))[[1, 1]]
    # Retrogression
    r <- (vr_data_sub %>%
            filter(Stage == "L" & !is.na(Retrogression)) %>%
            summarise(
              prop_prog = sum(Retrogression) / n()
            ))[[1, 1]]
    
    # Fertility
    rec_next_year <- nrow(vr_data %>%
      filter(Year == (yr + 1)) %>%
      filter(Is.recruit))
    f <- rec_next_year /
      nrow(filter(vr_data_sub, Stage == "L" & Died.this.census.final == 0))
    # Overall recruitment per large individual
    ro <- rec_fact_sub$N.new.tags.next.year[1] /
      nrow(filter(vr_data_sub, Stage == "L" & Died.this.census.final == 0))
    
    # Number of living small individuals
    ns <- nrow(filter(vr_data_sub, Stage == "S" & Died.this.census.final == 0))
    # Numberr of living large individuals
    nl <- nrow(filter(vr_data_sub, Stage == "L" & Died.this.census.final == 0))
    
    
    # Append vital rates to the dataframe
    vr_df <- rbind(vr_df, tibble(
      Year = yr,
      Recruitment.factor = rf,
      Small.survival = ss,
      Large.survival = sl,
      Progression = g,
      Retrogression = r,
      Reproduction = f,
      Overall.recruitment = ro,
      N.small = ns,
      N.large = nl
    ))
  }
  
  # Return the df of vital rates
  vr_df
}


