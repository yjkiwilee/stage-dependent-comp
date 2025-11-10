# Functions for manipulating the dataset

# Function for tracking individuals to add future vital rates & age
track_indivs <- function(sub_data) {
  # Insert future survival column & ID column into data
  sub_data_temp <- sub_data %>%
    mutate(
      Survival.future = rep(NA, nrow(.)),
      Length..cm..next = rep(NA, nrow(.)),
      Growth.future = rep(NA, nrow(.)),
      Age = rep(NA, nrow(.)),
      ID = consecutive_id(Plot, Year, Tag)
    )
  
  # Get vector of individual tags
  tags <- unique(sub_data$Tag)
  
  # Iterate through individual tag
  for(tag in tags) {
    # Get all records of individual
    indiv_data <- sub_data_temp %>%
      filter(Tag == tag) %>%
      arrange(Year) # Sort by year
    
    # Make vector of future survival
    future_surv <- indiv_data$Survival
    future_surv <- if(length(future_surv) == 1) { c(NA) } else { c(future_surv[2:length(future_surv)], NA) }
    
    # Make vector of future length
    future_length <- indiv_data$Length..cm.
    future_length <- if(length(future_length) == 1) { c(NA) } else { c(future_length[2:length(future_length)], NA) }
    
    # Make vector of age ONLY IF seedling has been observed
    age <- if(indiv_data$Is.seedling[1] == 1) { indiv_data$Year - min(indiv_data$Year) }
      else { rep(NA, nrow(indiv_data)) }
    
    # Get ID, iterate through it inserting the future states & age
    ids <- indiv_data$ID
    for(i in 1:nrow(indiv_data)) {
      sub_data_temp[ids[i], "Survival.future"] <- future_surv[i]
      sub_data_temp[ids[i], "Length..cm..next"] <- future_length[i]
      sub_data_temp[ids[i], "Age"] <- age[i]
    }
  }
  
  # Calculate future growth
  sub_data_temp <- sub_data_temp %>%
    mutate(
      Growth.future = Length..cm..next - Length..cm.
    )
  
  # Return sub_data_temp
  sub_data_temp %>%
    select(!ID)
}


