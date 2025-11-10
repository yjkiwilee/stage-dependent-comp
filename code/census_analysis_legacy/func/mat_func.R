# Functions needed to delineate stage boundary & calculate vital rates based on these stage boundaries
# Young Jun Lee, Nov 2024

# ===== Function for calculating a single matrix based on the vital rates =====
calc_single_mat <- function(rf, ss, sl, g, r, f) {
  # Calculate matrix elements
  # Stasis in small stage
  s_stasis <- ss * (1 - g)
  # Stasis in large stage
  l_stasis <- sl * (1 - r)
  # Transition from small to large
  s_to_l <- ss * g
  # Transition from large to small (reproduction + retrogression)
  l_to_s <- f * rf + sl * r
  
  # Build matrix
  ann_mat <- matrix(c(
    s_stasis, l_to_s,
    s_to_l, l_stasis
  ), ncol = 2, byrow = TRUE)
  
  # Label matrix
  rownames(ann_mat) <- c("Small", "Large")
  colnames(ann_mat) <- c("Small", "Large")
  
  # Return matrix
  ann_mat
}

# ===== Function for calculating the annual matrices based on the vital rates =====
calc_ann_mat <- function(ann_vrs) {
  # List to store the annual matrices
  ann_mats <- list()
  
  # Iterate through the rows
  for(rowi in 1:nrow(ann_vrs)) {
    # Get year
    yr = ann_vrs[[rowi, "Year"]]
    
    # Get vital rates
    rf <- ann_vrs[[rowi, "Recruitment.factor"]]
    ss <- ann_vrs[[rowi, "Small.survival"]]
    sl <- ann_vrs[[rowi, "Large.survival"]]
    g <- ann_vrs[[rowi, "Progression"]]
    r <- ann_vrs[[rowi, "Retrogression"]]
    f <- ann_vrs[[rowi, "Reproduction"]]
    
    # Calculate annual matrix
    ann_mat <- calc_single_mat(
      rf = rf,
      ss = ss,
      sl = sl,
      g = g,
      r = r,
      f = f
    )
    
    # Store matrix
    ann_mats[[as.character(yr)]] <- ann_mat
  }
  
  # Return matrix list
  ann_mats
}

# ===== Function for calculating the 'pooled' matrix across years based on vital rates =====
calc_pooled_mat <- function(ann_vrs) {
  # Calculate total small & large indivs across the years
  small_tot <- sum(ann_vrs$N.small)
  large_tot <- sum(ann_vrs$N.large)
  
  # Final matrix
  pooled_mat <- matrix(rep(0, 4), ncol = 2)
  
  # Iterate through the rows
  for(rowi in 1:nrow(ann_vrs)) {
    # Get year
    yr <- ann_vrs[[rowi, "Year"]]
    
    # Get small & large indiv count
    nsmall <- ann_vrs[[rowi, "N.small"]]
    nlarge <- ann_vrs[[rowi, "N.large"]]
    
    # Get vital rates
    rf <- ann_vrs[[rowi, "Recruitment.factor"]]
    ss <- ann_vrs[[rowi, "Small.survival"]]
    sl <- ann_vrs[[rowi, "Large.survival"]]
    g <- ann_vrs[[rowi, "Progression"]]
    r <- ann_vrs[[rowi, "Retrogression"]]
    f <- ann_vrs[[rowi, "Reproduction"]]
    
    # Calculate annual matrix
    ann_mat <- calc_single_mat(
      rf = rf,
      ss = ss,
      sl = sl,
      g = g,
      r = r,
      f = f
    )
    
    # Weight annual matrix by the number of small and large individuals
    # Multiplier
    mult_mat <- matrix(c(
      nsmall, 0,
      0, nlarge
    ), ncol = 2, byrow = TRUE)
    ann_mat <- mult_mat %*% ann_mat
    
    # Add matrix
    pooled_mat <- pooled_mat + ann_mat
  }
  
  # Divide matrix by the total number of small & large indivs
  # Divider
  div_mat <- matrix(c(
    1 / small_tot, 0,
    0, 1 / large_tot
  ), ncol = 2, byrow = TRUE)
  pooled_mat <- div_mat %*% pooled_mat
  
  # Return average matrix
  pooled_mat
}

# ===== Function for summarising a single matrix & vital rates in a dataframe =====
calc_mat_stats <- function(mat) {
  # Populate the dataframe
  mat_df <- tibble(
    S.S = mat[1,1],
    S.L = mat[2,1],
    L.S = mat[1,2],
    L.L = mat[2,2],
    Lambda = lambda(mat),
    Generation.time = generation.time(mat)
  )
  
  # Return dataframe
  mat_df
}


# ===== Function for summarising the annual matrices & vital rates in a dataframe =====
summ_ann_mat_stats <- function(ann_mats) {
  # Dataframe to store the annual matrices
  ann_mat_df <- tibble(
    Year = numeric(),
    S.S = numeric(),
    S.L = numeric(),
    L.S = numeric(),
    L.L = numeric(),
    Lambda = numeric(),
    Generation.time = numeric()
  )
  
  # Iterate through annual matrices
  for(yr_char in names(ann_mats)) {
    # Get numeric year
    yr <- as.numeric(yr_char)
    
    # Calculate matrix summary
    ann_mat_row <- calc_mat_stats(ann_mats[[yr_char]])
    
    # Insert year
    ann_mat_row$Year <- c(yr)
    
    # Merge row into dataframe
    ann_mat_df <- bind_rows(ann_mat_df, ann_mat_row)
  }
  
  # Return dataframe
  ann_mat_df
}




