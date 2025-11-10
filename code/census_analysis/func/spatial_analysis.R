################################################################################
# Functions for spatial processing
#
# Written by Young Jun Lee
# Jun 2025

# Truncate value between two bounds
bound <- function(x, min_bound, max_bound) {
  min(max(x, min_bound), max_bound)
}

# Function for calculating the area of a circle within bounds
area_circle_in_bound_scalar <- function(x, y, r, xmin, xmax, ymin, ymax, divs = 100) {
  # Simply area of the circle if circle doesn't intersect any bounds
  if(x - r >= xmin & x + r <= xmax & y - r >= ymin & y + r <= ymax) {
    return(pi * r^2)
  } else { # Otherwise, do numerical method
    testgrad <- seq(0, 1 - 1/divs, length.out = divs) + 1/(2*divs) - 0.5
    # Populate points to test
    xtest <- rep(testgrad * 2*r + x, divs)
    ytest <- rep(testgrad * 2*r + y, each = divs)
    # For each test point:
    # Test if point is in circle
    in_area <- ifelse((xtest - x)^2 + (ytest - y)^2 > r^2, 0, 1)
    # Test if point is within bounds
    in_area <- ifelse(xtest < xmin | xtest > xmax |
                          ytest < ymin | ytest > ymax, 0, in_area)
    # Multiply number of points with '1' with unit area
    return(sum(in_area) * (2*r/divs)^2)
  }
}

area_circle_in_bound <- base::Vectorize(area_circle_in_bound_scalar)

# Function for calculating proportional overlap between two circles
prop_overlap_geom <- function(rf, rn, d) {
  if(d >= (rf + rn)) { return(0) }
  if((d + rf) <= rn) { return(1) }
  if((d + rn) <= rf) { return(rn^2 / rf^2) }
  
  h <- sqrt(2 * d^2 * rf^2 + 2 * d^2 * rn^2 + 2 * rf^2 * rn^2
            - d^4 - rf^4 - rn^4) / (2 * d)
  
  # Calculate angles
  angle_f <- if(h / rf >= 1) { pi } else { 2 * asin(h / rf) }
  angle_n <- if(h / rn >= 1) { pi } else { 2 * asin(h / rn) }
  
  # Calculate area of overlap
  a_overlap <- ((angle_f - sin(angle_f)) * rf^2 +
                  (angle_n - sin(angle_n)) * rn^2) / 2
  # Return proportional overlap
  return(a_overlap / (pi * rf^2))
}

# Function for calculating proportional overlap between two circles numerically
# The coordinates & radii of circles are given,
# along with the X and Y coordinates of the boundary.
prop_overlap_num <- function(xf, yf, rf, xn, yn, rn,
                             xmin, xmax, ymin, ymax,
                             divs = 100) {
  # Populate points to test
  xtest <- rep(seq(bound(xf - rf, xmin, xmax), bound(xf + rf, xmin, xmax),
                   length.out = divs), divs)
  ytest <- rep(seq(bound(yf - rf, ymin, ymax), bound(yf + rf, ymin, ymax),
                   length.out = divs), each = divs)
  # For each test point:
  # Test if point is in focal circle
  in_focal <- ifelse((xtest - xf)^2 + (ytest - yf)^2 > rf^2, 0, 1)
  # Test if point is in neighbour circle
  in_area <- ifelse(in_focal != 0 &
                      ((xtest - xn)^2 + (ytest - yn)^2 > rn^2), 0, in_focal)
  # Return proportional overlap
  return(sum(in_area) / sum(in_focal))
}

# Function for calculating proportional overlap between two circles
# contained within the boundary, with automatic selection between the
# two strategies
prop_overlap_scalar <- function(xf, yf, rf, xn, yn, rn,
                         xmin, xmax, ymin, ymax,
                         divs = 100) {
  # Put zero if there is no overlap
  if(sqrt((xf - xn)^2 + (yf - yn)^2) >= rf + rn) {
    0
  } else if(xf - rf < xmin | xf + rf > xmax | yf - rf < ymin | yf + rf > ymax |
     xn - rn < xmin | xn + rn > xmax | yn - rn < ymin | yn + rn > ymax) {
    # Choose numerical if any of the two circles would overlap with the boundary
    prop_overlap_num(xf, yf, rf, xn, yn, rn, xmin, xmax, ymin, ymax, divs)
  } else { # Choose analytical otherwise
    prop_overlap_geom(rf, rn, sqrt((xf - xn)^2 + (yf - yn)^2))
  }
}

# Vectorised version
prop_overlap <- base::Vectorize(prop_overlap_scalar)

# Function for modified 'density'
calc_dens_overlap_pairs <- function(f_data, n_data) {
  f_data_grouped <- f_data %>% group_by(Plot, Year)
  n_data_grouped <- n_data %>% group_by(Plot, Year)
  
  pair_data
  
  pair_data <- do.call(cbind.data.frame, Map(expand.grid,
                                             foc = f_data,
                                             nei = n_data)) %>%
    filter(
      # Select pairs that occur in the same plot in the same year
      (as.character(Year.foc) == as.character(Year.nei) &
         as.character(Plot.foc) == as.character(Plot.nei)) &
        # And are not self-self pairs
        (as.character(Tag.foc) != as.character(Tag.nei))
    )
  
  cat("Data filtered!\n")
  
  pair_data <- pair_data %>%
    mutate(
      Distance..cm. = sqrt(
        (X..cm..foc - X..cm..nei)^2 + (Y..cm..foc - Y..cm..nei)^2
      )
    ) %>%
    mutate( # Calculate 'density'; note that focal and neighbour are swapped
      Prop.overlap = prop_overlap(
        X..cm..nei, Y..cm..nei, Length..cm..nei,
        X..cm..foc, Y..cm..foc, Length..cm..foc,
        0, 200, 0, 200
      ) /
        area_circle_in_bound(X..cm..foc, Y..cm..foc, Length..cm..foc,
                               0, 200, 0, 200) * # Divide by area
        100^2 # Convert to density per square metre
    )
  
  return(pair_data)
}

# Function for calculating total overlap from neighbours given pairwise info
calc_overlap_tot_frompair <- function(pair_data, f_data) {
  # Sum overlap across neighbours
  overlap_data <- pair_data %>%
    ungroup() %>%
    group_by(Year.foc, Tag.foc) %>%
    summarise(
      Total.prop.overlap = sum(Prop.overlap)
    ) %>%
    ungroup()
  
  # Map overlap data to focal dataset
  temp_data <- f_data %>%
    left_join(overlap_data, by = join_by(Year == Year.foc, Tag == Tag.foc)) %>%
    mutate( # Replace NA overlap as this is due to solitary plants
      Total.prop.overlap = ifelse(is.na(Total.prop.overlap), 0, Total.prop.overlap)
    )
  
  # Extract total overlap as vector
  tot_overlap <- temp_data$Total.prop.overlap
  
  return(tot_overlap)
}


# Function for calculating total overlap from neighbours
calc_overlap_tot <- function(f_data, n_data) {
  pair_data <- calc_overlap_pairs(f_data, n_data)
  
  f_data$Total.prop.overlap <- calc_overlap_tot_frompair(pair_data, f_data)
  
  return(f_data)
}

# Convert pairwise overlap data into matrix with row & col labels
calc_overlap_mat <- function(pair_data, f_ids, n_ids) {
  # Get dimensions of overlap matrix
  n_f <- length(f_ids)
  n_n <- length(n_ids)
  
  # Overlap matrix
  m_o <- matrix(0, n_f, n_n)
  rownames(m_o) <- f_ids
  colnames(m_o) <- n_ids
  
  # Insert proportional overlap into matrix
  m_o[cbind(pair_data$ID.foc, pair_data$ID.nei)] <- pair_data$Prop.overlap
  
  # Return matrix
  return(m_o)
}







