# ===== Functions for spatial processing =====

# Function to visualist a specific plot in a given year
plt_yr_plot <- function(data_sub) {
  # Plot points
  ggplot() +
    geom_circle(data = data_sub,
                aes(x0 = X..cm., y0 = Y..cm., r = Length..cm., color = Taxon, fill = Taxon),
                alpha = 0.5) +
    scale_radius() +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200), expand = FALSE) +
    theme_bw() +
    theme(aspect.ratio = 1)
}

# Function to visualist a specific plot in a given year
plt_yr_plot_sqrt <- function(data_sub) {
  # Plot points
  ggplot() +
    geom_circle(data = data_sub,
                aes(x0 = X..cm., y0 = Y..cm., r = sqrt(Length..cm.), color = Taxon, fill = Taxon),
                alpha = 0.5) +
    scale_radius() +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200), expand = FALSE) +
    theme_bw() +
    theme(aspect.ratio = 1)
}

# Function to visualist a specific plot in a given year with constant radius
plt_yr_plot_const_r <- function(data_sub, r) {
  # Plot points
  ggplot() +
    geom_circle(data = data_sub,
                aes(x0 = X..cm., y0 = Y..cm., r = r, color = Taxon, fill = Taxon),
                alpha = 0.5) +
    scale_radius() +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200), expand = FALSE) +
    theme_bw() +
    theme(aspect.ratio = 1)
}

# Function to construct a matrix of distances in a given plot
calc_dist_mat <- function(data_sub) {
  # Extract coordinates only
  data_coords <- data_sub %>%
    dplyr::select(X..cm., Y..cm.)
  
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(data_coords, diag = TRUE, upper = TRUE), labels = TRUE)
  
  # Rename rows and columns with Tag
  colnames(dist_mat) <- rownames(dist_mat) <- data_sub$Tag
  
  # Return matrix
  dist_mat
}

# Function to construct distance matrix for Bayesian fitting
calc_dist_mat_across <- function(data_row, data_col) {
  # Get dimensions of distance matrix
  n_r <- nrow(data_row)
  n_c <- nrow(data_col)
  
  # Add unique key column
  data_row <- data_row %>%
    mutate(Key = paste(Year, Plot, Tag, sep = "/"))
  data_col <- data_col %>%
    mutate(Key = paste(Year, Plot, Tag, sep = "/"))
  
  # Distance matrix
  dm_res <- matrix(Inf, n_r, n_c)
  rownames(dm_res) <- data_row$Key
  colnames(dm_res) <- data_col$Key
  
  # Iterate through years
  for(yr in unique(data_row$Year)) {
    # Subset data
    data_r_yr <- data_row %>%
      filter(Year == yr)
    data_c_yr <- data_col %>%
      filter(Year == yr)
    
    # Iterate through plots
    for(plt in unique(data_r_yr$Plot)) {
      # Subset data
      data_r_plt <- data_r_yr %>%
        filter(Plot == plt)
      data_c_plt <- data_c_yr %>%
        filter(Plot == plt)
      
      # Skip if either of the two subsets are empty
      if(nrow(data_r_plt) == 0 || nrow(data_c_plt) == 0) { next }
      
      # Merge subset data & calculate distance matrix
      dm_plt <- calc_dist_mat(bind_rows(data_r_plt, data_c_plt))
      
      # Select only rows in data_r_plt
      dm_plt <- dm_plt[data_r_plt$Tag, , drop = FALSE]
      # Select only cols in data_c_plt
      dm_plt <- dm_plt[, data_c_plt$Tag, drop = FALSE]
      
      # Set self-distance as Inf only if there are common tags between row and col
      if(length(intersect(data_r_plt$Tag, data_c_plt$Tag)) > 0) {
        dm_plt[cbind(
          data_r_plt$Tag,
          intersect(data_r_plt$Tag, data_c_plt$Tag)
        )] <- Inf
      }
      
      # Convert row & colnames into Key
      rownames(dm_plt) <- paste(yr, plt, rownames(dm_plt), sep = "/")
      colnames(dm_plt) <- paste(yr, plt, colnames(dm_plt), sep = "/")
      
      # Insert distance values into main distance matrix
      row_map <- match(rownames(dm_plt), rownames(dm_res))
      col_map <- match(colnames(dm_plt), colnames(dm_res))
      dm_res[row_map, col_map] <- dm_plt
    }
  }
  
  # Return distance matrix
  return(dm_res)
}

# Function to convert the distance matrix into a matrix of neighbourhood status with CONSTANT R
calc_neighbour_mat_constr <- function(dist_mat, r) {
  # Set elements with distance less than or equal to 2 * r as neighbours
  dist_mat <- ifelse(dist_mat <= 2 * r, yes = 1, no = 0)
  # Return matrix
  dist_mat
}

# Function to convert the distance matrix into a matris of neighbourhood status
# where individuals are neighbours if d <= r1 + r2 where r1 and r2 are the lengths of each individual
calc_neighbour_mat_dynr <- function(dist_mat, data_sub) {
  # Extract plant lengths from data_sub as column vector
  plant_lens <- matrix(data_sub$Length..cm., ncol = 1)
  
  # Calculate matrix of r1 + r2
  one <- matrix(rep(1, nrow(plant_lens)), nrow = 1)
  temp_mat <- plant_lens %*% one
  rsum_mat <- temp_mat + t(temp_mat)
  
  # Compare dist_mat against rsum_mat to determine neighbourhood status
  res_mat <- ifelse(dist_mat <= rsum_mat, yes = 1, no = 0)
  # Return resulting matrix
  res_mat
}

# Internal function to collapse the crowding matrix
calc_crowd_row <- function(crowd_mat, tags) {
  if(length(tags) == 0) { 0 }
  else if(length(tags) == 1) { crowd_mat[tags, , drop = FALSE] }
  else { colSums(crowd_mat[tags, , drop = FALSE]) }
}

# Function to convert the distance matrix into a matrix of crowding factor given some shape coefficient for the decaying kernel
# sp1_s_tags, sp1_l_tags, sp2_s_tags, sp2_l_tags are vectors storing the tags corresponding to a particular focal species & size class
# The default function (decay_fun) is the gaussian function
calc_crowd_mat <- function(dist_mat, alpha, sp1_s_tags, sp1_l_tags, sp2_s_tags, sp2_l_tags,
                                     decay_fun) {
  # Set diagonal distances as infinite to remove effect of self
  diag(dist_mat) <- Inf
  
  # Calculate crowding factor with the given function
  crowd_mat <- decay_fun(dist_mat, alpha)
  
  # Determine tags that don't belong to any of the categories
  other_tags <- setdiff(rownames(dist_mat), c(sp1_s_tags, sp1_l_tags, sp2_s_tags, sp2_l_tags))
  
  # Calculate species/stage-specific crowding factors
  sp1_s_c <- calc_crowd_row(crowd_mat, sp1_s_tags) # Select rows & collapse by taking the sum
  sp1_l_c <- calc_crowd_row(crowd_mat, sp1_l_tags)
  sp2_s_c <- calc_crowd_row(crowd_mat, sp2_s_tags)
  sp2_l_c <- calc_crowd_row(crowd_mat, sp2_l_tags)
  other_c <- calc_crowd_row(crowd_mat, other_tags)
  
  # Merge results in a matrix
  res_crowd_mat <- rbind(sp1_s_c, sp1_l_c, sp2_s_c, sp2_l_c, other_c)
  
  # Set rownames
  rownames(res_crowd_mat) <- c("Sp1.S.crowding", "Sp1.L.crowding", "Sp2.S.crowding", "Sp2.L.crowding", "Other.crowding")
  
  # Return result
  res_crowd_mat
}

# Function to convert the neighbourhood status matrix into a dataframe summarising the species identity of neighbours
summ_nei_sp <- function(nei_mat, data_sub) {
  # Copy neighbourhood matrix
  nei_mat_temp <- nei_mat
  
  # Dataframe to store neighbour instances
  neis <- tibble(sp1 = character(), sp2 = character())
  
  # Iterate through matrix & add neighbour species
  for(r_i in data_sub$Tag) {
    for(c_i in data_sub$Tag) {
      # Skip cases where row == column
      if(r_i == c_i) { next }
      # If individuals are neighbours AND this pair hasn't been accounted for yet
      if(nei_mat_temp[r_i,c_i] == 1 & !is.na(nei_mat_temp[c_i,r_i])) {
        # Get species from dataset
        sp1 <- as.character(data_sub[data_sub$Tag == r_i, "Taxon"])
        sp2 <- as.character(data_sub[data_sub$Tag == c_i, "Taxon"])
        # Add species to dataset
        neis <- bind_rows(neis,
                          tibble(sp1 = c(sp1), sp2 = c(sp2)))
        # Set nei_mat_temp value to NA
        nei_mat_temp[r_i,c_i] <- NA
      }
    }
  }
  
  # Get species counts in the dataset
  sp_counts <- data_sub %>%
    group_by(Taxon) %>%
    count() %>%
    pivot_wider(names_from = Taxon, values_from = n) %>%
    as.list()
  
  # Run count() and return result
  neis %>%
    group_by(sp1, sp2) %>% # Grouping ensures that unique combinations aren't merged
    count() %>%
    mutate(
      sp1_n = ifelse(nrow(neis) > 0, sp_counts[[sp1]], numeric()),
      sp2_n = ifelse(nrow(neis) > 0, sp_counts[[sp2]], numeric())
    )
}

# Functions that convert data for a plot in a given year into a summary of species neighbours

# Dynamic radius
summ_nei_sp_dyn <- function(data_sub) {
  # Calculate the distance matrix
  dist_mat <- calc_dist_mat(data_sub)
  
  # Calcualte matrix of neighbourhood status
  nei_mat <- calc_neighbour_mat_dynr(dist_mat, data_sub)
  
  # Get summary of neighbourhood status
  summ_nei <- summ_nei_sp(nei_mat, data_sub)
  
  # Return the summary
  summ_nei
}

# Constant radius
summ_nei_sp_const <- function(data_sub, r) {
  # Calculate the distance matrix
  dist_mat <- calc_dist_mat(data_sub)
  
  # Calcualte matrix of neighbourhood status
  nei_mat <- calc_neighbour_mat_constr(dist_mat, r)
  
  # Get summary of neighbourhood status
  summ_nei <- summ_nei_sp(nei_mat, data_sub)
  
  # Return the summary
  summ_nei
}

# Function to add neighbourhood line layer to plot
add_nei_lines <- function(plt, data_sub, nei_mat) {
  # Iterate through elements in the neighbourhood matrices to add lines where necessary
  # Dataframe to store line coordinates
  line_coords <- tibble(x1 = numeric(), y1 = numeric(), x2 = numeric(), y2 = numeric())
  # Add neighbourhood lines to dataframe
  for(r_i in data_sub$Tag) {
    for(c_i in data_sub$Tag) {
      # Skip cases where row == column
      if(r_i == c_i) { next }
      # If individuals are neighbours
      if(nei_mat[r_i,c_i] == 1) {
        # Get coordinates from dataset
        r_coords <- as.numeric(data_sub[data_sub$Tag == r_i, c("X..cm.", "Y..cm.")])
        c_coords <- as.numeric(data_sub[data_sub$Tag == c_i, c("X..cm.", "Y..cm.")])
        coords <- c(r_coords, c_coords)
        names(coords) <- c("x1", "y1", "x2", "y2")
        # Add coordinates to dataframe
        line_coords <- bind_rows(line_coords, coords)
      }
    }
  }
  # Add line segment layer to plot
  plt +
    geom_segment(data = line_coords, aes(x = x1, y = y1, xend = x2, yend = y2), color = "black", alpha = 0.3, linewidth = 0.5)
}

# Function for calculating the crowding coefficient for individuals in a plot
# NB: THE DATASET MUST ALREADY BE OF A PLOT IN A GIVEN YEAR
# Default function used is exponential
calc_crowd_coeff <- function(sub_data, sp1, sp2, sp1_thres, sp2_thres, alpha,
                             decay_fun = function(x, a) exp(1)^(-a*x)) {
  # Calculate distance matrix from the plot data
  dist_mat <- calc_dist_mat(sub_data)
  
  # Subset data to only contain either species 1 or 2
  sub_sp1_data <- sub_data %>% filter(Taxon == sp1)
  sub_sp2_data <- sub_data %>% filter(Taxon == sp2)
  
  # Extract the tags of live individuals of focal species & non-focal species
  sp1_s_tags <- unique(filter(sub_sp1_data, Died.this.census.final == 0 & Length..cm. <= sp1_thres))$Tag
  sp1_l_tags <- unique(filter(sub_sp1_data, Died.this.census.final == 0 & Length..cm. > sp1_thres))$Tag
  sp2_s_tags <- unique(filter(sub_sp2_data, Died.this.census.final == 0 & Length..cm. <= sp2_thres))$Tag
  sp2_l_tags <- unique(filter(sub_sp2_data, Died.this.census.final == 0 & Length..cm. > sp2_thres))$Tag
  
  # Calculate crowding matrix between individuals
  crowd_mat <- calc_crowd_mat(dist_mat, alpha,
                              sp1_s_tags, sp1_l_tags,
                              sp2_s_tags, sp2_l_tags,
                              decay_fun)
  
  # Result df to store the crowding factor alongside the existing variables
  res_data <- sub_data %>%
    mutate( # Add columns to store crowding factor
      Sp1.crowding = rep(NA, nrow(.)),
      Sp1.S.crowding = rep(NA, nrow(.)),
      Sp1.L.crowding = rep(NA, nrow(.)),
      Sp2.crowding = rep(NA, nrow(.)),
      Sp2.S.crowding = rep(NA, nrow(.)),
      Sp2.L.crowding = rep(NA, nrow(.)),
      Other.crowding = rep(NA, nrow(.)),
      Total.crowding = rep(NA, nrow(.))
    )
  
  # Insert crowding factor into the dataset
  # Iterate over columns of the crowding factor matrix
  for(tag in colnames(crowd_mat)) {
    # Move over the crowding factor to the appropriate element in the df
    res_data[res_data$Tag == tag, "Sp1.crowding"] <- crowd_mat["Sp1.S.crowding", tag] + crowd_mat["Sp1.L.crowding", tag]
    res_data[res_data$Tag == tag, "Sp1.S.crowding"] <- crowd_mat["Sp1.S.crowding", tag]
    res_data[res_data$Tag == tag, "Sp1.L.crowding"] <- crowd_mat["Sp1.L.crowding", tag]
    res_data[res_data$Tag == tag, "Sp2.crowding"] <- crowd_mat["Sp2.S.crowding", tag] + crowd_mat["Sp2.L.crowding", tag]
    res_data[res_data$Tag == tag, "Sp2.S.crowding"] <- crowd_mat["Sp2.S.crowding", tag]
    res_data[res_data$Tag == tag, "Sp2.L.crowding"] <- crowd_mat["Sp2.L.crowding", tag]
    res_data[res_data$Tag == tag, "Other.crowding"] <- crowd_mat["Other.crowding", tag]
    res_data[res_data$Tag == tag, "Total.crowding"] <- sum(c(crowd_mat[, tag]))
  }
  
  # Return df
  res_data
}

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

# Function for calculating pairwise overlap of neighbours on focal individuals
calc_overlap_pairs <- function(f_data, n_data) {
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
  
  pair_data <- pair_data %>%
    mutate(
      Distance..cm. = sqrt(
        (X..cm..foc - X..cm..nei)^2 + (Y..cm..foc - Y..cm..nei)^2
      )
    ) %>%
    mutate( # Calculate proportional overlap between indivs
      Prop.overlap = prop_overlap(
        X..cm..foc, Y..cm..foc, Length..cm..foc,
        X..cm..nei, Y..cm..nei, Length..cm..nei,
        0, 200, 0, 200
      )
    )
  
  return(pair_data)
}

# Function for calculating binary overlap of neighbours on focal individuals
calc_bin_overlap_pairs <- function(f_data, n_data) {
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
  
  pair_data <- pair_data %>%
    mutate(
      Distance..cm. = sqrt(
        (X..cm..foc - X..cm..nei)^2 + (Y..cm..foc - Y..cm..nei)^2
      )
    ) %>%
    mutate( # Calculate overlap between indivs
      Prop.overlap = ifelse(
        Distance..cm. <= Length..cm..foc + Length..cm..nei,
        1,
        0
      )
    )
  
  return(pair_data)
}

# Function for calculating distance-scaled binary overlap
calc_dist_overlap_pairs <- function(f_data, n_data) {
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
  
  pair_data <- pair_data %>%
    mutate(
      Distance..cm. = sqrt(
        (X..cm..foc - X..cm..nei)^2 + (Y..cm..foc - Y..cm..nei)^2
      )
    ) %>%
    mutate( # Calculate 'overlap' between indivs
      Prop.overlap = ifelse(
        Distance..cm. < (Length..cm..foc + Length..cm..nei),
        1 - (Distance..cm. / (Length..cm..foc + Length..cm..nei)),
        0
      )
    )
  
  return(pair_data)
}

# Function for calculating normalised distance-scaled binary overlap
calc_dist_overlap_pairs_norm <- function(f_data, n_data) {
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
  
  pair_data <- pair_data %>%
    mutate(
      Distance..cm. = sqrt(
        (X..cm..foc - X..cm..nei)^2 + (Y..cm..foc - Y..cm..nei)^2
      )
    ) %>%
    mutate( # Calculate 'overlap' between indivs
      Prop.overlap = ifelse(
        Distance..cm. < (Length..cm..foc + Length..cm..nei),
        (1 - (Distance..cm. / (Length..cm..foc + Length..cm..nei))) / 
          (Length..cm..foc + Length..cm..nei) ^ 2,
        0
      )
    )
  
  return(pair_data)
}

# Function for calculating 'overlap' assuming individuals have the same size
calc_fixed_dist_overlap_pairs <- function(f_data, n_data, assum_size) {
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
  
  pair_data <- pair_data %>%
    mutate(
      Distance..cm. = sqrt(
        (X..cm..foc - X..cm..nei)^2 + (Y..cm..foc - Y..cm..nei)^2
      )
    ) %>%
    mutate( # Calculate 'overlap' between indivs
      Prop.overlap = ifelse(
        Distance..cm. < 2 * assum_size,
        1 - (Distance..cm. / (2 * assum_size)),
        0
      )
    )
  
  return(pair_data)
}

# Function for modified 'density'
calc_dens_overlap_pairs <- function(f_data, n_data) {
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







