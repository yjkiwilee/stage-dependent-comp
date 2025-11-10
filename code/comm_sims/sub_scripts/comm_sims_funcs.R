################################################################################
#' Add-on functions used throughout the scripts
#'
#' Written by Young Jun Lee
#' Mar 2025
#' 

# Linear interpolation function
lerp <- function(x, a, b) {
  a + x * (b - a)
}

# Function for calculating fertility given vital rates & lambda
fert_for_lambda <- function(sj, sa, g, r, l) {
  (l^2 - l*sj + g*l*sj - l*sa + l*r*sa + sj*sa - g*sj*sa - r*sj*sa) /
    (g*sj)
}

# Function for calculating lambda given vital rates
lambda_from_vr <- function(sj, sa, g, r, f) {
  return(
    1/2 * (
      sj - g*sj + sa - r*sa +
        sqrt(
          (-sj + g*sj - sa + r*sa)^2 -
            4 * (
              -f*g*sj + sj*sa - g*sj*sa - r*sj*sa
            )
        )
    )
  )
}

# Function for calculating the sensitivity of lambda to a particular vital rate
d_lambda <- function(foc_vr, vrs) {
  sj <- vrs$s[1]
  sa <- vrs$s[2]
  g <- vrs$g
  r <- vrs$r
  f <- vrs$f
  
  if(foc_vr == "sj") {
    return(
      1/2 * (
        1 - g +
          (-4 * (-f*g + sa - g*sa - r*sa) + 2 * (-1 + g) * (-sj + g*sj - sa + r*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "sa") {
    return(
      1/2 * (
        1 - r +
          (-4 * (sj - g*sj - r*sj) + 2 * (-1 + r) * (-sj + g*sj - sa + r*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "g") {
    return(
      1/2 * (
        -sj +
          (2 * sj * (-sj + g*sj - sa + r*sa) - 4 * (-f*sj - sj*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "r") {
    return(
      1/2 * (
        -sa +
          (4*sj*sa + 2*sa*(-sj + g*sj - sa + r*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "f") {
    return(
      g*sj /
        sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa))
    )
  }
}

d_lambda_vec <- function(foc_vr, sj, sa, g, r, f) {
  if(foc_vr == "sj") {
    return(
      1/2 * (
        1 - g +
          (-4 * (-f*g + sa - g*sa - r*sa) + 2 * (-1 + g) * (-sj + g*sj - sa + r*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "sa") {
    return(
      1/2 * (
        1 - r +
          (-4 * (sj - g*sj - r*sj) + 2 * (-1 + r) * (-sj + g*sj - sa + r*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "g") {
    return(
      1/2 * (
        -sj +
          (2 * sj * (-sj + g*sj - sa + r*sa) - 4 * (-f*sj - sj*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "r") {
    return(
      1/2 * (
        -sa +
          (4*sj*sa + 2*sa*(-sj + g*sj - sa + r*sa)) /
          (2 * sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa)))
      )
    )
  } else if(foc_vr == "f") {
    return(
      g*sj /
        sqrt((-sj + g*sj - sa + r*sa)^2 - 4 * (-f*g*sj + sj*sa - g*sj*sa - r*sj*sa))
    )
  }
}

logistic <- function(x) {
  return(1 / (1 + exp(-x)))
}

logit <- function(x) {
  return(log(x / (1 - x)))
}

# Function for numerically calculating the mean of a logit-normal distribution
mean_logit_norm <- function(norm_mean, norm_sd, n = 1000000, seed = 1) {
  set.seed(seed)
  x <- rnorm(n, norm_mean, norm_sd)
  
  return(mean(logistic(x)))
}

# Function for calculating the mean of the exponential normal
mean_exp_norm <- function(norm_mean, norm_sd) {
  return(exp(norm_mean + norm_sd^2 / 2))
}

# Function for converting vital rates after effect of stochasticity
convert_stoch_vr <- function(vrs, vr_sds) {
  vrs_stoch <- vrs
  
  vrs_stoch$s[1] <- mean_logit_norm(logit(vrs$s[1]), vr_sds$s[1])
  vrs_stoch$s[2] <- mean_logit_norm(logit(vrs$s[2]), vr_sds$s[2])
  vrs_stoch$g <- mean_logit_norm(logit(vrs$g), vr_sds$g)
  vrs_stoch$r <- mean_logit_norm(logit(vrs$r), vr_sds$r)
  vrs_stoch$f <- mean_exp_norm(log(vrs$f), vr_sds$f)
  
  return(vrs_stoch)
}

# Function for calculating the value of a vital rate for a stationary lambda
dd_vr_for_stationary <- function(vrs, dd_vr) {
  sj <- vrs$s[1]
  sa <- vrs$s[2]
  g <- vrs$g
  r <- vrs$r
  f <- vrs$f
  
  dd_vr_val <- NA
  
  if(dd_vr == "sj") {
    dd_vr_val <- (1 + (-1 + r)*sa) /
      (1 + (-1 + r)*sa + g*(-1 + f +sa))
  } else if(dd_vr == "sa") {
    dd_vr_val <- (1 + (-1 + g - f*g)*sj) /
      (1 + r*(-1 + sj) + (-1 + g)*sj)
  } else if(dd_vr == "g") {
    dd_vr_val <- -((1 + (-1 + r)*sa) * (-1 + sj)) /
      (-1 + f + sa)*sj
  } else if(dd_vr == "r") {
    dd_vr_val <- (1 + (-1 + g - f*g)*sj + sa*(-1 + sj - g*sj)) /
      (sa * (-1 + sj))
  } else if(dd_vr == "f") {
    dd_vr_val <- (1 + (-1 + g)*sj + sa*(-1 + r + sj - g*sj - r*sj)) /
      (g * sj)
  }
  
  return(dd_vr_val)
}

build_base_model <- function(sp1_vrs, sp2_vrs,
                             ddvr, dd_sign, w_intra, w_inter,
                             pop_coeffs, vr_sds) {
  
  # Calculate baseline weights
  w_c1_base <- dd_sign * pop_coeffs[1] * w_intra
  w_h1_base <- dd_sign * pop_coeffs[1] * w_inter
  w_c2_base <- dd_sign * pop_coeffs[2] * w_intra
  w_h2_base <- dd_sign * pop_coeffs[2] * w_inter
  
  # Insert into k_mat
  k_mat <- list(
    s = matrix(c(
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0
    ), ncol = 4, byrow = TRUE),
    g = matrix(c( 
      0, 0, 0, 0,
      0, 0, 0, 0
    ), ncol = 4, byrow = TRUE),
    r = matrix(c(
      0, 0, 0, 0,
      0, 0, 0, 0
    ), ncol = 4, byrow = TRUE),
    f = matrix(c(
      0, 0, 0, 0,
      0, 0, 0, 0
    ), ncol = 4, byrow = TRUE)
  )
  
  if(ddvr == "sj") {
    k_mat$s[1,c(1,2)] <- w_c1_base
    k_mat$s[1,c(3,4)] <- w_h1_base
    k_mat$s[3,c(1,2)] <- w_h2_base
    k_mat$s[3,c(3,4)] <- w_c2_base
  } else if(ddvr == "sa") {
    k_mat$s[2,c(1,2)] <- w_c1_base
    k_mat$s[2,c(3,4)] <- w_h1_base
    k_mat$s[4,c(1,2)] <- w_h2_base
    k_mat$s[4,c(3,4)] <- w_c2_base
  } else {
    k_mat[[ddvr]][1,c(1,2)] <- w_c1_base
    k_mat[[ddvr]][1,c(3,4)] <- w_h1_base
    k_mat[[ddvr]][2,c(1,2)] <- w_h2_base
    k_mat[[ddvr]][2,c(3,4)] <- w_c2_base
  }
  
  sp1_vr_sds <- vr_sds
  sp1_vr_sds$s[1] <- sp1_vr_sds$s[1] * w_intra * pop_coeffs[1]
  sp1_vr_sds$s[2] <- sp1_vr_sds$s[2] * w_intra * pop_coeffs[1]
  sp1_vr_sds$g[1] <- sp1_vr_sds$g[1] * w_intra * pop_coeffs[1]
  sp1_vr_sds$r[1] <- sp1_vr_sds$r[1] * w_intra * pop_coeffs[1]
  sp1_vr_sds$f[1] <- sp1_vr_sds$f[1] * w_intra * pop_coeffs[1]
  
  sp2_vr_sds <- vr_sds
  sp2_vr_sds$s[1] <- sp2_vr_sds$s[1] * w_intra * pop_coeffs[2]
  sp2_vr_sds$s[2] <- sp2_vr_sds$s[2] * w_intra * pop_coeffs[2]
  sp2_vr_sds$g[1] <- sp2_vr_sds$g[1] * w_intra * pop_coeffs[2]
  sp2_vr_sds$r[1] <- sp2_vr_sds$r[1] * w_intra * pop_coeffs[2]
  sp2_vr_sds$f[1] <- sp2_vr_sds$f[1] * w_intra * pop_coeffs[2]
  
  # Make community model
  comm_model <-
    def_sdcomp_model(
      2,
      sp1_vrs, sp2_vrs,
      k_mat,
      sp1_vr_sds, sp2_vr_sds
    )
  
  return(comm_model)
}

# Function to calculate population size coefficient given vital rates
calc_pop_coeff <- function(vrs, dd_vr, nt, wc, wh) {
  sj <- vrs$s[1]
  sa <- vrs$s[2]
  g <- vrs$g
  r <- vrs$r
  f <- vrs$f
  
  phi <- NA
  
  if(dd_vr == "sj") { # Juvenile survival is density dependent
    phi <- - ((1 - sa + r*sa) * (-1 + sj)) /
      (g * (-1 + f + sa) * sj)
  } else if(dd_vr == "sa") { # Adult survival is density dependent
    phi <- ((-1 + sa) * (-1 + sj - g*sj + f*g*sj)) /
      (sa * (-r + f*g*sj + r*sj))
  } else if(dd_vr == "g") { # Progression is density dependent
    phi <- ((-1 + g) * (1 - sa + r*sa) * (-1 + sj)) /
      (g * (-1 + sa - r*sa + f*sj + r*sa*sj))
  } else if(dd_vr == "r") { # Retrogression is density dependent
    phi <- ((-1 + r) * (-1 + sa + sj - g*sj + f*g*sj - sa*sj + g*sa*sj)) /
      (r * (-1 + sj - g*sj + f*g*sj + g*sa*sj))
  } else if(dd_vr == "f") { # Fecundity is density dependent
    phi <- (1 + (-1 + g)*sj + sa*(-1 + r + sj - g*sj - r*sj)) /
      (f*g*sj)
  }
  
  return(abs(log(phi) / (nt * (wc + wh))))
}

# Function for generating triangular distribution
rtri <- function(n, min = 0, max = 1, mode = 0.5) {
  mode_trans <- (mode - min) / (max - min)
  c <- sqrt(1 - mode_trans)
  x <- runif(n)
  x <- ifelse(
    x < mode_trans,
    sqrt(mode_trans * x),
    (c - (1 - mode_trans) * sqrt(1 - x)) / c
  )
  x <- x * (max - min) + min
  
  return(x)
}

calc_pop_coeff_num <- function(sp1_vrs, sp2_vrs, dd_vr, dd_sign, nt, wc, wh,
                               vr_sds, l_init = 1.5, r_init = 0.3,
                               n_iter = 15, n_rep = 100, tmin = 300, tmax = 700,
                               sp1_0 = c(nt/2, nt/2), sp2_0 = c(nt/2, nt/2),
                               signif_level = 0.05,
                               sd_fact = 0.005) {
  # Get baseline population coefficients
  sp1_coeff_base <- calc_pop_coeff(sp1_vrs, dd_vr, nt, wc, wh)
  sp2_coeff_base <- calc_pop_coeff(sp2_vrs, dd_vr, nt, wc, wh)
  
  # Calculate 'left' and 'right' baseline coefficients
  sp1_coeff_l <- sp1_coeff_base * l_init
  sp1_coeff_r <- sp1_coeff_base * r_init
  sp2_coeff_l <- sp2_coeff_base * l_init
  sp2_coeff_r <- sp2_coeff_base * r_init
  
  sp1_coeff_m <- NULL
  sp2_coeff_m <- NULL
  sp1_coeff_m_best <- NULL
  sp2_coeff_m_best <- NULL
  sp1_coeff_m_best_dev <- NULL
  sp2_coeff_m_best_dev <- NULL
  sp1_coeff_m_best_n <- 1
  sp2_coeff_m_best_n <- 1
  
  l_mean_sp1 <- NULL
  l_mean_sp2 <- NULL
  r_mean_sp1 <- NULL
  r_mean_sp2 <- NULL
  m_mean_sp1 <- NULL
  m_mean_sp2 <- NULL
  
  t_val <- qt(1 - signif_level / 2, df = n_rep - 1)
  
  for(i in 1:n_iter) {
    if(is.null(l_mean_sp1) | is.null(l_mean_sp2)) {
      l_mean_sp1 <- c()
      l_mean_sp2 <- c()
      
      model_l <- build_base_model(sp1_vrs, sp2_vrs,
                                  dd_vr, dd_sign, wc, wh,
                                  c(sp1_coeff_l, sp2_coeff_l),
                                  vr_sds)
      
      for(j in 1:n_rep) {
        project_l <- sdcomp_project(
          model_l,
          sp1_0, sp2_0,
          timestep = tmax
        )
        l_mean_sp1 <- c(l_mean_sp1, mean(project_l$sp1_s1[tmin:tmax] + project_l$sp1_s2[tmin:tmax]))
        l_mean_sp2 <- c(l_mean_sp2, mean(project_l$sp2_s1[tmin:tmax] + project_l$sp2_s2[tmin:tmax]))
      }
    }
    
    if(is.null(r_mean_sp1) | is.null(r_mean_sp2)) {
      r_mean_sp1 <- c()
      r_mean_sp2 <- c()
      
      model_r <- build_base_model(sp1_vrs, sp2_vrs,
                                  dd_vr, dd_sign, wc, wh,
                                  c(sp1_coeff_r, sp2_coeff_r),
                                  vr_sds)
      
      for(j in 1:n_rep) {
        project_r <- sdcomp_project(
          model_r,
          sp1_0, sp2_0,
          timestep = tmax
        )
        r_mean_sp1 <- c(r_mean_sp1, mean(project_r$sp1_s1[tmin:tmax] + project_r$sp1_s2[tmin:tmax]))
        r_mean_sp2 <- c(r_mean_sp2, mean(project_r$sp2_s1[tmin:tmax] + project_r$sp2_s2[tmin:tmax]))
      }
    }
    
    # Determine if the confidence intervals of l_mean and r_mean overlap with nt
    l_mean_sp1_se <- sd(l_mean_sp1) / sqrt(n_rep)
    l_mean_sp1_margin <- t_val * l_mean_sp1_se
    l_mean_sp1_est <- mean(l_mean_sp1)
    r_mean_sp1_se <- sd(r_mean_sp1) / sqrt(n_rep)
    r_mean_sp1_margin <- t_val * r_mean_sp1_se
    r_mean_sp1_est <- mean(r_mean_sp1)
    l_mean_sp2_se <- sd(l_mean_sp2) / sqrt(n_rep)
    l_mean_sp2_margin <- t_val * l_mean_sp2_se
    l_mean_sp2_est <- mean(l_mean_sp2)
    r_mean_sp2_se <- sd(r_mean_sp2) / sqrt(n_rep)
    r_mean_sp2_margin <- t_val * r_mean_sp2_se
    r_mean_sp2_est <- mean(r_mean_sp2)
    
    sp1_overlaps <- abs(l_mean_sp1_est - nt) < l_mean_sp1_margin |
      abs(r_mean_sp1_est - nt) < r_mean_sp1_margin
    sp2_overlaps <- abs(l_mean_sp2_est - nt) < l_mean_sp2_margin |
      abs(r_mean_sp2_est - nt) < r_mean_sp2_margin
    
    l_dev_sp1 <- (l_mean_sp1_est - nt) / sd(l_mean_sp1)
    l_dev_sp2 <- (l_mean_sp2_est - nt) / sd(l_mean_sp2)
    r_dev_sp1 <- (r_mean_sp1_est - nt) / sd(r_mean_sp1)
    r_dev_sp2 <- (r_mean_sp2_est - nt) / sd(r_mean_sp2)
    
    # Assume linearity between l and r to estimate m
    # sp1_coeff_m_est <-
    #   sp1_coeff_l -
    #   (nt - l_mean_sp1_est) / (r_mean_sp1_est - l_mean_sp1_est) *
    #   (sp1_coeff_l - sp1_coeff_r)
    # sp2_coeff_m_est <-
    #   sp2_coeff_l -
    #   (nt - l_mean_sp2_est) / (r_mean_sp2_est - l_mean_sp2_est) *
    #   (sp2_coeff_l - sp2_coeff_r)
    sp1_coeff_m_est <- (sp1_coeff_r + sp1_coeff_l) / 2
    sp2_coeff_m_est <- (sp2_coeff_r + sp2_coeff_l) / 2
    
    # Get actual m deterministically if ci doesn't overlap, as random otherwise
    
    if(sp1_overlaps) {
      cat("sp1 overlaps!\n")
      
      sp1_coeff_m <- rnorm(
        1,
        sp1_coeff_m_best,
        abs(sp1_coeff_m_best * sd_fact * sp1_coeff_m_best_dev)
      )
    } else {
      sp1_coeff_m <- sp1_coeff_m_est
    }
    if(sp2_overlaps) {
      cat("sp2 overlaps!\n")
      
      sp2_coeff_m <- rnorm(
        1,
        sp2_coeff_m_best,
        abs(sp2_coeff_m_best * sd_fact * sp2_coeff_m_best_dev)
      )
    } else {
      sp2_coeff_m <- sp2_coeff_m_est
    }
    
    model_m <- build_base_model(sp1_vrs, sp2_vrs,
                                dd_vr, dd_sign, wc, wh,
                                c(sp1_coeff_m, sp2_coeff_m),
                                vr_sds)
    
    m_mean_sp1 <- c()
    m_mean_sp2 <- c()
    
    
    for(j in 1:n_rep) {
      project_m <- sdcomp_project(
        model_m,
        sp1_0, sp2_0,
        timestep = tmax
      )
      m_mean_sp1 <- c(m_mean_sp1, mean(project_m$sp1_s1[tmin:tmax] + project_m$sp1_s2[tmin:tmax]))
      m_mean_sp2 <- c(m_mean_sp2, mean(project_m$sp2_s1[tmin:tmax] + project_m$sp2_s2[tmin:tmax]))
    }
    
    m_mean_sp1_est <- mean(m_mean_sp1)
    m_mean_sp2_est <- mean(m_mean_sp2)
    
    m_dev_sp1 <- (m_mean_sp1_est - nt) / sd(m_mean_sp1)
    m_dev_sp2 <- (m_mean_sp2_est - nt) / sd(m_mean_sp2)
    
    m_dev_sp1_best <- NULL
    m_dev_sp2_best <- NULL
    
    if(sp1_overlaps | sp2_overlaps) {
      model_m_best <- build_base_model(sp1_vrs, sp2_vrs,
                                       dd_vr, dd_sign, wc, wh,
                                       c(sp1_coeff_m_best, sp2_coeff_m_best),
                                       vr_sds)
      
      m_mean_sp1_best <- c()
      m_mean_sp2_best <- c()
      
      for(j in 1:n_rep) {
        project_m_best <- sdcomp_project(
          model_m,
          sp1_0, sp2_0,
          timestep = tmax
        )
        m_mean_sp1_best <- c(m_mean_sp1_best,
                             mean(project_m_best$sp1_s1[tmin:tmax] + project_m_best$sp1_s2[tmin:tmax]))
        m_mean_sp2_best <- c(m_mean_sp2_best,
                             mean(project_m_best$sp2_s1[tmin:tmax] + project_m_best$sp2_s2[tmin:tmax]))
      }
      
      m_mean_sp1_est_best <- mean(m_mean_sp1_best)
      m_mean_sp2_est_best <- mean(m_mean_sp2_best)
      
      m_dev_sp1_best <- (m_mean_sp1_est_best - nt) / sd(m_mean_sp1_best)
      m_dev_sp2_best <- (m_mean_sp2_est_best - nt) / sd(m_mean_sp2_best)
    }
    
    if(!sp1_overlaps) { # Update bounds only if CIs don't overlap
      if(m_dev_sp1 > 0 & abs(m_dev_sp1) < abs(r_dev_sp1)) {
        sp1_coeff_r <- sp1_coeff_m
        r_mean_sp1 <- m_mean_sp1
      } else if(m_dev_sp1 < 0 & abs(m_dev_sp1) < abs(l_dev_sp1)) {
        sp1_coeff_l <- sp1_coeff_m
        l_mean_sp1 <- m_mean_sp1
      }
      sp1_coeff_m_best <- sp1_coeff_m
      sp1_coeff_m_best_dev <- m_dev_sp1
    } else {
      if(is.null(sp1_coeff_m_best)) {
        sp1_coeff_m_best <- sp1_coeff_m
        sp1_coeff_m_best_dev <- m_dev_sp1
      } else {
        sp1_coeff_m_best_n <- sp1_coeff_m_best_n + 1
        sp1_coeff_m_best_dev <-
          (sp1_coeff_m_best_dev * (sp1_coeff_m_best_n - 1) +
          m_dev_sp1_best) / sp1_coeff_m_best_n
        
        if(abs(m_dev_sp1) < abs(sp1_coeff_m_best_dev)) {
          cat("sp1 improved!\n")
          sp1_coeff_m_best <- sp1_coeff_m
          sp1_coeff_m_best_dev <- m_dev_sp1
          sp1_coeff_m_best_n <- 1
        }
      }
    }
    
    if(!sp2_overlaps) {
      if(m_dev_sp2 > 0 & abs(m_dev_sp2) < abs(r_dev_sp2)) {
        sp2_coeff_r <- sp2_coeff_m
        r_mean_sp2 <- m_mean_sp2
      } else if(m_dev_sp2 < 0 & abs(m_dev_sp2) < abs(l_dev_sp2)) {
        sp2_coeff_l <- sp2_coeff_m
        l_mean_sp2 <- m_mean_sp2
      }
      sp2_coeff_m_best <- sp2_coeff_m
      sp2_coeff_m_best_dev <- m_dev_sp2
    } else {
      if(is.null(sp2_coeff_m_best)) {
        sp2_coeff_m_best <- sp2_coeff_m
        sp2_coeff_m_best_dev <- m_dev_sp2
      } else {
        sp2_coeff_m_best_n <- sp2_coeff_m_best_n + 1
        sp2_coeff_m_best_dev <-
          (sp2_coeff_m_best_dev * (sp2_coeff_m_best_n - 1) +
             m_dev_sp2_best) / sp2_coeff_m_best_n
        
        if(abs(m_dev_sp2) < abs(sp2_coeff_m_best_dev)) {
          cat("sp2 improved!\n")
          sp2_coeff_m_best <- sp2_coeff_m
          sp2_coeff_m_best_dev <- m_dev_sp2
          sp2_coeff_m_best_n <- 1
        }
      }
    }
    
    cat(paste(
      sp1_coeff_m_best, sp1_coeff_m_best_dev,
      sp2_coeff_m_best, sp2_coeff_m_best_dev,
      "\n"
    ))
  }
  
  if(is.null(sp1_coeff_m_best)) {
    sp1_coeff_m_best <- sp1_coeff_m
  }
  if(is.null(sp2_coeff_m_best)) {
    sp2_coeff_m_best <- sp2_coeff_m
  }
  
  return(c(
    sp1_coeff_m_best,
    sp2_coeff_m_best
  ))
}

# Function to analytically alter stage-dependence
alter_stdep <- function(base_model, sd_fact, ddvr,
                        sp1_eq_struct, sp2_eq_struct) {
  k_mat <- base_model$k_mat
  sd_fact_interm <- (sd_fact + 1) / 2
  sd_fact_trans <- sd_fact_interm / (1 - sd_fact_interm)
  
  w_c1_base <- NA
  w_h1_base <- NA
  w_c2_base <- NA
  w_h2_base <- NA
  
  # Recover base weights
  if(ddvr == "sj") {
    w_c1_base <- k_mat$s[1,1]
    w_h1_base <- k_mat$s[1,3]
    w_h2_base <- k_mat$s[3,1]
    w_c2_base <- k_mat$s[3,3]
  } else if(ddvr == "sa") {
    w_c1_base <- k_mat$s[2,1]
    w_h1_base <- k_mat$s[2,3]
    w_h2_base <- k_mat$s[4,1]
    w_c2_base <- k_mat$s[4,3]
  } else {
    w_c1_base <- k_mat[[ddvr]][1,1]
    w_h1_base <- k_mat[[ddvr]][1,3]
    w_h2_base <- k_mat[[ddvr]][2,1]
    w_c2_base <- k_mat[[ddvr]][2,3]
  }
  
  # Calculate stage-dependent weights
  w_c1_a <- w_c1_base / (sd_fact_trans * sp1_eq_struct[1] + sp1_eq_struct[2])
  w_h1_a <- w_h1_base / (sd_fact_trans * sp2_eq_struct[1] + sp2_eq_struct[2])
  w_c2_a <- w_c2_base / (sd_fact_trans * sp2_eq_struct[1] + sp2_eq_struct[2])
  w_h2_a <- w_h2_base / (sd_fact_trans * sp1_eq_struct[1] + sp1_eq_struct[2])
  
  w_c1_j <- w_c1_base / (sp1_eq_struct[1] + sp1_eq_struct[2] / sd_fact_trans)
  w_h1_j <- w_h1_base / (sp2_eq_struct[1] + sp2_eq_struct[2] / sd_fact_trans)
  w_c2_j <- w_c2_base / (sp2_eq_struct[1] + sp2_eq_struct[2] / sd_fact_trans)
  w_h2_j <- w_h2_base / (sp1_eq_struct[1] + sp1_eq_struct[2] / sd_fact_trans)
  
  # Modify k_mat accordingly
  if(ddvr == "sj") {
    k_mat$s[1,1] <- w_c1_j
    k_mat$s[1,2] <- w_c1_a
    k_mat$s[1,3] <- w_h1_j
    k_mat$s[1,4] <- w_h1_a
    k_mat$s[3,1] <- w_h2_j
    k_mat$s[3,2] <- w_h2_a
    k_mat$s[3,3] <- w_c2_j
    k_mat$s[3,4] <- w_c2_a
  } else if(ddvr == "sa") {
    k_mat$s[2,1] <- w_c1_j
    k_mat$s[2,2] <- w_c1_a
    k_mat$s[2,3] <- w_h1_j
    k_mat$s[2,4] <- w_h1_a
    k_mat$s[4,1] <- w_h2_j
    k_mat$s[4,2] <- w_h2_a
    k_mat$s[4,3] <- w_c2_j
    k_mat$s[4,4] <- w_c2_a
  } else {
    k_mat[[ddvr]][1,1] <- w_c1_j
    k_mat[[ddvr]][1,2] <- w_c1_a
    k_mat[[ddvr]][1,3] <- w_h1_j
    k_mat[[ddvr]][1,4] <- w_h1_a
    k_mat[[ddvr]][2,1] <- w_h2_j
    k_mat[[ddvr]][2,2] <- w_h2_a
    k_mat[[ddvr]][2,3] <- w_c2_j
    k_mat[[ddvr]][2,4] <- w_c2_a
  }
  
  # Insert k_mat
  base_model$k_mat <- k_mat
  
  return(base_model)
}

# Function for generating pulse disturbances
rpulse <- function(n, p, mag) {
  rand_p <- runif(n)
  pulses <- ifelse(
    rand_p < (p/2),
    mag,
    ifelse(
      rand_p > (1 - p/2),
      -mag,
      0
    )
  )
  
  return(pulses)
}

# Function to calculate the SSE of a linear model
sse <- function(mod) {
  return(sum(mod$residuals^2))
}

# Function for calculating partial R squared
partial_rsq <- function(full_mod, red_mod) {
  (
    sse(red_mod) -
      sse(full_mod)
  ) / sse(red_mod)
}

# Function for wrapping plots together for simulation results
wrap_plots_custom <- function(plots) {
  wrap_plots(
    plots$sj, plots$sa,
    plots$g, plots$r,
    plots$f,
    ncol = 2
  ) +
    plot_layout(guides = "collect") +
    guide_area() +
    plot_annotation(
      tag_levels = "A"
    ) &
    theme(plot.tag = element_text(face = 'bold'))
}


# Function for calculating statistical summary at each level of stage dependence
calc_stat_summs <- function(tseries_sub, key_df) {
  # Statistical models to fit to each group
  group_mods_exp <- list(
    ab = "sp1_n + sp2_n",
    struct = "sp1_struct + sp2_struct",
    ab_struct = "sp1_n + sp2_n + sp1_struct + sp2_struct",
    full = "sp1_n + sp2_n + sp1_struct + sp2_struct + sp1_stoch_ddvr + sp2_stoch_ddvr",
    st_ab = "sp1_s1 + sp1_s2 + sp2_s1 + sp2_s2"
  )
  group_mods_resp <- list(
    sp1_dd_vr = "sp1_dd_vr",
    sp2_dd_vr = "sp2_dd_vr",
    sp1_lambda = "sp1_lambda",
    sp2_lambda = "sp2_lambda",
    sp1_real_lambda = "sp1_real_lambda",
    sp2_real_lambda = "sp2_real_lambda"
  )
  stat_names <- c("rsq", "regdf", "resdf", "f", "p")
  cor_stat_names <- c("cor_r", "cor_p")
  
  mod_summs_list <- lapply(
    names(group_mods_resp),
    function(resp_name) {
      res_df_list <- lapply(
        names(group_mods_exp),
        function(exp_name) {
          mod_form <- as.formula(paste0(
            group_mods_resp[[resp_name]],
            " ~ ",
            group_mods_exp[[exp_name]]
          ))
          
          fit_mod <- lm(mod_form, tseries_sub)
          
          mod_summ <- summary(fit_mod)
          
          rsq <- mod_summ$r.squared
          regdf <- mod_summ$fstatistic[2]
          resdf <- mod_summ$fstatistic[3]
          fstat <- mod_summ$fstatistic[1]
          p_val <- pf(fstat, regdf, resdf, lower.tail = FALSE)
          
          stats_vec <- c(rsq, regdf, resdf, fstat, p_val)
          
          names(stats_vec) <- paste0(
            resp_name, "_", exp_name, "_", stat_names
          )
          
          return(as.data.frame(as.list(stats_vec)))
        }
      )
      
      res_df <- bind_cols(res_df_list)
      
      return(res_df)
    }
  )
  
  mod_summs_row_df <- bind_cols(mod_summs_list)
  
  # # Calculate partial R squared
  # partial_names <- c()
  # partial_vals <- c()
  # for(resp_var in names(group_mods_resp)) {
  #   partial_vals <- c(partial_vals, c(
  #     partial_rsq(fit_mods[[resp_var]]$ab_struct, fit_mods[[resp_var]]$struct),
  #     partial_rsq(fit_mods[[resp_var]]$ab_struct, fit_mods[[resp_var]]$ab)
  #   ))
  #   
  #   partial_names <- c(partial_names, paste0(
  #     resp_var, c("_ab_pr", "_st_pr")
  #   ))
  # }
  # names(partial_vals) <- partial_names
  # partial_row_df <- data.frame(as.list(partial_vals))
  
  # Calculate correlations
  cor_tests <- list(
    sp1_ab_struct = cor.test(tseries_sub$sp1_n, tseries_sub$sp1_struct),
    sp1_rl_struct = cor.test(tseries_sub$sp1_real_lambda, tseries_sub$sp1_struct),
    sp1_rl_ab = cor.test(tseries_sub$sp1_real_lambda, tseries_sub$sp1_n),
    sp2_ab_struct = cor.test(tseries_sub$sp2_n, tseries_sub$sp2_struct),
    sp2_rl_struct = cor.test(tseries_sub$sp2_real_lambda, tseries_sub$sp2_struct),
    sp2_rl_ab = cor.test(tseries_sub$sp2_real_lambda, tseries_sub$sp2_n)
  )

  cor_tests_summ <- lapply(names(cor_tests), function(ct_name) {
    ct <- cor_tests[[ct_name]]
    
    cor_vec <- c(ct$estimate, ct$p.value)
    
    names(cor_vec) <- paste0(
      ct_name, "_", cor_stat_names
    )
    
    return(as.data.frame(as.list(cor_vec)))
  })
  
  cor_row_df <- bind_cols(
    tibble(cor_df = nrow(tseries_sub) - 2),
    cor_tests_summ
  )
  
  # Merge dfs
  var_row_df <- tibble(
    n = nrow(tseries_sub)
  ) %>%
    bind_cols(mod_summs_row_df) %>%
    bind_cols(cor_row_df)
  
  return(var_row_df)
}

# Function for getting slope and intercept summaries
calc_ddvr_lm <- function(tseries_sub, key_df) {
  # Statistical models to fit to each group
  group_mods_exp <- list(
    ab = "sp1_n + sp2_n",
    st_ab = "sp1_s1 + sp1_s2 + sp2_s1 + sp2_s2",
    env_only = "",
    ab_only = "sp1_n + sp2_n",
    st_ab_only = "sp1_s1 + sp1_s2 + sp2_s1 + sp2_s2"
  )
  group_mods_resp <- list(
    sp1_dd_vr = "sp1_dd_vr",
    sp2_dd_vr = "sp2_dd_vr"
  )
  
  mod_summs_list <- lapply(
    names(group_mods_resp),
    function(resp_name) {
      lapply(
        names(group_mods_exp),
        function(exp_name) {
          mod_form <- as.formula(paste0(
            group_mods_resp[[resp_name]],
            " ~ ",
            group_mods_exp[[exp_name]],
            ifelse(exp_name %in% c("ab_only", "st_ab_only"),
                   "",
                   ifelse(resp_name == "sp1_dd_vr",
                          " + sp1_stoch_ddvr",
                          " + sp2_stoch_ddvr"))
          ))
          
          fit_mod <- lm(mod_form, tseries_sub)
          
          mod_summ <- summary(fit_mod)
          
          stats_part <- tibble(
            exp_name = exp_name,
            resp_name = resp_name,
            var_name = rownames(mod_summ$coefficients),
            var_est = mod_summ$coefficients[,1],
            var_se = mod_summ$coefficients[,2],
            var_tval = mod_summ$coefficients[,3],
            var_pval = mod_summ$coefficients[,4]
          )
          
          return(as.data.frame(as.list(stats_part)))
        }
      )
    }
  )
  
  mod_summs_row_df <- bind_rows(mod_summs_list)
  
  # Merge dfs
  var_row_df <- mod_summs_row_df %>%
    mutate(
      n = nrow(tseries_sub)
    )
  
  return(var_row_df)
}

# Function for getting estimates & p-values from lm
get_lm_est <- function(mod, coeff_name) {
  mod_summ <- summary(mod)
  mod_summ$coefficients[coeff_name, 1]
}
get_lm_pval <- function(mod, coeff_name) {
  mod_summ <- summary(mod)
  mod_summ$coefficients[coeff_name, 4]
}

# Function for extracting interaction coefficients from community model
get_kmat_section <- function(comm_model, dd_vr, idx) {
  comm_k_mat <- comm_model$k_mat[[dd_vr]]
  n_vrs <- nrow(comm_k_mat)
  return(
    comm_k_mat[c(idx, n_vrs / 2 + idx),]
  )
}

# Function for calculating slopes of logistic & exponential functions given intercept
logistic_slope <- function(x) {
  x_temp <- logistic(x)
  return((1 - temp) * temp)
}
exp_slope <- function(x) {
  return(exp(x))
}

# Function for converting unpacked vrs to list of mpms
vrs_to_mpms_vec <- function(sj, sa, g, r, f) {
  return(lapply(
    1:length(sj),
    function(i) {
      return(matrix(c(
        (1 - g[i]) * sj[i], f[i] + r[i] * sa[i],
        g[i] * sj[i], (1 - r[i]) * sa[i]
      ), ncol = 2, byrow = TRUE))
    }
  ))
}
# Function for converting list of mpms to damping ratios
damping_ratio_list <- function(mpms_list) {
  res_list <- lapply(
    mpms_list,
    function(mpm) { 
      if(NA %in% mpm) {
        return(NA)
      } else {
        return(damping.ratio(mpm))
      }
    }
  )
  
  return(unlist(res_list, use.names = FALSE))
}
damping_ratio_vr <- function(sj, sa, g, r, f) {
  a <- sj - g*sj + sa - r*sa
  b <- sqrt(a^2 - 4*(-f*g*sj + sj*sa - g*sj*sa - r*sj*sa))
  
  return(
    (a + b) / abs(a - b)
  )
}

# Function for calculating projected juvenile and adult abundances from vrs
project_juv_vr <- function(nj, na, sj, sa, g, r, f) {
  return(
    nj * (1 - g) * sj + na * (f + r * sa)
  )
}
project_ad_vr <- function(nj, na, sj, sa, g, r, f) {
  return(
    nj * g * sj + na * (1 - r) * sa
  )
}

# # Function for calculating fertility for some target equilibrium population size
# update_f_for_eq <- function(base_vrs, target_n, dd_vr_name, w_intra, w_inter) {
#   # Translate target population size to target density-dependent effect at equilibrium
#   target_dd <- target_n * (w_intra + w_inter)
#   
#   # Resulting fertility
#   res_f <- NULL
#   res_vrs <- base_vrs
#   
#   # If dd vr is not fertility
#   if(dd_vr_name != "f") {
#     # Calculate the density-dependent vital rate
#     # at the target strength of density-dependent effects
#     if(dd_vr_name == "sj") {
#       res_vrs$s[1] <- logistic(logit(res_vrs$s[1]) +
#                                  dd_sign[[dd_vr_name]] *
#                                  target_dd)
#     } else if(dd_vr_name == "sa") {
#       res_vrs$s[2] <- logistic(logit(res_vrs$s[2]) +
#                                  dd_sign[[dd_vr_name]] *
#                                  target_dd)
#     } else {
#       res_vrs[[dd_vr_name]][1] <- logistic(logit(res_vrs[[dd_vr_name]][1]) +
#                                              dd_sign[[dd_vr_name]] *
#                                              target_dd)
#     }
#     
#     # Calculate fertility such that lambda is equal to 1
#     res_f <- fert_for_lambda(
#       res_vrs$s[1],
#       res_vrs$s[2],
#       res_vrs$g[1],
#       res_vrs$r[1],
#       1
#     )
#   } else { # If dd vr is fertility
#     # Calculate the value of fertility such that
#     # lambda will be 1 at the target population size
#     sj <- res_vrs$s[1]
#     sa <- res_vrs$s[2]
#     g <- res_vrs$g[1]
#     r <- res_vrs$r[1]
#     
#     exp_dd <- exp(target_dd)
#     res_f <- -(
#       exp_dd*(-1 + sa - r*sa + sj - g*sj - sa*sj + g*sa*sj + r*sa*sj)
#     )/(
#       g*sj
#     )
#   }
#   
#   # Return final vital rates
#   res_vrs <- base_vrs
#   res_vrs$f[1] <- res_f
#   return(res_vrs)
# }

# Function for calculating fertility for some target equilibrium population size
update_f_for_eq <- function(base_vrs, target_n, dd_vr_name, w_intra, w_inter) {
  # Translate target population size to target density-dependent effect at equilibrium
  target_dd <- target_n * (w_intra + w_inter)
  
  # Resulting fertility
  res_f <- NULL
  res_vrs <- base_vrs
  
  # If dd vr is not fertility
  if(dd_vr_name != "f") {
    # Calculate the density-dependent vital rate
    # at the target strength of density-dependent effects
    if(dd_vr_name == "sj") {
      res_vrs$s[1] <- mod_logistic(mod_logit(res_vrs$s[1]) +
                                 dd_sign[[dd_vr_name]] *
                                 target_dd)
    } else if(dd_vr_name == "sa") {
      res_vrs$s[2] <- mod_logistic(mod_logit(res_vrs$s[2]) +
                                 dd_sign[[dd_vr_name]] *
                                 target_dd)
    } else {
      res_vrs[[dd_vr_name]][1] <- mod_logistic(mod_logit(res_vrs[[dd_vr_name]][1]) +
                                             dd_sign[[dd_vr_name]] *
                                             target_dd)
    }
    
    # Calculate fertility such that lambda is equal to 1
    res_f <- fert_for_lambda(
      res_vrs$s[1],
      res_vrs$s[2],
      res_vrs$g[1],
      res_vrs$r[1],
      1
    )
  } else { # If dd vr is fertility
    # Calculate the value of fertility such that
    # lambda will be 1 at the target population size
    sj <- res_vrs$s[1]
    sa <- res_vrs$s[2]
    g <- res_vrs$g[1]
    r <- res_vrs$r[1]
    
    # exp_dd <- exp(target_dd)
    # res_f <- (
    #   exp_dd*(-1 + sa - r*sa + sj - g*sj - sa*sj + g*sa*sj + r*sa*sj)
    # )/(
    #   g*sj
    # )
    res_f <- (
      1 - sa + r*sa - sj + g*sj + sa*sj - g*sa*sj - r*sa*sj + g*sj*target_dd
    )/(
      g*sj
    )
  }
  
  # Return final vital rates
  res_vrs <- base_vrs
  res_vrs$f[1] <- res_f
  return(res_vrs)
}

mod_logistic <- function(x) {
  return(pmax(pmin(x, 1), 0))
}
mod_logit <- function(y) {
  return(ifelse(y < 0, NA, ifelse(y > 1, NA, y)))
}
mod_exp <- function(x) {
  return(pmax(x, 0))
}
mod_log <- function(y) {
  return(ifelse(y < 0, NA, y))
}

# Function for calculating the sensitivity of proportion of juveniles to changes in vital rates
d_prop_juv <- function(sj, sa, g, r, f, dd_vr) {
  d_val <- NA
  
  if(dd_vr == "sj") {
    d_val <- (-2*g*((-1 + r)**2*sa**2 + 2*f*g*sj + (-1 + g + r + g*r)*sa*sj + 
                      -      (-1 + r)*sa*sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))))/
      -  (sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))*
            -    ((-1 + r)*sa + sj + g*sj + sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj)))**2)
  } else if(dd_vr == "sa") {
    d_val <- (2*g*sj*((-1 + r)**2*sa + (-1 + g + r + g*r)*sj + 
                        -      (-1 + r)*sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))))/
      -  (sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))*
            -    ((-1 + r)*sa + sj + g*sj + sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj)))**2)
  } else if(dd_vr == "g") {
    d_val <- (-2*sj*((-1 + r)**2*sa**2 + (-2 + g + 2*r + g*r)*sa*sj + 
                       -      (-1 + r)*sa*sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj)) + 
                       -      sj*(2*f*g + sj - g*sj + sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj)))))/
      -  (sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))*
            -    ((-1 + r)*sa + sj + g*sj + sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj)))**2)
  } else if(dd_vr == "r") {
    d_val <-  (2*g*sa*sj)/(sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))*
                             -    ((-1 + r)*sa + sj + g*sj + sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))))
  } else if(dd_vr == "f") {
    d_val <- (4*g**2*sj**2)/(sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj))*
                               -    ((-1 + r)*sa + sj + g*sj + sqrt((-1 + r)**2*sa**2 + 2*(-1 + g + r + g*r)*sa*sj + sj*(4*f*g + (-1 + g)**2*sj)))**2)
  }
  
  return(d_val)
}


filter_na <- function(x, if_na) {
  return(
    ifelse(is.na(x), if_na, x)
  )
}

