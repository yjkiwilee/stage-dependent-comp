################################################################################
# Functions for fitting statistical models to the data
#
# Written by Young Jun Lee
# Nov 2024

# Function for fitting statistical data
fit_vr_glmm <- function(model_list, stage, vr_data, fec_family = "poisson", maxit = 7e4) {
  
  fit_models <- lapply(names(model_list), function(vr) {
    # Skip if the vital rate is not applicable to the life stage
    if(stage == "S" & vr %in% c("Retrogression","X..Capitulescences")) { return(NA) }
    if(stage == "L" & vr %in% c("Progression")) { return(NA) }
    
    # Extract models
    models <- model_list[[vr]]
    
    # Determine whether variable is binary or continuous
    var_type <- if(vr == "X..Capitulescences") {fec_family} 
      else if(vr == "Growth.future") {"gaussian"}
      else {"binomial"}
    
    # If vital rate is Growth.future, remove individuals that die in the next year
    if(vr == "Growth.future") {
      vr_data <- vr_data %>%
        filter(
          Survival.future == 1
        )
    }
    
    # Negative binomial fitting
    if(var_type == "nb") {
      # Fit glmm
      models_temp <- lapply(names(models), function(stat_model_name) {
        cat(paste0("Fitting ", stat_model_name, " to ", stage, ", ", vr, "..."))
        
        stat_model <- models[[stat_model_name]]
        
        vr_data_temp <- if(grepl("^M[0-9].LPREV", stat_model_name)) {
          filter(vr_data, !is.na(Length..cm..prev))
        } else {vr_data}
        
        res_model <- tryCatch({
          glmer.nb(
            stat_model,
            data = vr_data_temp,
            control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = maxit))
          )},
          error = function(e) {
            print(e)
            NA
          })
        
        cat(" done!\n")
        
        return(res_model)
      })
      
      # Put names back & return list
      names(models_temp) <- names(models)
      models_temp
    } else if(var_type == "gaussian") {
      # Fit lmm
      models_temp <- lapply(names(models), function(stat_model_name) {
        cat(paste0("Fitting ", stat_model_name, " to ", stage, ", ", vr, "..."))
        
        stat_model <- models[[stat_model_name]]
        
        vr_data_temp <- if(grepl("^M[0-9].LPREV", stat_model_name)) {
          filter(vr_data, !is.na(Length..cm..prev))
        } else {vr_data}
        
        res_model <- tryCatch({
          lmer(
            stat_model,
            data = vr_data_temp,
            control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = maxit))
          )},
          error = function(e) {
            print(e)
            NA
          })
        
        cat(" done!\n")
        
        return(res_model)
      })
      
      # Put names back & return list
      names(models_temp) <- names(models)
      models_temp
    } else {
      # Fit glmm
      models_temp <- lapply(names(models), function(stat_model_name) {
        cat(paste0("Fitting ", stat_model_name, " to ", stage, ", ", vr, "..."))
        
        stat_model <- models[[stat_model_name]]
        
        vr_data_temp <- if(grepl("^M[0-9].LPREV", stat_model_name)) {
          filter(vr_data, !is.na(Length..cm..prev))
        } else {vr_data}
        
        res_model <- tryCatch({
          glmer(
            stat_model,
            family = var_type,
            data = vr_data_temp,
            control = glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = maxit))
          )},
          error = function(e) {
            print(e)
            NA
          })
        
        cat(" done!\n")
        
        return(res_model)
      })
      
      # Put names back & return list
      names(models_temp) <- names(models)
      models_temp
    }
  })
  # Set names
  names(fit_models) <- names(model_list)
  # Return
  fit_models
}

# Function for fitting LMM
fit_vr_lmm <- function(model_list, stage, vr_data, maxit = 7e4) {
  
  fit_models <- lapply(names(model_list), function(vr) {
    # Skip if the vital rate is not applicable to the life stage
    if(stage == "S" & vr %in% c("Retrogression","X..Capitulescences")) { return(NA) }
    if(stage == "L" & vr %in% c("Progression")) { return(NA) }
    
    # Extract models
    models <- model_list[[vr]]
    
    # If vital rate is Growth.future, remove individuals that die in the next year
    if(vr == "Growth.future") {
      vr_data <- vr_data %>%
        filter(
          Survival.future == 1
        )
    }
    
    # Fit lmm
    models_temp <- lapply(names(models), function(stat_model_name) {
      cat(paste0("Fitting ", stat_model_name, " to ", stage, ", ", vr, "..."))
      
      stat_model <- models[[stat_model_name]]
      
      vr_data_temp <- if(grepl("^M[0-9].LPREV", stat_model_name)) {
        filter(vr_data, !is.na(Length..cm..prev))
      } else {vr_data}
      
      res_model <- tryCatch({
        lmer(
          stat_model,
          data = vr_data_temp,
          control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = maxit))
        )},
        error = function(e) {
          print(e)
          NA
        })
      
      cat(" done!\n")
      
      return(res_model)
    })
    
    # Put names back & return list
    names(models_temp) <- names(models)
    models_temp
  })
  # Set names
  names(fit_models) <- names(model_list)
  # Return
  fit_models
}

# Function for fitting statistical data without random effects
# fit_vr_models_worand <- function(model_list, stage, vr_data, fec_family = "poisson") {
#   fit_models <- lapply(names(model_list), function(vr) {
#     # Skip if the vital rate is not applicable to the life stage
#     if(stage == "S" & vr %in% c("Retrogression","X..Capitulescences")) { return(NA) }
#     if(stage == "L" & vr %in% c("Progression")) { return(NA) }
#     
#     # Extract models
#     models <- model_list[[vr]]
#     
#     # Determine whether variable is binary or continuous
#     var_type <- if(vr == "X..Capitulescences") {fec_family} 
#       else if(vr == "Growth.future") {"gaussian"}
#       else {"binomial"}
#     
#     # Negtaive binomial fitting
#     if(var_type == "nb") {
#       # Fit glmm
#       lapply(models, function(stat_model) {
#         glm.nb(
#           stat_model,
#           data = vr_data
#         )
#       })
#     } else if(var_type == "gaussian") {
#       lapply(models, function(stat_model) {
#         lm(
#           stat_model,
#           data = vr_data
#         )
#       })
#     } else {
#       # Fit glmm
#       lapply(models, function(stat_model) {
#         glm(
#           stat_model,
#           family = var_type,
#           data = vr_data
#         )
#       })
#     }
#   })
#   # Set names
#   names(fit_models) <- names(model_list)
#   # Return
#   fit_models
# }

# Function for fitting data while accounting for spatial structure
# fit_vr_models_spat <- function(model_list, stage, vr_data, fec_family = "poisson") {
#   fit_models <- lapply(names(model_list), function(vr) {
#     # Skip if the vital rate is not applicable to the life stage
#     if(stage == "S" & vr %in% c("Retrogression","X..Capitulescences")) { return(NA) }
#     if(stage == "L" & vr %in% c("Progression")) { return(NA) }
#     
#     # Extract models
#     models <- model_list[[vr]]
#     
#     # Determine whether variable is binary or continuous
#     var_type <- if(vr == "X..Capitulescences") {fec_family} 
#     else if(vr == "Growth.future") {"gaussian"}
#     else {"binomial"}
#     
#     # Negative binomial fitting
#     if(var_type == "nb") {
#       # Fit glmm
#       lapply(models, function(stat_model) {
#         lme(
#           stat_model,
#           data = vr_data,
#           random = ~ 1|Plot,
#           control = glmerControl(optCtrl = list(maxfun = 3e4))
#         )
#       })
#     } else if(var_type == "gaussian") {
#       lapply(models, function(stat_model) {
#         lmer(
#           stat_model,
#           data = vr_data,
#           control = lmerControl(optCtrl = list(maxfun = 3e4))
#         )
#       })
#     } else {
#       # Fit glmm
#       lapply(models, function(stat_model) {
#         glmer(
#           stat_model,
#           family = var_type,
#           data = vr_data,
#           control = glmerControl(optCtrl = list(maxfun = 3e4))
#         )
#       })
#     }
#   })
#   # Set names
#   names(fit_models) <- names(model_list)
#   # Return
#   fit_models
# }

# Function for Bayesian fitting
fit_vr_models_b <- function(model_list, stage, vr_data, fec_family = "poisson") {
  fit_models <- lapply(names(model_list), function(vr) {
    # Skip if the vital rate is not applicable to the life stage
    if(stage == "S" & vr %in% c("Retrogression","X..Capitulescences")) { return(NA) }
    if(stage == "L" & vr %in% c("Progression")) { return(NA) }
    
    # Extract models
    models <- model_list[[vr]]
    
    # Determine whether variable is binary or continuous
    var_type <- if(vr == "X..Capitulescences") {
        if(fec_family == "poisson") {poisson}
        else if(fec_family == "nb") {negbinomial()}
      } 
      else if(vr == "Growth.future") {gaussian()}
      else {bernoulli()}
    
    # Fit glmm
    lapply(models, function(stat_model) {
      brm(
        stat_model,
        family = var_type,
        data = vr_data,
        chains = 2, # nb of chains
        iter = 5000, # nb of iterations, including burnin
        warmup = 1000, # burnin
        thin = 1
      )
    })
  })
  # Set names
  names(fit_models) <- names(model_list)
  # Return
  fit_models
}

# Function for determining whether model has converged
# https://stackoverflow.com/questions/72090177/how-can-i-know-whether-the-model-is-converged-or-failed-to-converge-in-lme4-with
has_converged <- function(mm) {
  if ( !inherits(mm, "merMod")) stop("Error: must pass a lmerMod object")
  
  retval <- NULL
  
  if(is.null(unlist(mm@optinfo$conv$lme4))) {
    retval = 1
  }
  else {
    if (isSingular(mm)) {
      retval = 0
    } else {
      retval = -1
    }
  }
  return(retval)
}


# Function to calculate the summary statistics of each variable
calc_var_summ_stat <- function(df, cols) {
  # Summary df
  std_op_df <- NULL
  
  # For each column, calculate summary statistic
  for(var_col in cols) {
    std_op_df_row <- tibble(
      Variable = var_col,
      Mean = mean(df[[var_col]]),
      Std.dev = sd(df[[var_col]])
    )
  }
}


# Function for extracting summary values from the fitted models
get_stat_summ_df <- function(fit_model_list, confint_method = "Wald") {
  # Resulting df
  stat_summ_df <- NULL
  
  for(sp_stage in names(fit_model_list)) {
    for(vr in names(fit_model_list[[sp_stage]])) {
      for(model_name in names(fit_model_list[[sp_stage]][[vr]])) {
        cat("Processing ", sp_stage, " ", vr, " ", model_name, "\n")
        # Get model
        focal_model <- fit_model_list[[sp_stage]][[vr]][[model_name]]
        
        # Expected list of variables for each model
        exp_vars <- list(
          M1.L.YR = c("(Intercept)", "Length..cm."),
          M2.L.YR = c("(Intercept)", "Length..cm.", "Total.crowding"),
          M3.L.YR = c("(Intercept)", "Length..cm.", "Sp1.crowding", "Sp2.crowding"),
          M4.L.YR = c("(Intercept)", "Length..cm.",
                      "Sp1.L.crowding", "Sp1.S.crowding",
                      "Sp2.L.crowding", "Sp2.S.crowding")
        )
        
        # Determine if model failed to fit or column dropping occurred
        model_failed <- is.na(focal_model) || length(setdiff(
          names(fixef(focal_model)), exp_vars[[model_name]]
        )) > 0
        
        # Exit if model fitting has failed (indicated by NA)
        if(model_failed) {
          # Row of the df
          stat_summ_row <- tibble(
            Taxon = str_split(sp_stage, "[.]")[[1]][1],
            Stage = str_split(sp_stage, "[.]")[[1]][2],
            Vital.rate = vr,
            Model.name = model_name,
            Has.converged = NA,
            Is.singular = NA,
            AIC = NA,
            BIC = NA,
            Random.variance = NA
          )
          # Cross-join into single df
          stat_summ_part <- cross_join(stat_summ_row, fixefs_df)
          # Add to df
          if(is.null(stat_summ_df)) {
            stat_summ_df <- stat_summ_part
          } else {
            stat_summ_df <- bind_rows(stat_summ_df, stat_summ_part)
          }
          
          # Continue to next model
          next
        }
        
        # Get confidence intervals based on Wald's Z test
        confints <- confint.merMod(focal_model, method = confint_method)
        # Get coefficient estimates
        fixefs <- fixef(focal_model)
        fixefs_df <- tibble(
          Variable = names(fixefs),
          Coefficient = fixefs,
          Lower.95CI = confints[(nrow(confints)-length(fixefs)+1):nrow(confints),1],
          Upper.95CI = confints[(nrow(confints)-length(fixefs)+1):nrow(confints),2],
          Confint.method = confint_method
        ) %>%
          mutate(
            Effect.direction = ifelse(Lower.95CI > 0, rep("+", nrow(.)),
                                      ifelse(Upper.95CI < 0, rep("-", nrow(.)), rep("0", nrow(.)))
            )
          )
        # Row of the df
        stat_summ_row <- tibble(
          Taxon = str_split(sp_stage, "[.]")[[1]][1],
          Stage = str_split(sp_stage, "[.]")[[1]][2],
          Vital.rate = vr,
          Model.name = model_name,
          Has.converged = has_converged(focal_model) != -1,
          Is.singular = has_converged(focal_model) == 0,
          AIC = AIC(focal_model),
          BIC = BIC(focal_model),
          Random.variance = as.data.frame(VarCorr(focal_model))$vcov[1]
        )
        # Cross-join into single df
        stat_summ_part <- cross_join(stat_summ_row, fixefs_df)
        # Add to df
        if(is.null(stat_summ_df)) {
          stat_summ_df <- stat_summ_part
        } else {
          stat_summ_df <- bind_rows(stat_summ_df, stat_summ_part)
        }
      }
    }
  }
  
  # Return df
  stat_summ_df
}


# Function for extracting summary values from the fitted models without random effects
get_stat_summ_df_worand <- function(fit_model_list) {
  # Resulting df
  stat_summ_df <- NULL
  
  for(sp_stage in names(fit_model_list)) {
    for(vr in names(fit_model_list[[sp_stage]])) {
      for(model_name in names(fit_model_list[[sp_stage]][[vr]])) {
        cat("Processing ", sp_stage, " ", vr, " ", model_name, "\n")
        # Get model
        focal_model <- fit_model_list[[sp_stage]][[vr]][[model_name]]
        # Get confidence intervals based on Wald's Z test
        confints <- confint(focal_model)
        # Get coefficient estimates
        coefs <- coef(focal_model)
        coefs_df <- tibble(
          Variable = names(coefs),
          Coefficient = coefs,
          Lower.95CI = confints[(nrow(confints)-length(coefs)+1):nrow(confints),1],
          Upper.95CI = confints[(nrow(confints)-length(coefs)+1):nrow(confints),2]
        ) %>%
          mutate(
            Effect.direction = ifelse(Lower.95CI > 0, rep("+", nrow(.)),
              ifelse(Upper.95CI < 0, rep("-", nrow(.)), rep("0", nrow(.)))
            )
          )
        # Row of the df
        stat_summ_row <- tibble(
          Taxon = str_split(sp_stage, "[.]")[[1]][1],
          Stage = str_split(sp_stage, "[.]")[[1]][2],
          Vital.rate = vr,
          Model.name = model_name,
          AIC = AIC(focal_model),
          BIC = BIC(focal_model)
        )
        # Cross-join into single df
        stat_summ_part <- cross_join(stat_summ_row, coefs_df)
        # Add to df
        if(is.null(stat_summ_df)) {
          stat_summ_df <- stat_summ_part
        } else {
          stat_summ_df <- bind_rows(stat_summ_df, stat_summ_part)
        }
      }
    }
  }
  
  # Return df
  stat_summ_df
}
