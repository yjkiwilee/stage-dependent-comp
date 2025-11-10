################################################################################
#' Run script for processing statistical summaries
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

# ===== Load raw statistical summaries =====

meanvar_df <- read_csv(
  here("StageDependentComp","result_data","wp3","meanvar.csv")
)

stat_summ_df <- read_csv(
  here("StageDependentComp","result_data","wp3","stat_summ.csv"),
  col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
)

mse_df <- read_csv(
  here("StageDependentComp","result_data","wp3","mse_summ.csv"),
  col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
)

mse_long_df <- read_csv(
  here("StageDependentComp","result_data","wp3","mse_long_summ.csv"),
  col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
)

# ====== Functions =======

# Make DFs longer
pivot_summ_longer <- function(df) {
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("sp") & (!c(sp1_name, sp2_name)),
      names_to = c("sp_index", ".value"),
      names_pattern = "sp([0-9])_(.+)"
    ) %>%
    mutate(
      foc_sp_name = ifelse(sp_index == 1, sp1_name, sp2_name),
      comp_sp_name = ifelse(sp_index == 1, sp2_name, sp1_name)
    ) %>%
    select(
      !c(sp1_name, sp2_name, sp_index)
    ) %>%
    relocate(
      dd_vr, foc_sp_name, comp_sp_name, st_dep_level, rand_seed, .before = everything()
    )
  
  return(df_long)
}

# Calculate 95% confidence intervals across replicates
stat_summ_calc_ci <- function(df, group_vars = c(), cols = NULL) {
  df %>%
    ungroup() %>%
    group_by(dd_vr, foc_sp_name, comp_sp_name, st_dep_level) %>%
    group_by(
      pick(all_of(group_vars)), .add = TRUE
    ) %>%
    select(!rand_seed) %>%
    summarise(
      across(
        if(is.null(cols)) { everything() } else { all_of(cols) },
        list(
          avg = function(x) { mean(x, na.rm = TRUE) },
          se = function(x) { sd(x, na.rm = TRUE) / sqrt(length(x)) },
          upper_ci = function(x) {
            mean(x, na.rm = TRUE) +
              qt(0.975, df = length(x) - 1) *
              (sd(x, na.rm = TRUE) / sqrt(length(x)))
          },
          lower_ci = function(x) {
            mean(x, na.rm = TRUE) -
              qt(0.975, df = length(x) - 1) *
              (sd(x, na.rm = TRUE) / sqrt(length(x)))
          }
        )
      )
    )
}

# Summarise across stage dependence and replicates with function
apply_across_sd_rep <- function(df, fun, group_vars = c()) {
  df %>%
    group_by(dd_vr, foc_sp_name, comp_sp_name) %>%
    group_by(
      across(any_of(group_vars)),
      .add = TRUE
    ) %>%
    select(!c(st_dep_level, rand_seed)) %>%
    summarise(
      across(
        everything(),
        fun
      )
    )
}

# ====== Process mean & sds =====

meanvar_df_long <- pivot_summ_longer(meanvar_df) %>%
  mutate(
    i = row_number()
  ) %>%
  pivot_longer(
    !c(dd_vr, foc_sp_name, comp_sp_name, st_dep_level, rand_seed, i),
    names_to = c("var_name", ".value"),
    names_pattern = "^(.+)_(mean|var)$",
  ) %>%
  mutate(
    sd = sqrt(var)
  ) %>%
  mutate(
    cv = sd / mean
  ) %>%
  select(!`var`) %>%
  pivot_wider(
    names_from = var_name,
    names_glue = "{var_name}_{.value}",
    values_from = c(`mean`, `sd`, `cv`)
  ) %>%
  select(!i)

# Calculate confidence intervals
meanvar_df_ci <- stat_summ_calc_ci(meanvar_df_long)

# Store results in df
write_csv(
  meanvar_df_ci,
  here("StageDependentComp","result_data","wp3","meanvar_ci.csv")
)

# Calculate mean
meanvar_df_avg <- apply_across_sd_rep(meanvar_df_long, mean)

write_csv(
  meanvar_df_avg,
  here("StageDependentComp","result_data","wp3","meanvar_avg.csv")
)

# ===== Calculate partial r squared =====

stat_summ_df_long <- pivot_summ_longer(stat_summ_df)

resp_vars <- c(
  "dd_vr", "lambda", "real_lambda"
)
resp_vars_str <- paste0("(", paste(resp_vars, collapse = "|"), ")")

mod_names <- c(
  "ab", "struct", "ab_struct", "full", "st_ab"
)
mod_names_str <- paste0("(", paste(mod_names, collapse = "|"), ")")

partial_df <- stat_summ_df_long %>%
  select(
    c(dd_vr, foc_sp_name, comp_sp_name, st_dep_level, rand_seed) |
      ends_with("_rsq")
  ) %>%
  pivot_longer(
    ends_with("_rsq"),
    names_to = c("resp_var", ".value"),
    names_pattern = paste0("^", resp_vars_str, "_", mod_names_str, "_rsq$")
  ) %>%
  mutate(
    partial_ab = (ab_struct - struct) / (1 - struct),
    partial_struct = (ab_struct - ab) / (1 - ab),
    partial_abraw = ab,
    partial_stabraw = st_ab
  ) %>%
  select(!any_of(mod_names)) %>%
  pivot_longer(
    starts_with("partial_"),
    names_to = "exp_var",
    names_prefix = "partial_",
    values_to = "pr"
  )

write_csv(
  partial_df,
  here("StageDependentComp","result_data","wp3","partial.csv")
)

partial_df_ci <- stat_summ_calc_ci(
  partial_df,
  c("resp_var", "exp_var"),
  c("pr")
)

# Store results in df
write_csv(
  partial_df_ci,
  here("StageDependentComp","result_data","wp3","partial_ci.csv")
)

# Calculate mean
partial_df_avg <- apply_across_sd_rep(partial_df, mean, c("resp_var", "exp_var"))

write_csv(
  partial_df_avg,
  here("StageDependentComp","result_data","wp3","partial_avg.csv")
)

# ===== Process MSE =====

mse_diff <- mse_df %>%
  pivot_wider(
    names_from = "exp_name",
    values_from = ends_with("_mse")
  ) %>%
  mutate(
    s1_mse_diff = (s1_mse_ab - s1_mse_st_ab) / s1_mse_ab * 100,
    s2_mse_diff = (s2_mse_ab - s2_mse_st_ab) / s2_mse_ab * 100,
    n_mse_diff = (n_mse_ab - n_mse_st_ab) / n_mse_ab * 100
  )

mse_df_long <- mse_diff %>%
  mutate(
    foc_sp_name = ifelse(foc_sp == "sp1", sp1_name, sp2_name),
    comp_sp_name = ifelse(foc_sp == "sp1", sp2_name, sp1_name),
  ) %>%
  select(!c(foc_sp, sp1_name, sp2_name)) %>%
  pivot_longer(
    starts_with(c("s1_", "s2_", "n_")),
    names_to = c(".value", "exp_name"),
    names_pattern = "((?:s1|s2|n)_mse)_(.+)"
  )

mse_df_ci <- stat_summ_calc_ci(
  mse_df_long,
  c("exp_name"),
  c("s1_mse", "s2_mse", "n_mse")
) %>%
  pivot_longer(
    !1:5,
    names_to = c("resp_name", ".value"),
    names_pattern = "((?:s1|s2|n)_mse)_(.+)"
  )

write_csv(
  mse_df_ci,
  here("StageDependentComp","result_data","wp3","mse_ci.csv")
)

# Long projection MSE

mse_long_df_long <- mse_long_df %>%
  mutate(
    foc_sp_name = ifelse(foc_sp == "sp1", sp1_name, sp2_name),
    comp_sp_name = ifelse(foc_sp == "sp1", sp2_name, sp1_name),
  ) %>%
  select(!c(foc_sp, sp1_name, sp2_name))

mse_long_df_ci <- stat_summ_calc_ci(
  mse_long_df_long,
  c("exp_name"),
  c("s1_mape", "s2_mape", "n_mape")
) %>%
  pivot_longer(
    !1:5,
    names_to = c("resp_name", ".value"),
    names_pattern = "((?:s1|s2|n)_mape)_(.+)"
  )

write_csv(
  mse_long_df_ci,
  here("StageDependentComp","result_data","wp3","mse_long_ci.csv")
)

# Merge competitor identities
mse_long_df_ci_mcomp <- mse_long_df_long %>%
  ungroup() %>%
  group_by(dd_vr, foc_sp_name, st_dep_level, exp_name) %>%
  select(!rand_seed) %>%
  summarise(
    across(
      all_of(c("s1_mape", "s2_mape", "n_mape")),
      list(
        avg = function(x) { mean(x, na.rm = TRUE) },
        se = function(x) { sd(x, na.rm = TRUE) / sqrt(length(x)) },
        upper_ci = function(x) {
          mean(x, na.rm = TRUE) +
            qt(0.975, df = length(x) - 1) *
            (sd(x, na.rm = TRUE) / sqrt(length(x)))
        },
        lower_ci = function(x) {
          mean(x, na.rm = TRUE) -
            qt(0.975, df = length(x) - 1) *
            (sd(x, na.rm = TRUE) / sqrt(length(x)))
        }
      )
    )
  ) %>%
  pivot_longer(
    !1:4,
    names_to = c("resp_name", ".value"),
    names_pattern = "((?:s1|s2|n)_mape)_(.+)"
  )

write_csv(
  mse_long_df_ci_mcomp,
  here("StageDependentComp","result_data","wp3","mse_long_ci_mcomp.csv")
)
