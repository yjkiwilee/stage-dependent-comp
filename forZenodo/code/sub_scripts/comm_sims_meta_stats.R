################################################################################
#' Run script for calculating meta-statistics on variables
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

# ==== Load statistical data ======

meanvar_df <- read_csv(
  file.path("result_data","meanvar.csv")
)

stat_summ_df <- read_csv(
  file.path("result_data","stat_summ.csv")
)

partial_df <- read_csv(
  file.path("result_data","partial.csv"),
  col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
)

spp_vrs_df <- read_csv(
  file.path("result_data","spp_vrs.csv")
)

mse_long_df <- read_csv(
  file.path("result_data","mse_long_summ.csv"),
  col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
)

# ===== Get mean vital rates of each species =====

dd_vr_df <- meanvar_df %>%
  filter(st_dep_level == 0) %>%
  select(dd_vr, sp1_name, sp2_name, rand_seed, sp1_dd_vr_mean, sp2_dd_vr_mean) %>%
  pivot_longer(
    c(sp1_name, sp2_name, sp1_dd_vr_mean, sp2_dd_vr_mean),
    names_to = c("species", ".value"),
    names_pattern = "^sp([0-9]+)_(.+)$"
  ) %>%
  group_by(dd_vr, name) %>%
  summarise(
    dd_vr_mean = mean(dd_vr_mean)
  )

realised_vr_df <- spp_vrs_df %>%
  left_join(dd_vr_df, by = join_by(dd_vr, sp_name == name))

for(dd_vr in unique(realised_vr_df$dd_vr)) {
  dd_vr_ids <- realised_vr_df$dd_vr == dd_vr
  realised_vr_df[[dd_vr]][dd_vr_ids] <- realised_vr_df$dd_vr_mean[dd_vr_ids]
}

realised_vr_df <- realised_vr_df %>%
  select(!dd_vr_mean)

realised_vr_df_new <- NULL

for(row_i in 1:nrow(realised_vr_df)) {
  spp_vrs_row <- realised_vr_df[row_i,]
  
  vrs <- list(
    s = c(spp_vrs_row$sj[1], spp_vrs_row$sa[1]),
    g = c(spp_vrs_row$g[1]),
    r = c(spp_vrs_row$r[1]),
    f = c(spp_vrs_row$f[1])
  )
  
  sp_mpm <- vrs_to_mpm(2, vrs)
  
  spp_vrs_row_new <- spp_vrs_row %>%
    mutate(
      m_1_1 = sp_mpm$matA[1,1],
      m_1_2 = sp_mpm$matA[1,2],
      m_2_1 = sp_mpm$matA[2,1],
      m_2_2 = sp_mpm$matA[2,2],
      gen_time = generation.time(sp_mpm$matA),
      net_rep_rate = net.reproductive.rate(sp_mpm$matA),
      life_expect = life_expect_mean(sp_mpm$matU),
      rep_val_j = reproductive.value(sp_mpm$matA)[1],
      rep_val_a = reproductive.value(sp_mpm$matA)[2],
      sens_dd_vr = d_lambda(dd_vr, vrs),
      juv_sens_dd_vr = d_prop_juv(sj, sa, g, r, f, dd_vr),
      max_real_lambda = max(
        f + sa + (sj - sa - f) * 1,
        f + sa + (sj - sa - f) * 0
      ),
      min_real_lambda = min(
        f + sa + (sj - sa - f) * 1,
        f + sa + (sj - sa - f) * 0
      ),
      diff_real_lambda = abs(
        sj - sa - f
      )
    ) %>%
    select(
      !(starts_with("sens_") & (!sens_dd_vr))
    )
  
  if(is.null(realised_vr_df_new)) {
    realised_vr_df_new <- spp_vrs_row_new
  } else {
    realised_vr_df_new <- bind_rows(realised_vr_df_new, spp_vrs_row_new)
  }
}

realised_vr_df <- realised_vr_df_new

write_csv(
  realised_vr_df,
  file.path("result_data","realised_spp_vrs.csv")
)

# ===== Minimum and maximum partial R ======

partial_df

partial_df %>%
  filter(resp_var == "lambda") %>%
  group_by(exp_var) %>%
  summarise(
    min_pr = min(pr, na.rm = TRUE),
    max_pr = max(pr, na.rm = TRUE)
  )

partial_df %>%
  filter(resp_var == "real_lambda") %>%
  group_by(exp_var) %>%
  summarise(
    min_pr = min(pr, na.rm = TRUE),
    max_pr = max(pr, na.rm = TRUE)
  )

# ===== Test for difference between partial R of different drivers =====

t_test_res_df <- partial_df %>%
  group_by(dd_vr, foc_sp_name, comp_sp_name, st_dep_level, resp_var) %>%
  group_modify(
    function(df, key) {
      ab_pr <- (df %>%
        filter(exp_var == "ab"))$pr
      struct_pr <- (df %>%
        filter(exp_var == "struct"))$pr
      test_res <- t.test(ab_pr, struct_pr, paired = TRUE)
      
      return(tibble(
        t_stat = test_res$statistic,
        p_val = test_res$p.value,
        df = test_res$parameter
      ))
    }
  )

print(t_test_res_df %>%
  filter(p_val >= 0.05), n = 65)

test <- t.test(c(1,2,3), c(1,2,3))
test$statistic
test$p.value
test$parameter

# === Does the slope of partial Rsq against stage-dependence vary with species? =====

## ==== Fit LMs =====

# Column names to use in lm summary
res_colnames <- expand.grid(c("slope"), c("neg", "pos", "all"), c("est", "p_val"))
res_colnames <- paste(res_colnames[,1], res_colnames[,2], res_colnames[,3],
                      sep = "_")

# Function to fit lm to positive and negative sections
fit_partial_lm <- function(df, key) {
  partial_mod_pos <- lm(
    pr ~ st_dep_level + rand_seed,
    df %>% filter(st_dep_level >= 0)
  )
  
  partial_mod_neg <- lm(
    pr ~ st_dep_level + rand_seed,
    df %>% filter(st_dep_level <= 0)
  )
  
  partial_mod_all <- lm(
    pr ~ st_dep_level + rand_seed,
    df
  )
  
  res_row <- c(
    get_lm_est(partial_mod_neg, c("st_dep_level")),
    get_lm_est(partial_mod_pos, c("st_dep_level")),
    get_lm_est(partial_mod_all, c("st_dep_level")),
    get_lm_pval(partial_mod_neg, c("st_dep_level")),
    get_lm_pval(partial_mod_pos, c("st_dep_level")),
    get_lm_pval(partial_mod_all, c("st_dep_level"))
  )

  names(res_row) <- res_colnames

  return(as.data.frame(as.list(res_row)))
}

# Fit lms
partial_lm_summ_df <- partial_df %>%
  mutate(
    st_dep_level = as.numeric(as.character(st_dep_level)),
    rand_seed = as.factor(rand_seed)
  ) %>%
  group_by(dd_vr, foc_sp_name, comp_sp_name, resp_var) %>%
  group_modify(fit_partial_lm)

partial_lm_summ_df

partial_lm_summ_df %>%
  filter(slope_pos_p_val < 0.05) %>%
  filter(resp_var == "lambda")

partial_lm_summ_df %>%
  filter(slope_neg_p_val < 0.05) %>%
  filter(resp_var == "lambda")

print(partial_lm_summ_df %>%
  filter(slope_all_p_val < 0.05) %>%
  filter(resp_var == "real_lambda"), n = 40)


# Save summary
write_csv(
  partial_lm_summ_df,
  file.path("result_data","partial_lm_summ.csv")
)

# ===== Maximum MAPE? =====

mse_long_df %>%
  group_by(exp_name) %>%
  summarise(
    across(
      !1:7,
      list(
        min = min,
        max = max
      )
    )
  )

# ====== Summarise interaction coefficients =====

kmat_df <- read_csv(
  file.path("result_data","inter_coeffs.csv"),
  col_types = cols(dd_vr = col_character(), st_dep_level = col_factor())
)

kmat_summ_df <- kmat_df %>%
  mutate(
    dd_vr = factor(dd_vr, levels = c("sj", "sa", "g", "r", "f")),
    sp1_name = factor(sp1_name, levels = c("F2", "F1", "M", "S1", "S2"))
  ) %>%
  select(
    dd_vr, sp1_name, st_dep_level, sp1_cj, sp1_ca, sp1_hj, sp1_ha
  ) %>%
  group_by(
    dd_vr, sp1_name, st_dep_level
  ) %>%
  summarise(
    across(
      1:4,
      mean
    )
  ) %>%
  group_by(
    dd_vr, sp1_name
  ) %>%
  group_modify(
    function(df, key) {
      return(
        tibble(
          c_base = df[df$st_dep_level == "0",]$sp1_cj,
          c_max_j = df[df$st_dep_level == "1",]$sp1_cj,
          c_max_a = df[df$st_dep_level == "-1",]$sp1_ca,
          h_base = df[df$st_dep_level == "0",]$sp1_hj,
          h_max_j = df[df$st_dep_level == "1",]$sp1_hj,
          h_max_a = df[df$st_dep_level == "-1",]$sp1_ha
        )
      )
    }
  ) %>%
  mutate(
    across(
      starts_with(c("c_", "h_")),
      function(x) {
        formatC(x, format = "e", digits = 2)
      }
    )
  )

kmat_summ_df

write_csv(
  kmat_summ_df,
  file.path("result_data","inter_coeffs_summ.csv")
)







