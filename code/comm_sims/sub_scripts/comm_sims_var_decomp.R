################################################################################
#' Run script for Q2.1: Decomposition of variance
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

# Timestep to start sampling
t_start <- 600
# Timestep to end sampling
t_end <- 2000

# DF to store results
var_res_df <- NULL

# DF to store correlation between abundance and structure
cor_df <- NULL

# Iterate through density-dependent vital rates, species combinations and stage-dependence
# dd_vrs <- c("sj", "sa", "g", "r", "f")
dd_vrs <- c("sj", "sa", "g", "r", "f")
sp_names <- c("F2", "F1", "M", "S1", "S2")
sd_level <- seq(-1, 1, 0.2)

for(ddvr in dd_vrs) {
  complete_tseries <- read_csv(
    here("StageDependentComp","result_data","wp3",
         paste0("sim_tseries_", ddvr, ".csv")),
    col_types = cols(dd_vr = col_character())
  )
  
  rand_seeds <- unique(complete_tseries$rand_seed)
  
  for(sp1_name_i in 1:length(sp_names)) {
    for(sp2_name_i in sp1_name_i:length(sp_names)) {
      sp1_n <- sp_names[sp1_name_i]
      sp2_n <- sp_names[sp2_name_i]
      cat(paste("Testing",ddvr,sp1_n,sp2_n,"\n"))
      
      start_time <- Sys.time()
      cat(paste0(format(start_time, "%c"), "\n"))
      
      for(sd_level in sd_levels) {
        for(rand_sd in rand_seeds) {
          cat(paste("\tTesting",ddvr,sp1_n,sp2_n,sd_level,rand_sd,"\n"))
          tseries_sub <- complete_tseries %>%
            filter(
              sp1_name == sp1_n &
                sp2_name == sp2_n &
                st_dep_level > (sd_level - 0.01) &
                st_dep_level < (sd_level + 0.01) &
                rand_seed == rand_sd &
                time >= t_start & time < t_end
            ) %>%
            mutate(
              sp1 = sp1_s1 + sp1_s2,
              sp2 = sp2_s1 + sp2_s2,
              sp1_struct = sp1_s1 / (sp1_s1 + sp1_s2),
              sp2_struct = sp2_s1 / (sp2_s1 + sp2_s2),
              lambda_sp1 = lambda_from_vr(
                vr_sp1_s1, vr_sp1_s2,
                vr_sp1_g1,
                vr_sp1_r1,
                vr_sp1_f1
              ),
              lambda_sp2 = lambda_from_vr(
                vr_sp2_s1, vr_sp2_s2,
                vr_sp2_g1,
                vr_sp2_r1,
                vr_sp2_f1
              )
            )
          
          # Fit full and reduced models
          mods <- list(
            ab_sp1 = lm(
              lambda_sp1 ~
                sp1 + sp2,
              tseries_sub
            ),
            ab_sp2 = lm(
              lambda_sp2 ~
                sp1 + sp2,
              tseries_sub
            ),
            struct_sp1 = lm(
              lambda_sp1 ~
                sp1_struct + sp2_struct,
              tseries_sub
            ),
            struct_sp2 = lm(
              lambda_sp2 ~
                sp1_struct + sp2_struct,
              tseries_sub
            ),
            ab_struct_sp1 = lm(
              lambda_sp1 ~
                sp1 + sp2 +
                sp1_struct + sp2_struct,
              tseries_sub
            ),
            ab_struct_sp2 = lm(
              lambda_sp2 ~
                sp1 + sp2 +
                sp1_struct + sp2_struct,
              tseries_sub
            ),
            env_sp1 = lm(
              lambda_sp1 ~
                stoch_sj_1 + stoch_sa_1 + stoch_g_1 + stoch_r_1 + stoch_f_1 +
                stoch_sj_2 + stoch_sa_2 + stoch_g_2 + stoch_r_2 + stoch_f_2,
              tseries_sub
            ),
            env_sp2 = lm(
              lambda_sp2 ~
                stoch_sj_1 + stoch_sa_1 + stoch_g_1 + stoch_r_1 + stoch_f_1 +
                stoch_sj_2 + stoch_sa_2 + stoch_g_2 + stoch_r_2 + stoch_f_2,
              tseries_sub
            )
          )
          mod_summs <- lapply(mods, function(mod) { summary(mod) })
          
          # Calculate R squared
          stats_vec <- c()
          mod_names <- c("ab_sp1", "ab_sp2", "struct_sp1", "struct_sp2",
                         "env_sp1", "env_sp2")
          for(mod_name in mod_names) {
            rsq <- mod_summs[[mod_name]]$r.squared
            regdf <- mod_summs[[mod_name]]$fstatistic[2]
            resdf <- mod_summs[[mod_name]]$fstatistic[3]
            fstat <- mod_summs[[mod_name]]$fstatistic[1]
            p_val <- pf(fstat, regdf, resdf, lower.tail = FALSE)
            
            stats_vec <- c(stats_vec, rsq, regdf, resdf, fstat, p_val)
          }
          stat_names <- c("r", "regdf", "resdf", "f", "p")
          col_names <- (expand_grid(mod_name = mod_names, stat_name = stat_names) %>%
                          mutate(
                            col_n = paste(mod_name, stat_name, sep = "_")
                          ))$col_n
          names(stats_vec) <- col_names
          stats_row_df <- data.frame(as.list(stats_vec))
          
          # Calculate partial R squared
          partial_row_df <- tibble(
            partial_r_ab_sp1 = (
              sse(mods$struct_sp1) -
                sse(mods$ab_struct_sp1)
            ) / sse(mods$struct_sp1),
            partial_r_struct_sp1 = (
              sse(mods$ab_sp1) -
                sse(mods$ab_struct_sp1)
            ) / sse(mods$ab_sp1),
            partial_r_ab_sp2 = (
              sse(mods$struct_sp2) -
                sse(mods$ab_struct_sp2)
            ) / sse(mods$struct_sp2),
            partial_r_struct_sp2 = (
              sse(mods$ab_sp2) -
                sse(mods$ab_struct_sp2)
            ) / sse(mods$ab_sp2)
          )
          
          # Calculate R squared
          var_row_df <- tibble(
            dd_vr = ddvr,
            sp1_name = sp1_n,
            sp2_name = sp2_n,
            st_dep_level = sd_level,
            rand_seed = rand_sd,
            n = nrow(tseries_sub)
          ) %>%
            bind_cols(stats_row_df) %>%
            bind_cols(partial_row_df)
          
          sp1_ab_struct_cor_test <- cor.test(tseries_sub$sp1, tseries_sub$sp1_struct)
          sp2_ab_struct_cor_test <- cor.test(tseries_sub$sp2, tseries_sub$sp2_struct)
          
          cor_row_df <- tibble(
            dd_vr = ddvr,
            sp1_name = sp1_n,
            sp2_name = sp2_n,
            st_dep_level = sd_level,
            rand_seed = rand_sd,
            sp1_ab_struct_cor = sp1_ab_struct_cor_test$estimate,
            sp1_ab_struct_cor_p = sp1_ab_struct_cor_test$p.value,
            sp2_ab_struct_cor = sp2_ab_struct_cor_test$estimate,
            sp2_ab_struct_cor_p = sp2_ab_struct_cor_test$p.value,
            n = nrow(tseries_sub),
            df = nrow(tseries_sub) - 2
          )
          
          # Merge to results DF
          if(is.null(var_res_df)) {
            var_res_df <- var_row_df
          } else {
            var_res_df <- bind_rows(var_res_df, var_row_df)
          }
          if(is.null(cor_df)) {
            cor_df <- cor_row_df
          } else {
            cor_df <- bind_rows(cor_df, cor_row_df)
          }
        }
      }
      
      end_time <- Sys.time()
      cat(paste("Done in",
                as.numeric(difftime(end_time, start_time, units = c("secs"))),
                "\n"))
    }
  }
}

# Store results in DF
write_csv(
  var_res_df,
  here("StageDependentComp","result_data","wp3","var_decomp.csv")
)
write_csv(
  cor_df,
  here("StageDependentComp","result_data","wp3","ab_struct_cor.csv")
)

# ========= Visualise results ==========

var_res_df <- read_csv(
  here("StageDependentComp","result_data","wp3","var_decomp.csv")
) %>%
  mutate(
    sp1_name = factor(sp1_name, levels = c("F2", "F1", "M", "S1", "S2")),
    sp2_name = factor(sp2_name, levels = c("F2", "F1", "M", "S1", "S2"))
  )

cor_df <- read_csv(
  here("StageDependentComp","result_data","wp3","ab_struct_cor.csv")
) %>%
  mutate(
    sp1_name = factor(sp1_name, levels = c("F2", "F1", "M", "S1", "S2")),
    sp2_name = factor(sp2_name, levels = c("F2", "F1", "M", "S1", "S2"))
  )

## ====== Exploratory tests =====

summary((var_res_df %>%
  filter(dd_vr == "sj" & sp1_name == "S2" & sp2_name == "S2"))$env_sp1_r)

test_lm <- lm(
  env_sp1_r ~ sp1_name:dd_vr,
  var_res_df
)

summary(test_lm)

min(var_res_df$ab_sp1_r, var_res_df$ab_sp2_r) # 0.0066969
max(var_res_df$ab_sp1_r, var_res_df$ab_sp2_r) # 0.1399748
min(var_res_df$struct_sp1_r, var_res_df$struct_sp2_r) # < 0.0001
max(var_res_df$struct_sp1_r, var_res_df$struct_sp2_r) # 0.1106574
min(var_res_df$env_sp1_r, var_res_df$env_sp2_r) # < 0.0001
max(var_res_df$env_sp1_r, var_res_df$env_sp2_r) # 0.1106574

min(var_res_df$ab_sp1_f, var_res_df$ab_sp2_f)
max(var_res_df$ab_sp1_p, var_res_df$ab_sp2_p)

min(var_res_df$struct_sp1_f, var_res_df$struct_sp2_f)
max(var_res_df$struct_sp1_p, var_res_df$struct_sp2_p)

min(var_res_df$env_sp1_f, var_res_df$env_sp2_f)
max(var_res_df$env_sp1_p, var_res_df$env_sp2_p)

min((var_res_df %>%
  filter(struct_sp1_p < 0.05))$struct_sp1_r)

summary(cor_df)
cor_df %>%
  filter(sp1_ab_struct_cor_p >= 0.05)
print(cor_df %>%
  filter(sp2_ab_struct_cor_p >= 0.05), n = 30)

min_struct_r_summ <- var_res_df %>%
  group_by(dd_vr, sp1_name, sp2_name) %>%
  summarise(
    sd_at_min_struct_r_sp1 = .data$st_dep_level[
      .data$struct_sp1_r == min(.data$struct_sp1_r)
    ],
    min_struct_r_sp1 = min(.data$struct_sp1_r),
    sd_at_min_struct_r_sp2 = .data$st_dep_level[
      .data$struct_sp2_r == min(.data$struct_sp2_r)
    ],
    min_struct_r_sp2 = min(.data$struct_sp2_r)
  )

print(min_struct_r_summ, n = 75)

## ====== Plot pairs =======

# List for storing the plots
var_decomp_plots = list()
var_decomp_plots_woenv = list()
partial_plots = list()
var_signif_plots = list()
cor_plots = list()
min_struct_var_plots = list()

dd_vr_labels = c(
  "sj" = "Juvenile survival",
  "sa" = "Adult survival",
  "g" = "Progression",
  "r" = "Retrogression",
  "f" = "Fertility"
)

for(ddvr in dd_vrs) {
  dd_vr_df <- var_res_df %>%
    filter(
      dd_vr == ddvr
    )
  
  het_only_df <- dd_vr_df %>%
    filter(
      sp1_name != sp2_name
    )
  
  het_org_df <- het_only_df %>%
    mutate(
      r_ab = ab_sp1_r,
      r_struct = struct_sp1_r,
      r_env = env_sp1_r,
      partial_r_ab = partial_r_ab_sp1,
      partial_r_struct = partial_r_struct_sp1
    ) %>%
    select(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, r_ab, r_struct, r_env,
           partial_r_ab, partial_r_struct)
  
  het_swap_df <- het_only_df
  het_swap_df$sp1_name <- het_only_df$sp2_name
  het_swap_df$sp2_name <- het_only_df$sp1_name
  het_swap_df <- het_swap_df %>%
    mutate(
      r_ab = ab_sp2_r,
      r_struct = struct_sp2_r,
      r_env = env_sp2_r,
      partial_r_ab = partial_r_ab_sp2,
      partial_r_struct = partial_r_struct_sp2
    ) %>%
    select(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, r_ab, r_struct, r_env,
           partial_r_ab, partial_r_struct)
  
  con_only_df <- dd_vr_df %>%
    filter(sp1_name == sp2_name) %>%
    mutate(
      r_ab = ab_sp1_r,
      r_struct = struct_sp1_r,
      r_env = env_sp1_r,
      partial_r_ab = partial_r_ab_sp1,
      partial_r_struct = partial_r_struct_sp1
    ) %>%
    select(dd_vr, sp1_name, sp2_name, st_dep_level, rand_seed, r_ab, r_struct, r_env,
           partial_r_ab, partial_r_struct)
  
  plot_df <- bind_rows(con_only_df, het_org_df, het_swap_df) %>%
    pivot_longer(
      cols = c(r_ab, r_struct, r_env,
               partial_r_ab, partial_r_struct),
      names_to = "r_type",
      values_to = "r"
    ) %>%
    mutate(
      r_type = factor(
        r_type,
        levels = c("r_ab", "r_struct", "r_env", "partial_r_ab", "partial_r_struct"),
        labels = c("ab", "struct", "env", "partial_ab", "partial_struct")
      ),
      sp1_name = factor(sp1_name, levels = c("F2", "F1", "M", "S1", "S2"),
                        labels = paste("F:", c("F2", "F1", "M", "S1", "S2"))),
      sp2_name = factor(sp2_name, levels = c("F2", "F1", "M", "S1", "S2"),
                        labels = paste("C:", c("F2", "F1", "M", "S1", "S2")))
    )
  
  decomp_df <- plot_df %>%
    filter(r_type %in% c("ab", "struct", "env"))
  var_decomp_plots[[ddvr]] <- ggplot(decomp_df, aes(x = st_dep_level, y = r, color = r_type)) +
    geom_point() +
    facet_grid(cols = vars(sp1_name), rows = vars(sp2_name)) +
    scale_color_brewer(
      palette = "Set1",
      labels = c("Abundances", "Structures", "Environment"),
      guide = guide_legend(
        ncol = 1
      )
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_jun1() +
    labs(
      x = bquote(r[Δ]),
      y = bquote(R^2~"in explaining λ of focal population"),
      color = "Focal explanatory variables",
      title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
    )
  
  decomp_df_wo_env <- decomp_df %>%
    filter(r_type != "env")
  
  var_decomp_plots_woenv[[ddvr]] <- ggplot(decomp_df_wo_env, aes(x = st_dep_level, y = r, color = r_type)) +
    geom_point() +
    facet_grid(cols = vars(sp1_name), rows = vars(sp2_name)) +
    scale_color_brewer(
      palette = "Set1",
      labels = c("Abundances", "Structures"),
      guide = guide_legend(
        ncol = 1
      )
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(n.breaks = 3) +
    theme_jun1() +
    labs(
      x = bquote(r[Δ]),
      y = bquote(R^2~"in explaining λ of focal population"),
      color = "Focal explanatory variables",
      title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
    )
  
  partial_df <- plot_df %>%
    filter(r_type %in% c("partial_ab", "partial_struct"))
  partial_plots[[ddvr]] <- ggplot(partial_df, aes(x = st_dep_level, y = r, color = r_type)) +
    geom_point() +
    facet_grid(cols = vars(sp1_name), rows = vars(sp2_name)) +
    scale_color_brewer(
      palette = "Set1",
      labels = c("Abundances", "Structures"),
      guide = guide_legend(
        ncol = 1
      )
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(n.breaks = 3) +
    theme_jun1() +
    labs(
      x = bquote(r[Δ]),
      y = bquote("Partial"~R^2~"in explaining λ of focal population"),
      color = "Focal explanatory variables",
      title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
    )
  
  plot_df_struct <- dd_vr_df %>%
    select(dd_vr, sp1_name, sp2_name, rand_seed, st_dep_level, struct_sp1_p, struct_sp2_p) %>%
    pivot_longer(
      cols = c(struct_sp1_p, struct_sp2_p),
      names_to = "foc_sp",
      values_to = "p_value"
    ) %>%
    mutate(
      foc_sp = factor(foc_sp, levels = c("struct_sp1_p", "struct_sp2_p"),
                           labels = c("sp1", "sp2")),
      sp1_name = factor(sp1_name, levels = c("F2", "F1", "M", "S1", "S2"),
                        labels = paste("1:", c("F2", "F1", "M", "S1", "S2"))),
      sp2_name = factor(sp2_name, levels = c("F2", "F1", "M", "S1", "S2"),
                        labels = paste("2:", c("F2", "F1", "M", "S1", "S2"))),
      is_significant = p_value < 0.05
    )
  
  var_signif_plots[[ddvr]] <- ggplot(plot_df_struct, aes(x = st_dep_level, y = is_significant, color = foc_sp)) +
    geom_point() +
    facet_grid(cols = vars(sp1_name), rows = vars(sp2_name)) +
    scale_color_brewer(
      palette = "Set1",
      labels = c("Population 1", "Population 2"),
      guide = guide_legend(
        ncol = 1
      )
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    theme_jun1() +
    labs(
      x = bquote(r[Δ]),
      y = "Significance at α = 0.05",
      color = "Focal population",
      title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
    )
  
  plot_cor <- cor_df %>%
    filter(dd_vr == ddvr) %>%
    select(dd_vr, sp1_name, sp2_name, rand_seed, st_dep_level, sp1_ab_struct_cor, sp2_ab_struct_cor) %>%
    pivot_longer(
      cols = c(sp1_ab_struct_cor, sp2_ab_struct_cor),
      names_to = "foc_sp",
      values_to = "cor_value"
    ) %>%
    mutate(
      foc_sp = factor(foc_sp, levels = c("sp1_ab_struct_cor", "sp2_ab_struct_cor"),
                      labels = c("sp1", "sp2")),
      sp1_name = factor(sp1_name, levels = c("F2", "F1", "M", "S1", "S2"),
                        labels = paste("1:", c("F2", "F1", "M", "S1", "S2"))),
      sp2_name = factor(sp2_name, levels = c("F2", "F1", "M", "S1", "S2"),
                        labels = paste("2:", c("F2", "F1", "M", "S1", "S2")))
    )
  
  cor_plots[[ddvr]] <- ggplot(plot_cor, aes(x = st_dep_level, y = cor_value, color = foc_sp)) +
    geom_point() +
    facet_grid(cols = vars(sp1_name), rows = vars(sp2_name)) +
    scale_color_brewer(
      palette = "Set1",
      labels = c("Population 1", "Population 2"),
      guide = guide_legend(
        ncol = 1
      )
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(n.breaks = 3) +
    theme_jun1() +
    labs(
      x = bquote(r[Δ]),
      y = "r between abundance and structure",
      color = "Focal population",
      title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
    )
  # 
  # min_struct_r_summ <- decomp_df %>%
  #   filter(r_type == "struct") %>%
  #   group_by(dd_vr, sp1_name, sp2_name) %>%
  #   summarise(
  #     sd_at_min_struct_r = .data$st_dep_level[
  #       .data$r == min(.data$r)
  #     ],
  #     min_struct_r = min(.data$r)
  #   ) %>%
  #   ungroup()
  # 
  # min_struct_var_plots[[ddvr]] <- ggplot(min_struct_r_summ, aes(
  #   x = sp1_name, y = sp2_name, fill = sd_at_min_struct_r
  # )) +
  #   geom_tile() +
  #   scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) +
  #   scale_x_discrete(labels = c(
  #     "F: F2" = "F2",
  #     "F: F1" = "F1",
  #     "F: M" = "M",
  #     "F: S1" = "S1",
  #     "F: S2" = "S2"
  #   ), position = "top") +
  #   scale_y_discrete(labels = c(
  #     "C: F2" = "F2",
  #     "C: F1" = "F1",
  #     "C: M" = "M",
  #     "C: S1" = "S1",
  #     "C: S2" = "S2"
  #   ), limits = rev) +
  #   theme_jun1() +
  #   labs(
  #     x = "Focal population",
  #     y = "Competitor population",
  #     fill = bquote(r[Δ]~"at minimum"~R^2),
  #     title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
  #   )
}
# var_decomp_plots

var_decomp_plots_woenv$f

var_decomp_plot <- wrap_plots_custom(var_decomp_plots)

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "var_decomp_all.png"
  ),
  var_decomp_plot,
  width = 9,
  height = 12
)

var_decomp_plot_woenv <- wrap_plots_custom(var_decomp_plots_woenv)

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "var_decomp_woenv.png"
  ),
  var_decomp_plot_woenv,
  width = 9,
  height = 12
)

var_signif_plot <- wrap_plots_custom(var_signif_plots)

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "var_signif.png"
  ),
  var_signif_plot,
  width = 9,
  height = 12
)


partial_plot <- wrap_plots_custom(partial_plots)

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "partial.png"
  ),
  partial_plot,
  width = 9,
  height = 12
)

cor_plot <- wrap_plots_custom(cor_plots)

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "ab_struct_cor.png"
  ),
  cor_plot,
  width = 9,
  height = 12
)

min_plot <- wrap_plots_custom(min_struct_var_plots)

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "min_r.png"
  ),
  min_plot,
  width = 8,
  height = 11
)
