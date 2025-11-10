################################################################################
#' Run script for generating summary plots
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

dd_vrs <- c("sj", "sa", "g", "r", "f")
sp_names <- c("F2", "F1", "M", "S1", "S2")
sd_levels <- seq(-1, 1, 0.2)
rand_seeds <- c(1:10, -1 * (1:10))

# ===== Load processed data =====

meanvar_df_ci <- read_csv(
  here("StageDependentComp","local_data","wp3","meanvar_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

partial_df_ci <- read_csv(
  here("StageDependentComp","local_data","wp3","partial_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

partial_df <- read_csv(
  here("StageDependentComp","local_data","wp3","partial.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

mse_df_ci <- read_csv(
  here("StageDependentComp","local_data","wp3","mse_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

mse_long_df_ci <- read_csv(
  here("StageDependentComp","local_data","wp3","mse_long_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

mse_long_df_ci_mcomp <- read_csv(
  here("StageDependentComp","local_data","wp3","mse_long_ci_mcomp.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names),
    dd_vr = factor(dd_vr, levels = dd_vrs)
  )

# ===== Produce plots =======

dd_vr_labels <- c(
  "sj" = "Juvenile survival",
  "sa" = "Adult survival",
  "g" = "Progression",
  "r" = "Retrogression",
  "f" = "Fertility"
)

## ====== Plot means & variances =======

var_names <- c(
  "n", "s1", "s2", "struct", "dd_vr", "dl",
  "lambda", "real_lambda", "damping_r"
)
var_labels <- c(
  "n" = "abundance",
  "s1" = "juvenile abundance",
  "s2" = "adult abundance",
  "struct" = "proportion of juveniles",
  "dd_vr" = "D-D vital rate",
  "dl" = "sensitivity",
  "lambda" = "λ",
  "real_lambda" = "N(t+1)/N(t)",
  "damping_r" = "damping ratio"
)
type_labels <- c(
  "mean" = "Mean",
  "sd" = "SD of",
  "cv" = "CV of"
)

# Pooled long term


mse_resp <- c("s1_mape", "s2_mape", "n_mape")

mse_resp_var_labels <- c(
  "s1_mape" = "small abundance",
  "s2_mape" = "large abundance",
  "n_mape" = "abundance"
)

mse_exp <- c("ab", "st_ab")

mse_exp_var_labels <- c(
  "Model A",
  "Model AS"
)

ddvr_labs <- c(
  "sj" = "σ",
  "sa" = "σ",
  "g" = "γ",
  "r" = "ρ",
  "f" = "φ"
)

ddvr_sub_labs <- c(
  "sj" = "j",
  "sa" = "a",
  "g" = " ",
  "r" = " ",
  "f" = " "
)

for(rv in mse_resp) {
  mse_sub_df <- mse_long_df_ci_mcomp %>%
    filter(resp_name == rv) %>%
    filter(exp_name == "ab") %>%
    mutate(
      dd_vr = factor(dd_vr, levels = dd_vrs)
    )
  
  mse_plot <- ggplot(mse_sub_df, aes(
    x = st_dep_level, y = avg,
    ymin = lower_ci, ymax = upper_ci
  )) +
    geom_ribbon(alpha = 0.3, linetype = 0) +
    geom_line() +
    facet_grid(cols = vars(foc_sp_name), rows = vars(dd_vr),
               scales = "free_y",
               labeller = label_bquote(
                 cols = bold("FOC. SP: "~.(as.character(foc_sp_name))),
                 rows = bold("D-D. VR: "~.(ddvr_labs[dd_vr])[.(ddvr_sub_labs[dd_vr])])
               )) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(n.breaks = 3) +
    # coord_cartesian(ylim = c(0, NA)) +
    theme_jun1() +
    theme(
      plot.title = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    ) +
    labs(
      x = "Stage dependence",
      y = paste("Mean absolute % difference\n",
                "between projections"),
      title = "Model error in predicting population abundance"
    )
  
  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0(rv, "_long_merge_presentation.png")
    ),
    mse_plot,
    width = 6.5,
    height = 5.5,
    scale = 0.6
  )
}

# Long-term, example plot

mse_resp <- c("s1_mape", "s2_mape", "n_mape")

mse_resp_var_labels <- c(
  "s1_mape" = "small abundance",
  "s2_mape" = "large abundance",
  "n_mape" = "abundance"
)

mse_exp <- c("ab", "st_ab")

mse_exp_var_labels <- c(
  "Model A",
  "Model AS"
)

mse_sub_df <- mse_long_df_ci_mcomp %>%
  filter(resp_name == "n_mape" &
           exp_name %in% c("ab", "st_ab") &
           foc_sp_name == "M" &
           dd_vr == "sa") %>%
  filter(exp_name %in% c("ab", "st_ab")) %>%
  mutate(
    dd_vr = factor(dd_vr, levels = dd_vrs)
  )

mse_sub_df

mse_sub_df <- mse_sub_df %>%
  filter(exp_name == "ab")

mse_example_plot <- ggplot(mse_sub_df, aes(
  x = st_dep_level, y = avg,
  ymin = lower_ci, ymax = upper_ci
)) +
  geom_ribbon(alpha = 0.3, linetype = 0) +
  geom_line() +
  # scale_y_continuous(n.breaks = 3) +
  # coord_cartesian(ylim = c(0, NA)) +
  theme_jun1() +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) +
  labs(
    x = "Stage dependence",
    y = paste("Mean absolute % difference\n",
              "between projected population sizes"),
    title = "Accuracy of model in predicting abundance of M\nwith density-dependent adult survival"
  )

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "mape_example_plot_presentation.png"
  ),
  mse_example_plot,
  width = 9,
  height = 8,
  scale = 0.4
)

# Long-term, example plot

mse_resp <- c("s1_mape", "s2_mape", "n_mape")

mse_resp_var_labels <- c(
  "s1_mape" = "small abundance",
  "s2_mape" = "large abundance",
  "n_mape" = "abundance"
)

mse_exp <- c("ab", "st_ab")

mse_exp_var_labels <- c(
  "Model A",
  "Model AS"
)

mse_sub_df <- mse_long_df_ci_mcomp %>%
  filter(resp_name == "n_mape" &
           exp_name %in% c("ab", "st_ab") &
           foc_sp_name == "M" &
           dd_vr == "sa") %>%
  filter(exp_name %in% c("ab", "st_ab")) %>%
  mutate(
    dd_vr = factor(dd_vr, levels = dd_vrs)
  )

mse_sub_df

mse_sub_df <- mse_sub_df %>%
  filter(exp_name == "ab")

mse_example_plot <- ggplot(mse_sub_df, aes(
  x = st_dep_level, y = avg,
  ymin = lower_ci, ymax = upper_ci
)) +
  geom_ribbon(alpha = 0.3, linetype = 0) +
  geom_line() +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(0, 100)) +
  # coord_cartesian(ylim = c(0, NA)) +
  theme_jun1() +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) +
  labs(
    x = "Stage dependence",
    y = paste("Mean absolute % difference\n",
              "between projecttions"),
    title = "Accuracy of model in predicting abundance of M\nwith density-dependent adult survival"
  )

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "mape_example_plot_presentation_noscale.png"
  ),
  mse_example_plot,
  width = 9,
  height = 8,
  scale = 0.3
)
