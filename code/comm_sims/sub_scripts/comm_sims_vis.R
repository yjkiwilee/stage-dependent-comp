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
  here("StageDependentComp","result_data","wp3","meanvar_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

partial_df_ci <- read_csv(
  here("StageDependentComp","result_data","wp3","partial_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

partial_df <- read_csv(
  here("StageDependentComp","result_data","wp3","partial.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

mse_df_ci <- read_csv(
  here("StageDependentComp","result_data","wp3","mse_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

mse_long_df_ci <- read_csv(
  here("StageDependentComp","result_data","wp3","mse_long_ci.csv")
) %>%
  mutate(
    foc_sp_name = factor(foc_sp_name, levels = sp_names,
                         labels = paste("C:", sp_names)),
    comp_sp_name = factor(comp_sp_name, levels = sp_names,
                          labels = paste("H:", sp_names))
  )

mse_long_df_ci_mcomp <- read_csv(
  here("StageDependentComp","result_data","wp3","mse_long_ci_mcomp.csv")
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

for(var_name in var_names) {
  for(var_type in names(type_labels)) {
    plots_list <- list()
    
    for(ddvr in dd_vrs) {
      sub_df <- meanvar_df_ci %>%
        filter(dd_vr == ddvr)
      
      sub_plt <- ggplot(sub_df, aes(x = st_dep_level,
                                    y = .data[[
                                      paste0(var_name, "_", var_type, "_avg")
                                    ]],
                                    ymin = .data[[
                                      paste0(var_name, "_", var_type, "_lower_ci")
                                    ]],
                                    ymax = .data[[
                                      paste0(var_name, "_", var_type, "_upper_ci")
                                    ]])) +
        geom_ribbon(alpha = 0.3) +
        geom_line() +
        facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
        scale_x_continuous(breaks = c(-1, 0, 1)) +
        scale_y_continuous(n.breaks = 3) +
        theme_jun1() +
        labs(
          x = bquote(r[Δ]),
          y = paste(type_labels[var_type], var_labels[var_name], "of focal pop."),
          title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
        )
      
      plots_list[[ddvr]] <- sub_plt
    }
    
    plots_wrapped <- wrap_plots_custom(plots_list)
    
    ggsave(
      here(
        "StageDependentComp",
        "figures",
        "comm_sims",
        paste0(var_name, "_", var_type, ".png")
      ),
      plots_wrapped,
      width = 9,
      height = 12
    )
  }
}

## ===== Plot statistical summaries ======

### ====== Partial rsq plots ======

partial_resp <- c("dd_vr", "lambda", "real_lambda")

resp_var_labels <- c(
  "dd_vr" = "D-D vital rate",
  "lambda" = "λ",
  "real_lambda" = "N(t+1)/N(t)"
)

# partial_exp <- c("ab", "struct")
# partial_exp <- c("struct")
partial_exp <- c("abraw", "stabraw")

exp_var_labels <- c(
  "Population abundance",
  "Stage-specific abundance"
)

for(rv in partial_resp) {
  partial_plots <- list()
  
  for(ddvr in dd_vrs) {
    partial_sub_df <- partial_df_ci %>%
      filter(resp_var == rv & dd_vr == ddvr) %>%
      filter(exp_var %in% partial_exp)
    
    # partial_plots[[ddvr]] <- ggplot(partial_sub_df, aes(
    #   x = st_dep_level, y = pr_avg,
    #   ymin = pr_lower_ci, ymax = pr_upper_ci
    # )) +
    #   geom_ribbon(alpha = 0.3, linetype = 0) +
    #   geom_line() +
    #   facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
    #   scale_x_continuous(breaks = c(-1, 0, 1)) +
    #   scale_y_continuous(n.breaks = 3) +
    #   coord_cartesian(ylim = c(0, NA)) +
    #   theme_jun1() +
    #   labs(
    #     x = bquote(r[Δ]),
    #     y = paste("Partial R² in explaining",
    #               resp_var_labels[rv],
    #               "of focal pop."),
    #     title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
    #   )

    partial_plots[[ddvr]] <- ggplot(partial_sub_df, aes(
        x = st_dep_level, y = pr_avg,
        ymin = pr_lower_ci, ymax = pr_upper_ci,
        color = exp_var, fill = exp_var
      )) +
      geom_ribbon(alpha = 0.3, linetype = 0) +
      geom_line() +
      facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
      scale_color_brewer(
        palette = "Set1",
        limits = partial_exp,
        labels = exp_var_labels,
        guide = guide_legend(
          ncol = 1
        )
      ) +
      scale_fill_brewer(
        palette = "Set1",
        limits = partial_exp,
        labels = exp_var_labels,
        guide = guide_legend(
          ncol = 1
        )
      ) +
      scale_x_continuous(breaks = c(-1, 0, 1)) +
      scale_y_continuous(n.breaks = 3) +
      # coord_cartesian(ylim = c(0, NA)) +
      theme_jun1() +
      labs(
        x = bquote(r[Δ]),
        y = paste("Partial R² in explaining",
                  resp_var_labels[rv],
                  "of focal pop."),
        color = "Focal explanatory variables",
        fill = "Focal explanatory variables",
        title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
      )
    
  #   partial_sub_df <- partial_df %>%
  #     filter(resp_var == rv & dd_vr == ddvr)
  #   
  #   partial_plots[[ddvr]] <- ggplot(partial_sub_df, aes(
  #     x = st_dep_level, y = pr,
  #     color = exp_var
  #   )) +
  #     geom_point(size = 0.3) +
  #     facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
  #     scale_color_brewer(
  #       palette = "Set1",
  #       labels = c("Abundances", "Structures"),
  #       guide = guide_legend(
  #         ncol = 1
  #       )
  #     ) +
  #     scale_x_continuous(breaks = c(-1, 0, 1)) +
  #     scale_y_continuous(n.breaks = 3) +
  #     coord_cartesian(ylim = c(0, NA)) +
  #     theme_jun1() +
  #     labs(
  #       x = bquote(r[Δ]),
  #       y = paste("Partial R² in explaining",
  #                 resp_var_labels[rv],
  #                 "of focal pop."),
  #       color = "Focal explanatory variables",
  #       fill = "Focal explanatory variables",
  #       title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
  #     )
  }

  partial_plot <- wrap_plots_custom(partial_plots)

  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0("partial_", rv, "_mod.png")
    ),
    partial_plot,
    width = 9,
    height = 12
  )
}

# ======= MSE plots ======

mse_resp <- c("s1_mse", "s2_mse", "n_mse")

mse_resp_var_labels <- c(
  "s1_mse" = "small abundance",
  "s2_mse" = "large abundance",
  "n_mse" = "abundance"
)

mse_exp <- c("ab", "st_ab")

mse_exp_var_labels <- c(
  "Population abundance",
  "Stage-specific abundance"
)

for(rv in mse_resp) {
  mse_plots <- list()
  
  for(ddvr in dd_vrs) {
    mse_sub_df <- mse_df_ci %>%
      filter(resp_name == rv & dd_vr == ddvr) %>%
      filter(exp_name %in% mse_exp)
    
    mse_plots[[ddvr]] <- ggplot(mse_sub_df, aes(
      x = st_dep_level, y = avg,
      ymin = lower_ci, ymax = upper_ci,
      color = exp_name, fill = exp_name
    )) +
      geom_ribbon(alpha = 0.3, linetype = 0) +
      geom_line() +
      facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
      scale_color_brewer(
        palette = "Set1",
        limits = mse_exp,
        labels = mse_exp_var_labels,
        guide = guide_legend(
          ncol = 1
        )
      ) +
      scale_fill_brewer(
        palette = "Set1",
        limits = mse_exp,
        labels = mse_exp_var_labels,
        guide = guide_legend(
          ncol = 1
        )
      ) +
      scale_x_continuous(breaks = c(-1, 0, 1)) +
      scale_y_continuous(n.breaks = 3) +
      # coord_cartesian(ylim = c(0, NA)) +
      theme_jun1() +
      labs(
        x = bquote(r[Δ]),
        y = paste("MSE of model in predicting",
                  mse_resp_var_labels[rv],
                  "of focal pop."),
        color = "Focal explanatory variables",
        fill = "Focal explanatory variables",
        title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
      )
  }
  
  mse_plot <- wrap_plots_custom(mse_plots)
  
  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0(rv, ".png")
    ),
    mse_plot,
    width = 9,
    height = 12
  )
}

# Difference between models
for(rv in mse_resp) {
  mse_plots <- list()
  
  for(ddvr in dd_vrs) {
    mse_sub_df <- mse_df_ci %>%
      filter(resp_name == rv & dd_vr == ddvr) %>%
      filter(exp_name == "diff")
    
    mse_plots[[ddvr]] <- ggplot(mse_sub_df, aes(
      x = st_dep_level, y = avg,
      ymin = lower_ci, ymax = upper_ci
    )) +
      geom_ribbon(alpha = 0.3, linetype = 0) +
      geom_line() +
      facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
      scale_x_continuous(breaks = c(-1, 0, 1)) +
      scale_y_continuous(n.breaks = 3) +
      # coord_cartesian(ylim = c(0, NA)) +
      theme_jun1() +
      labs(
        x = bquote(r[Δ]),
        y = paste("% reduction in MSE in predicting\n",
                  mse_resp_var_labels[rv],
                  "of focal pop."),
        title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
      )
  }
  
  mse_plot <- wrap_plots_custom(mse_plots)
  
  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0(rv, "_diff.png")
    ),
    mse_plot,
    width = 9,
    height = 12
  )
}


# Long-term

mse_resp <- c("s1_mape", "s2_mape", "n_mape")

mse_resp_var_labels <- c(
  "s1_mape" = "small abundance",
  "s2_mape" = "large abundance",
  "n_mape" = "abundance"
)

mse_exp <- c("ab", "st_ab")

mse_exp_var_labels <- c(
  "Population abundance",
  "Stage-specific abundance"
)

for(rv in mse_resp) {
  mse_plots <- list()
  
  for(ddvr in dd_vrs) {
    mse_sub_df <- mse_long_df_ci %>%
      filter(resp_name == rv & dd_vr == ddvr) %>%
      filter(exp_name %in% mse_exp)
    
    mse_plots[[ddvr]] <- ggplot(mse_sub_df, aes(
      x = st_dep_level, y = avg,
      ymin = lower_ci, ymax = upper_ci,
      color = exp_name, fill = exp_name
    )) +
      geom_ribbon(alpha = 0.3, linetype = 0) +
      geom_line() +
      facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
      scale_color_brewer(
        palette = "Set1",
        limits = mse_exp,
        labels = mse_exp_var_labels,
        guide = guide_legend(
          ncol = 1
        )
      ) +
      scale_fill_brewer(
        palette = "Set1",
        limits = mse_exp,
        labels = mse_exp_var_labels,
        guide = guide_legend(
          ncol = 1
        )
      ) +
      scale_x_continuous(breaks = c(-1, 0, 1)) +
      scale_y_continuous(n.breaks = 3) +
      # coord_cartesian(ylim = c(0, NA)) +
      theme_jun1() +
      labs(
        x = bquote(r[Δ]),
        y = paste("MAPE of model in predicting\n",
                  mse_resp_var_labels[rv],
                  "of focal sp. (%)"),
        color = "Predictive model",
        fill = "Predictive model",
        title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
      )
  }
  
  mse_plot <- wrap_plots_custom(mse_plots)
  
  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0(rv, "_long.png")
    ),
    mse_plot,
    width = 9,
    height = 12
  )
}

# Long-term, pooled across competitor species


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
    filter(exp_name %in% mse_exp) %>%
    mutate(
      dd_vr = factor(dd_vr, levels = dd_vrs)
    )
  
  mse_plot <- ggplot(mse_sub_df, aes(
    x = st_dep_level, y = avg,
    ymin = lower_ci, ymax = upper_ci,
    color = exp_name, fill = exp_name
  )) +
    geom_ribbon(alpha = 0.3, linetype = 0) +
    geom_line() +
    facet_grid(cols = vars(foc_sp_name), rows = vars(dd_vr),
               scales = "free_y",
               labeller = label_bquote(
                 cols = bold("FOC. SP: "~.(as.character(foc_sp_name))),
                 rows = bold("D-D. VR: "~.(ddvr_labs[dd_vr])[.(ddvr_sub_labs[dd_vr])])
               )) +
    scale_color_brewer(
      palette = "Set1",
      limits = mse_exp,
      labels = mse_exp_var_labels,
      guide = guide_legend(
        nrow = 1
      )
    ) +
    scale_fill_brewer(
      palette = "Set1",
      limits = mse_exp,
      labels = mse_exp_var_labels,
      guide = guide_legend(
        nrow = 1
      )
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(n.breaks = 3) +
    # coord_cartesian(ylim = c(0, NA)) +
    theme_jun1() +
    labs(
      x = bquote(r[Δ]),
      y = paste("MAPE of model in predicting\n",
                mse_resp_var_labels[rv],
                "of focal sp. (%)"),
      color = "Predictive model",
      fill = "Predictive model",
      title = "Model error in predicting population abundance"
    )
  
  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0(rv, "_long_merge.png")
    ),
    mse_plot,
    width = 6,
    height = 5.5
  )
}

# Long-term, env only

mse_resp <- c("s1_mape", "s2_mape", "n_mape")

mse_resp_var_labels <- c(
  "s1_mape" = "small abundance",
  "s2_mape" = "large abundance",
  "n_mape" = "abundance"
)

mse_exp <- c("env_only")

for(rv in mse_resp) {
  mse_plots <- list()
  
  for(ddvr in dd_vrs) {
    mse_sub_df <- mse_long_df_ci %>%
      filter(resp_name == rv & dd_vr == ddvr) %>%
      filter(exp_name %in% mse_exp)
    
    mse_plots[[ddvr]] <- ggplot(mse_sub_df, aes(
      x = st_dep_level, y = avg,
      ymin = lower_ci, ymax = upper_ci
    )) +
      geom_ribbon(alpha = 0.3, linetype = 0) +
      geom_line() +
      facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
      scale_x_continuous(breaks = c(-1, 0, 1)) +
      scale_y_continuous(n.breaks = 3) +
      # coord_cartesian(ylim = c(0, NA)) +
      theme_jun1() +
      labs(
        x = bquote(r[Δ]),
        y = paste("MAPE of model in predicting\n",
                  mse_resp_var_labels[rv],
                  "of focal sp. (%)"),
        title = paste("Density-dependent", tolower(dd_vr_labels[ddvr]))
      )
  }
  
  mse_plot <- wrap_plots_custom(mse_plots)
  
  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0(rv, "_long_envonly.png")
    ),
    mse_plot,
    width = 9,
    height = 12
  )
}

# Long-term, env only, merged


mse_resp <- c("s1_mape", "s2_mape", "n_mape")

mse_resp_var_labels <- c(
  "s1_mape" = "small abundance",
  "s2_mape" = "large abundance",
  "n_mape" = "abundance"
)

mse_exp <- c("env_only")

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
    filter(exp_name %in% mse_exp) %>%
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
    labs(
      x = bquote(r[Δ]),
      y = paste("MAPE of model in predicting\n",
                mse_resp_var_labels[rv],
                "of focal sp. (%)"),
      title = "Model error in predicting population abundance"
    )
  
  ggsave(
    here(
      "StageDependentComp",
      "figures",
      "comm_sims",
      paste0(rv, "_long_merge_envonly.png")
    ),
    mse_plot,
    width = 6,
    height = 5.5
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

mse_example_plot <- ggplot(mse_sub_df, aes(
  x = st_dep_level, y = avg,
  ymin = lower_ci, ymax = upper_ci,
  color = exp_name, fill = exp_name
)) +
  geom_ribbon(alpha = 0.3, linetype = 0) +
  geom_line() +
  # facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
  scale_color_brewer(
    palette = "Set1",
    limits = mse_exp,
    labels = mse_exp_var_labels,
    guide = guide_legend(
      ncol = 1
    )
  ) +
  scale_fill_brewer(
    palette = "Set1",
    limits = mse_exp,
    labels = mse_exp_var_labels,
    guide = guide_legend(
      ncol = 1
    )
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  # scale_y_continuous(n.breaks = 3) +
  # coord_cartesian(ylim = c(0, NA)) +
  theme_jun1() +
  labs(
    x = bquote(r[Δ]),
    y = paste("MAPE of model in predicting\n",
              "abundance of focal sp. (%)"),
    color = "Predictive model",
    fill = "Predictive model",
    title = "Accuracy of model in predicting abundance of M\nwith density-dependent adult survival"
  )

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "mape_example_plot.png"
  ),
  mse_example_plot,
  width = 7,
  height = 5.5,
  scale = 0.8
)

# Long-term, example plot, absolute y axis

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

mse_example_plot <- ggplot(mse_sub_df, aes(
  x = st_dep_level, y = avg,
  ymin = lower_ci, ymax = upper_ci,
  color = exp_name, fill = exp_name
)) +
  geom_ribbon(alpha = 0.3, linetype = 0) +
  geom_line() +
  # facet_grid(cols = vars(foc_sp_name), rows = vars(comp_sp_name)) +
  scale_color_brewer(
    palette = "Set1",
    limits = mse_exp,
    labels = mse_exp_var_labels,
    guide = guide_legend(
      ncol = 1
    )
  ) +
  scale_fill_brewer(
    palette = "Set1",
    limits = mse_exp,
    labels = mse_exp_var_labels,
    guide = guide_legend(
      ncol = 1
    )
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(0, 100)) +
  # coord_cartesian(ylim = c(0, NA)) +
  theme_jun1() +
  labs(
    x = bquote(r[Δ]),
    y = paste("MAPE of model in predicting\n",
              "abundance of focal sp. (%)"),
    color = "Predictive model",
    fill = "Predictive model",
    title = "Accuracy of model in predicting abundance of M\nwith density-dependent adult survival"
  )

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "mape_example_plot_noscale.png"
  ),
  mse_example_plot,
  width = 7,
  height = 5.5,
  scale = 0.8
)

