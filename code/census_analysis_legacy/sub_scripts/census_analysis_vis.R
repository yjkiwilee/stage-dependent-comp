################################################################################
#' Script for visualising density-dependence
#'
#' Written by Young Jun Lee
#' Apr 2025
#' 

# Merge bootstrap data
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_merge_bootstrap.R"))

# ======= Read bootstrap data =======

res_df <- read_csv(
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "all_bootstrap.csv"
  )
)

res_df

file.remove(
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "all_bootstrap.csv"
  )
)

# ======== Clean up bootstrap data =======

bootstrap_data <- res_df %>%
  filter(Focal.taxon != "Elymus lanceolatus" & Neighbour.taxon != "Elymus lanceolatus") %>%
  mutate(
    Focal.vital.rate.stage = paste(Focal.vital.rate, Focal.stage, sep = ","),
    Stage.dependence = Overlap.small.slope - Overlap.large.slope
  ) %>%
  mutate(
    Focal.vital.rate.stage = factor(
      Focal.vital.rate.stage,
      levels = c(
        "Growth.future,S",
        "Growth.future,L",
        "Survival.future,S",
        "Survival.future,L",
        "Is.flowering,L",
        "X..Capitulescences,L"
      )
    )
  ) %>%
  group_by(Focal.taxon, Neighbour.taxon, Focal.vital.rate.stage) %>%
  mutate(
    # Overlap.small.slope.sd = sd(Overlap.small.slope, na.rm = TRUE),
    # Overlap.large.slope.sd = sd(Overlap.large.slope, na.rm = TRUE),
    # Stage.dependence.sd = sd(Stage.dependence, na.rm = TRUE)
    # Overlap.small.slope.sd = 1,
    # Overlap.large.slope.sd = 1,
    # Stage.dependence.sd = 1
    Overlap.small.slope.mean = mean(Overlap.small.slope, na.rm = TRUE),
    Overlap.large.slope.mean = mean(Overlap.large.slope, na.rm = TRUE),
    Stage.dependence.mean = mean(Stage.dependence, na.rm = TRUE)
  ) %>%
  # mutate(
  #   Overlap.small.slope.norm = Overlap.small.slope / Overlap.small.slope.sd,
  #   Overlap.large.slope.norm = Overlap.large.slope / Overlap.large.slope.sd,
  #   Stage.dependence.norm = Stage.dependence / Stage.dependence.sd,
  # ) %>%
  mutate(
    Focal.taxon.abb = paste0(
      str_match(
        Focal.taxon,
        "^([A-Za-z]{1})[A-Za-z]* ([A-Za-z]{3})[A-Za-z]*$"
      )[,2],
      ". ",
      str_match(
        Focal.taxon,
        "^([A-Za-z]{1})[A-Za-z]* ([A-Za-z]{3})[A-Za-z]*$"
      )[,3]
    ),
    Neighbour.taxon.abb = paste0(
      str_match(
        Neighbour.taxon,
        "^([A-Za-z]{1})[A-Za-z]* ([A-Za-z]{3})[A-Za-z]*$"
      )[,2],
      ". ",
      str_match(
        Neighbour.taxon,
        "^([A-Za-z]{1})[A-Za-z]* ([A-Za-z]{3})[A-Za-z]*$"
      )[,3]
    )
  ) %>%
  mutate(
    Focal.taxon.abb = factor(
      Focal.taxon.abb,
      levels = c("L. arg", "I. gor", "E. umb"),
      labels = c("L. arg (recipient)", "P. gor (recipient)", "E. umb (recipient)")
    ),
    Neighbour.taxon.abb = factor(
      Neighbour.taxon.abb,
      levels = c("L. arg", "I. gor", "E. umb"),
      labels = c("L. arg (actor)", "P. gor (actor)", "E. umb (actor)")
    )
  )

# ====== Run statistical tests ======

p_correction <- 1

bootstrap_summ_stat <- bootstrap_data %>%
  group_by(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage) %>%
  summarise(
    N.converged = sum(.data$Has.converged),
    
    Overlap.small.slope.median = median(.data$Overlap.small.slope, na.rm = TRUE),
    Overlap.small.slope.lower = quantile(.data$Overlap.small.slope, c(0.025), na.rm = TRUE),
    Overlap.small.slope.upper = quantile(.data$Overlap.small.slope, c(0.975), na.rm = TRUE),
    Overlap.small.slope.p = ifelse(
      median(.data$Overlap.small.slope, na.rm = TRUE) > 0,
      2 * (p_correction + sum(.data$Overlap.small.slope < 0, na.rm = TRUE)) / (p_correction + sum(!is.na(.data$Overlap.small.slope))),
      2 * (p_correction + sum(.data$Overlap.small.slope > 0, na.rm = TRUE)) / (p_correction + sum(!is.na(.data$Overlap.small.slope)))
    ),
    Overlap.small.slope.mean = unique(Overlap.small.slope.mean),
    
    Overlap.large.slope.median = median(.data$Overlap.large.slope, na.rm = TRUE),
    Overlap.large.slope.lower = quantile(.data$Overlap.large.slope, c(0.025), na.rm = TRUE),
    Overlap.large.slope.upper = quantile(.data$Overlap.large.slope, c(0.975), na.rm = TRUE),
    Overlap.large.slope.p = ifelse(
      median(.data$Overlap.large.slope, na.rm = TRUE) > 0,
      2 * (p_correction + sum(.data$Overlap.large.slope < 0, na.rm = TRUE)) / (p_correction + sum(!is.na(.data$Overlap.large.slope))),
      2 * (p_correction + sum(.data$Overlap.large.slope > 0, na.rm = TRUE)) / (p_correction + sum(!is.na(.data$Overlap.large.slope)))
    ),
    Overlap.large.slope.mean = unique(Overlap.large.slope.mean),
    
    Stage.dependence.median = median(.data$Stage.dependence, na.rm = TRUE),
    Stage.dependence.lower = quantile(.data$Stage.dependence, c(0.025), na.rm = TRUE),
    Stage.dependence.upper = quantile(.data$Stage.dependence, c(0.975), na.rm = TRUE),
    Stage.dependence.p = ifelse(
      median(.data$Stage.dependence, na.rm = TRUE) > 0,
      2 * (p_correction + sum(.data$Stage.dependence < 0, na.rm = TRUE)) / (p_correction + sum(!is.na(.data$Stage.dependence))),
      2 * (p_correction + sum(.data$Stage.dependence > 0, na.rm = TRUE)) / (p_correction + sum(!is.na(.data$Stage.dependence)))
    ),
    Stage.dependence.mean = unique(Stage.dependence.mean)
  ) 

bootstrap_summ_stat

# ========== Test significance ============
signif_level <- 0.05

bootstrap_small_signif <- bootstrap_summ_stat %>%
  filter(!is.na(Overlap.small.slope.p)) %>%
  arrange(Overlap.small.slope.p) %>%
  rowid_to_column("i") %>%
  ungroup() %>%
  mutate(
    Signif.level = signif_level / (n() - i + 1)
  ) %>%
  mutate(
    Overlap.small.is.significant = Overlap.small.slope.p <= Signif.level
  ) %>%
  select(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage,
         Overlap.small.is.significant)

bootstrap_summ_stat <- bootstrap_summ_stat %>%
  left_join(bootstrap_small_signif,
            by = join_by(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage))

bootstrap_large_signif <- bootstrap_summ_stat %>%
  filter(!is.na(Overlap.large.slope.p)) %>%
  arrange(Overlap.large.slope.p) %>%
  rowid_to_column("i") %>%
  ungroup() %>%
  mutate(
    Signif.level = signif_level / (n() - i + 1)
  ) %>%
  mutate(
    Overlap.large.is.significant = Overlap.large.slope.p <= Signif.level
  ) %>%
  select(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage,
         Overlap.large.is.significant)

bootstrap_summ_stat <- bootstrap_summ_stat %>%
  left_join(bootstrap_large_signif,
            by = join_by(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage))

bootstrap_sd_signif <- bootstrap_summ_stat %>%
  filter(!is.na(Stage.dependence.p)) %>%
  arrange(Stage.dependence.p) %>%
  rowid_to_column("i") %>%
  ungroup() %>%
  mutate(
    Signif.level = signif_level / (n() - i + 1)
  ) %>%
  mutate(
    Stage.dependence.is.significant = Stage.dependence.p <= Signif.level
  ) %>%
  select(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage,
         Stage.dependence.is.significant)

bootstrap_summ_stat <- bootstrap_summ_stat %>%
  left_join(bootstrap_sd_signif,
            by = join_by(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage))

# Write summary of bootstrap runs
write_csv(
  bootstrap_summ_stat,
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "summ_bootstrap.csv"
  )
)

small_sig_stat <- bootstrap_summ_stat %>%
  filter(Overlap.small.slope.p < 0.05) %>%
  select(
    1:4 | starts_with("Overlap.small")
  )

write_csv(
  small_sig_stat,
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "summ_small.csv"
  )
)

large_sig_stat <- bootstrap_summ_stat %>%
  filter(Overlap.large.slope.p < 0.05) %>%
  select(
    1:4 | starts_with("Overlap.large")
  )

write_csv(
  large_sig_stat,
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "summ_large.csv"
  )
)

sd_sig_stat <- bootstrap_summ_stat %>%
  filter(Stage.dependence.p < 0.05) %>%
  select(
    1:4 | starts_with("Stage.dependence")
  )

write_csv(
  sd_sig_stat,
  here(
    "StageDependentComp",
    "result_data",
    "wp2",
    "bootstrap",
    "summ_sd.csv"
  )
)

bootstrap_summ_stat

# ========== Prepare for plotting ==========

# Vertical offset for the significance asterisk
signif_offset <- 0.5
# Vertical offset for failed convergence crosss
fail_offset <- 1.5

# bootstrap_boxplot <- bootstrap_summ_ %>%
#   group_by(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage) %>%
#   summarise(
#     Overlap.small.top =
#       Overlap.small.slope.,
#     Overlap.large.top =
#       max(boxplot.stats(Overlap.large.slope.norm)$stats[5], 0),
#     Stage.dependence.top =
#       max(boxplot.stats(Stage.dependence.norm)$stats[5], 0)
#   )

bootstrap_summ_stat

bootstrap_summ_plot <- bootstrap_summ_stat %>%
  group_by(Focal.taxon.abb, Neighbour.taxon.abb, Focal.vital.rate.stage) %>%
  mutate(
    Overlap.small.slope.signif = ifelse(
      Overlap.small.is.significant,
      max(Overlap.small.slope.upper / Overlap.small.slope.mean, 0) + signif_offset,
      NA
    ),
    Overlap.small.slope.fail = ifelse(
      is.na(Overlap.small.slope.median),
      fail_offset,
      NA
    ),
    Overlap.large.slope.signif = ifelse(
      Overlap.large.is.significant,
      max(Overlap.large.slope.upper / Overlap.large.slope.mean, 0) + signif_offset,
      NA
    ),
    Overlap.large.slope.fail = ifelse(
      is.na(Overlap.large.slope.median),
      fail_offset,
      NA
    ),
    Stage.dependence.signif = ifelse(
      Stage.dependence.is.significant,
      max(Stage.dependence.upper / Stage.dependence.mean, 0) + signif_offset,
      NA
    ),
    Stage.dependence.fail = ifelse(
      is.na(Stage.dependence.median),
      fail_offset,
      NA
    )
  )

# =========== Generate plots ================

bootstrap_summ_plot

yscale_lim_small <- 5.2
small_plot <- ggplot(data = bootstrap_summ_plot,
                     mapping = aes(
                       x = Focal.vital.rate.stage,
                       y = Overlap.small.slope.median / abs(Overlap.small.slope.mean),
                       ymin = Overlap.small.slope.lower / abs(Overlap.small.slope.mean),
                       ymax = Overlap.small.slope.upper / abs(Overlap.small.slope.mean),
                       color = Focal.vital.rate.stage
                     )) +
  geom_hline(yintercept = 0) +
  geom_errorbar() +
  geom_point() +
  # coord_cartesian(
  #   ylim = c(-yscale_lim_small, yscale_lim_small)
  # ) +
  # scale_y_continuous(position = "right") +
  # scale_x_discrete(position = "top") +
  # scale_x_discrete(
  #   limits = c(
  #     "Survival.future,S",
  #     "Survival.future,L",
  #     "Growth.future,S",
  #     "Growth.future,L",
  #     "Is.flowering,L",
  #     "X..Capitulescences,L"
  #   ),
  #   labels = c(
  #     bquote(sigma[j]),
  #     bquote(sigma[a]),
  #     bquote(gamma[j]),
  #     bquote(gamma[a]),
  #     bquote(phi[0]),
  #     bquote(phi[1])
  #   )
  # ) +
  # scale_color_brewer(
  #   palette = "Dark2",
  #   guide = NULL
  # ) +
  # facet_grid(
  #   cols = vars(Focal.taxon.abb),
  #   rows = vars(Neighbour.taxon.abb)
  # ) +
  # theme_jun1() +
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.title.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   panel.grid.major.x = element_blank(),
  #   panel.grid.minor.x = element_blank()
  # ) +
  labs(
    title = "Effect of juvenile density on vital rates",
    y = "Normalised estimated effect of\njuvenile density",
    x = "Focal vital rate"
  )

if(sum(!is.na(bootstrap_summ_plot$Overlap.small.slope.signif) > 0)) {
  small_plot <- small_plot + geom_text(
    mapping = aes(x = Focal.vital.rate.stage, y = Overlap.small.slope.signif),
    color = "black",
    label = "*",
    size = 5,
    fontface = "bold"
  )
}
if(sum(!is.na(bootstrap_summ_plot$Overlap.small.slope.fail) > 0)) {
  small_plot <- small_plot + geom_text(
    mapping = aes(x = Focal.vital.rate.stage, y = Overlap.small.slope.fail),
    color = "black",
    label = "×",
    size = 5
  )
}


# ggsave(
#   here(
#     "StageDependentComp",
#     "figures",
#     "census_analysis",
#     "small_dens_plot.png"
#   ),
#   small_plot,
#   width = 7.5,
#   height = 4.5,
#   scale = 0.9
# )

yscale_lim_large <- 5
large_plot <- ggplot(data = bootstrap_summ_plot,
                     mapping = aes(
                       x = Focal.vital.rate.stage,
                       y = Overlap.large.slope.median / abs(Overlap.large.slope.mean),
                       ymin = Overlap.large.slope.lower / abs(Overlap.large.slope.mean),
                       ymax = Overlap.large.slope.upper / abs(Overlap.large.slope.mean),
                       color = Focal.vital.rate.stage
                     )) +
  geom_hline(yintercept = 0) +
  geom_errorbar() +
  geom_point() +
  # coord_cartesian(
  #   ylim = c(-yscale_lim_large, yscale_lim_large)
  # ) +
  # scale_y_continuous(position = "right") +
  # scale_x_discrete(position = "top") +
  # scale_x_discrete(
  #   limits = c(
  #     "Survival.future,S",
  #     "Survival.future,L",
  #     "Growth.future,S",
  #     "Growth.future,L",
  #     "Is.flowering,L",
  #     "X..Capitulescences,L"
  #   ),
  #   labels = c(
  #     bquote(sigma[j]),
  #     bquote(sigma[a]),
  #     bquote(gamma[j]),
  #     bquote(gamma[a]),
  #     bquote(phi[0]),
  #     bquote(phi[1])
  #   )
  # ) +
  # scale_color_brewer(
  #   palette = "Dark2",
  #   guide = NULL
  # ) +
  # facet_grid(
  #   cols = vars(Focal.taxon.abb),
  #   rows = vars(Neighbour.taxon.abb)
  # ) +
  # theme_jun1() +
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.title.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   panel.grid.major.x = element_blank(),
  #   panel.grid.minor.x = element_blank()
  # ) +
  labs(
    title = "Effect of adult density on vital rates",
    y = "Normalised estimated effect of\nadult density",
    x = "Focal vital rate"
  )

if(sum(!is.na(bootstrap_summ_plot$Overlap.large.slope.signif) > 0)) {
  large_plot <- large_plot + geom_text(
    mapping = aes(x = Focal.vital.rate.stage, y = Overlap.large.slope.signif),
    color = "black",
    label = "*",
    size = 5,
    fontface = "bold"
  )
}
if(sum(!is.na(bootstrap_summ_plot$Overlap.large.slope.fail) > 0)) {
  large_plot <- large_plot + geom_text(
    mapping = aes(x = Focal.vital.rate.stage, y = Overlap.large.slope.fail),
    color = "black",
    label = "×",
    size = 5
  )
}

# large_plot

# ggsave(
#   here(
#     "StageDependentComp",
#     "figures",
#     "census_analysis",
#     "large_dens_plot.png"
#   ),
#   large_plot,
#   width = 7.5,
#   height = 5.5
# )

yscale_lim_sd <- 5
sd_plot <- ggplot(data = bootstrap_summ_plot,
                  mapping = aes(
                    x = Focal.vital.rate.stage,
                    y = Stage.dependence.median / abs(Stage.dependence.mean),
                    ymin = Stage.dependence.lower / abs(Stage.dependence.mean),
                    ymax = Stage.dependence.upper / abs(Stage.dependence.mean),
                    color = Focal.vital.rate.stage
                  )) +
  geom_hline(yintercept = 0) +
  geom_errorbar() +
  geom_point() +
  # coord_cartesian(
  #   ylim = c(-yscale_lim_sd, yscale_lim_sd)
  # ) +
  # scale_y_continuous(position = "right") +
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.title.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   panel.grid.major.x = element_blank(),
  #   panel.grid.minor.x = element_blank()
  # ) +
  labs(
    title = "Stage-dependence in density-dependent effects",
    y = "Normalised\n(Juvenile effect - Adult effect)",
    x = "Focal vital rate"
  )

if(sum(!is.na(bootstrap_summ_plot$Stage.dependence.signif) > 0)) {
  sd_plot <- sd_plot + geom_text(
    mapping = aes(x = Focal.vital.rate.stage, y = Stage.dependence.signif),
    color = "black",
    label = "*",
    size = 5,
    fontface = "bold"
  )
}
if(sum(!is.na(bootstrap_summ_plot$Stage.dependence.fail) > 0)) {
  sd_plot <- sd_plot + geom_text(
    mapping = aes(x = Focal.vital.rate.stage, y = Stage.dependence.fail),
    color = "black",
    label = "×",
    size = 5
  )
}

# sd_plot

# ggsave(
#   here(
#     "StageDependentComp",
#     "figures",
#     "census_analysis",
#     "stage_dep_plot.png"
#   ),
#   sd_plot,
#   width = 7.5,
#   height = 5.5
# )

all_plots <- wrap_plots(
    small_plot, large_plot, sd_plot,
    design = "
    12
    3#
    ",
    heights = c(0.5, 0.5)
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  scale_x_discrete(
    limits = c(
      "Survival.future,S",
      "Survival.future,L",
      "Growth.future,S",
      "Growth.future,L",
      "Is.flowering,L",
      "X..Capitulescences,L"
    ),
    labels = c(
      bquote(sigma[j]),
      bquote(sigma[a]),
      bquote(gamma[j]),
      bquote(gamma[a]),
      bquote(phi[0]),
      bquote(phi[1])
    )
  ) &
  scale_color_brewer(
    palette = "Dark2",
    guide = NULL
  ) &
  facet_grid(
    cols = vars(Focal.taxon.abb),
    rows = vars(Neighbour.taxon.abb)
  ) &
  theme_jun1() &
  theme(
    plot.tag = element_text(face = 'bold'),
    axis.text.x = element_text(size = rel(1.3))
  )


ggsave(
  here(
    "StageDependentComp",
    "figures",
    "census_analysis",
    "summ_plot_norm.png"
  ),
  all_plots,
  width = 13,
  height = 11,
  scale = 0.8
)
