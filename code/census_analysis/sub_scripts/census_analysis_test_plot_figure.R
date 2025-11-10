# Script for visualising a test plot

test_plot_data <- census_data %>%
  filter(
    Year == 2016 & Plot == "20"
  )

plt_yr_plot_temp <- plt_yr_plot(test_plot_data) +
  theme_jun1() +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    axis.title = element_text(size = rel(1.1)),
    title = element_text(size = rel(1.2))
  ) +
  guides(
    color = guide_legend(nrow = 5, direction = "horizontal"),
    fill = guide_legend(nrow = 5, direction = "horizontal")
  ) +
  labs(
    x = "X coordinate (cm)",
    y = "Y coordinate (cm)",
    title = "Distribution of plants in plot ID 20 in 2016"
  )

plt_yr_plot_temp

ggsave(
  here(
    "StageDependentComp",
    "figures",
    "census_analysis",
    "plot_2016_20.png"
  ),
  plt_yr_plot_temp,
  width = 6.8,
  height = 7
)
