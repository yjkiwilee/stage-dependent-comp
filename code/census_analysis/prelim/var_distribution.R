# This file contains code that summarises the distribution of different variables of individuals.

# Load tidyverse package
pacman::p_load("tidyverse")

# Load census data
census_data <- read_csv("StageDependentComp/code/census_analysis/demography_2014-2022.csv")

spec(census_data)

FIG_W = 6
FIG_H = 4

# ===== Plot length distribution =====

plt_len_distrib <- ggplot(census_data, aes(x = Length..cm., y = reorder(Taxon, Length..cm., median, order = TRUE))) +
  geom_boxplot() +
  labs(
    x = "Length (cm)",
    y = "Taxon",
    title = "Distribution of plant length"
  ) +
  theme_bw()

ggsave("StageDependentComp/figures/census_analysis/length_distrib.png", plt_len_distrib, width = FIG_W, height = FIG_H)

# ===== Plot length distribution of particular species =====

get_length_distrib_plot <- function(sp) {
  ggplot(filter(census_data, Taxon == sp), aes(x = Length..cm.)) +
    geom_histogram(binwidth = 2.5, boundary = 0) +
    labs(
      title = paste("Plant length distribution of", sp),
      x = "Plant length (cm)",
      y = "Count"
    ) +
    theme_bw()
}

# Elymus lanceolatus

plt_elan_len_distrib <- get_length_distrib_plot("Elymus lanceolatus")

ggsave("StageDependentComp/figures/census_analysis/elan_len_distrib.png", plt_elan_len_distrib, width = FIG_W, height = FIG_H)

# Lupinus argenteus

plt_larg_len_distrib <- get_length_distrib_plot("Lupinus argenteus")

ggsave("StageDependentComp/figures/census_analysis/larg_len_distrib.png", plt_larg_len_distrib, width = FIG_W, height = FIG_H)

# Carex siccata

plt_csic_len_distrib <- get_length_distrib_plot("Carex siccata")

ggsave("StageDependentComp/figures/census_analysis/csic_len_distrib.png", plt_csic_len_distrib, width = FIG_W, height = FIG_H)

# Eriogonum umbellatum

plt_eumb_len_distrib <- get_length_distrib_plot("Eriogonum umbellatum")

ggsave("StageDependentComp/figures/census_analysis/eumb_len_distrib.png", plt_eumb_len_distrib, width = FIG_W, height = FIG_H)

