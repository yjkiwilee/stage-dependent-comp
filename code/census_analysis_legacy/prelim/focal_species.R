# Script to analyse the characteristics of focal species given the dataset

# Load tidyverse package
pacman::p_load("tidyverse", "ggforce")

# Load census data
census_data <- read_csv("StageDependentComp/code/census_analysis/demography_2014-2022.csv")

count(census_data, Year)

spec(census_data)

# Constants; define global figure size
FIG_W <- 6
FIG_H <- 4

# Two focal species
SP1 <- "Eriogonum umbellatum"
SP2 <- "Arenaria congesta"

# Subset data
sub_data <- census_data %>%
  filter(
    (Taxon %in% c(SP1, SP2))
  )

# Extract individual species only
sp1_data <- sub_data %>%
  filter(Taxon == SP1)
sp2_data <- sub_data %>%
  filter(Taxon == SP2)

# Calculate median plant size
sp1_med_len <- median(sp1_data$Length..cm.)
sp2_med_len <- median(sp2_data$Length..cm.)
# Populate tibble containing median size
med_data <- tibble(
  Taxon = c(SP1, SP2),
  Length..cm. = c(sp1_med_len, sp2_med_len)
)

# ===== Size distribution plot =====

plt_size_distrib <- ggplot() +
  geom_histogram(data = sub_data, aes(x = Length..cm.), binwidth = 2, boundary = 0) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, dir = "v", scales = "free") +
  theme_bw()
plt_size_distrib
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_distrib.png", plt_size_distrib, width = FIG_W, height = FIG_H)

# ====== Size - reproductive status plot ======

plt_size_reprod <- ggplot() +
  geom_point(data = sub_data, aes(x = Length..cm., y = X..Capitulescences), size = 0.5) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  theme_bw() +
  labs(
    x = "Plant length (cm)",
    y = "Number of capitulescences",
    title = "Plant length - reproduction of focal species"
  )
plt_size_reprod
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_reprod.png", plt_size_reprod, width = FIG_W, height = FIG_H)

# ===== Size - growwth plot ====

plt_size_growth <- ggplot() +
  geom_point(data = sub_data %>% filter(Died.this.census.final == 0), aes(x = Length..cm., y = Length..cm. - Length..cm..prev), size = 0.5) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  theme_bw() +
  labs(
    x = "Plant length (cm)",
    y = "Annual growth (cm)",
    title = "Plant length - growth of focal species\n(live indiv only)"
  )
plt_size_growth
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_growth.png", plt_size_growth, width = FIG_W, height = FIG_H)

# ===== Size - survival plot ====

plt_size_surv <- ggplot() +
  geom_point(data = sub_data, aes(x = Length..cm..prev, y = 1 - Died.this.census.final), size = 1, alpha = 0.3) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  geom_smooth(data = sub_data, aes(x = Length..cm..prev, y = 1 - Died.this.census.final),
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE,
              linewidth = 0.5) +
  theme_bw() +
  labs(
    x = "Plant length (cm)",
    y = "Survival",
    title = "Plant length - survival of focal species"
  )
plt_size_surv
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_surv.png", plt_size_surv, width = FIG_W, height = FIG_H)

# ====== More plots to determine stage boundary =====

# Import function needed for boundary calculation
source("StageDependentComp/code/census_analysis/utils/stage_func.R")

# Calculate stage boundaries
sp1_thres_len <- calc_thres_len(0, sp1_data)
sp2_thres_len <- calc_thres_len(0, sp2_data)

sp1_thres_len
sp2_thres_len

# Determine proportion of individuals in juvenile stage
sp1_j_prop <- nrow(sp1_data %>% filter(Length..cm. <= sp1_thres_len)) / nrow(sp1_data)
sp2_j_prop <- nrow(sp2_data %>% filter(Length..cm. <= sp2_thres_len)) / nrow(sp2_data)

sp1_j_prop
sp2_j_prop

# Test instances of retrogression
sp1_retrog <- sp1_data %>%
  filter(Length..cm. > sp1_thres_len & Length..cm..prev <= sp1_thres_len)
sp2_retrog <- sp2_data %>%
  filter(Length..cm. > sp2_thres_len & Length..cm..prev <= sp2_thres_len)

nrow(sp1_retrog) / nrow(sp1_data)
nrow(sp2_retrog) / nrow(sp2_data)

sp1_retrog %>%
  count(Year)
sp2_retrog %>%
  count(Year)


