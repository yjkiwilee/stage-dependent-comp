################################################################################
# Script for plotting the vital rates of the focal species
#
# Written by Young Jun Lee
# Oct 2024

# ===== Dependencies =====

# Check RUN_STATE
if(!exists("RUN_STATE") || RUN_STATE == FALSE) {

  # Setup ggplot
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_setup_ggplot.R"))
  
  # Load modified census data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_load_proc_data.R"))
  
  # Choose focal species & subset data
  source(here("StageDependentComp", "code", "census_analysis", "sub_scripts", "census_analysis_set_sp.R"))
  
}

# ===== Exploratory analyses on plant vital rates =====

# # Size distribution across taxa
# ggplot(data = census_data, aes(x = Height..cm., y = Taxon)) +
#   geom_boxplot()
# 
# summary(census_data$Length..cm.)
# quantile(census_data$Length..cm., probs = c(0, 0.025, 0.95, 0.975, 1))
# 
# # Density of individuals
# dens_data <- census_data %>%
#   ungroup() %>%
#   group_by(Plot, Year) %>%
#   summarise(
#     N = n(),
#     Density = n() / 4
#   )
# 
# # Average density
# mean(dens_data[dens_data$Year == 2020,]$Density)
# 
# 
# ggplot(data = dens_data, aes(x = as.factor(Year), y = Density)) +
#   geom_boxplot()
# 
# # Proportion of different species
# prop_data <- census_data %>%
#   group_by(Taxon) %>%
#   summarise(
#     n = length(unique(Tag))
#   ) %>%
#   mutate(
#     prop = n / sum(.data$n)
#   ) %>%
#   arrange(desc(n))
# 
# print(prop_data, n = 30)

# ===== Calculate median plant size =====

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
  geom_histogram(data = sub_data %>% filter(Died.this.census.final == 0), aes(x = Length..cm.), binwidth = 2, boundary = 0) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, dir = "v", scales = "free") +
  theme_jun1() +
  labs(
     x = "Plant length (cm)",
     y = "Count",
     title = "Distribution of plant length"
  )
plt_size_distrib
ggsave(here("StageDependentComp","figures","census_analysis","focal_species","length_distrib.png"), plt_size_distrib, width = FIG_W, height = FIG_H)

# ====== Age - size plot =====

plt_age_size <- ggplot() +
  geom_point(data = sub_data %>% filter(Died.this.census.final == 0 & !is.na(Age)), aes(x = Age, y = Length..cm.), size = 0.5) +
  facet_wrap(~ Taxon, scales = "free") +
  theme_jun1() +
  labs(
    x = "Plant age (yr)",
    y = "Plant length (cm)",
    title = "Plant length - age of focal species"
  )
plt_age_size
ggsave(here("StageDependentComp","figures","census_analysis","focal_species","length_age.png"), plt_age_size, width = FIG_W, height = FIG_H)

# ====== Age - reproduction plot =====

plt_age_reprod <- ggplot() +
  geom_point(data = sub_data %>% filter(Died.this.census.final == 0 & !is.na(Age)), aes(x = Age, y = X..Capitulescences), size = 0.5) +
  facet_wrap(~ Taxon, scales = "free") +
  theme_jun1() +
  labs(
    x = "Plant age (yr)",
    y = "Number of capitulescences",
    title = "Plant age - reproduction of focal species"
  )
plt_age_reprod
ggsave(here("StageDependentComp","figures","census_analysis","focal_species","age_reprod.png"), plt_age_reprod, width = FIG_W, height = FIG_H)

# ====== Size - reproductive status plot ======

plt_size_reprod <- ggplot() +
  geom_point(data = sub_data %>% filter(Died.this.census.final == 0), aes(x = Length..cm., y = X..Capitulescences), size = 0.5) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  theme_jun1() +
  labs(
    x = "Plant length (cm)",
    y = "Number of capitulescences",
    title = "Plant length - reproduction of focal species"
  )
plt_size_reprod
ggsave(here("StageDependentComp","figures","census_analysis","focal_species","length_reprod.png"), plt_size_reprod, width = FIG_W, height = FIG_H)

# ===== Size - growwth plot ====

plt_size_growth <- ggplot() +
  geom_point(data = sub_data %>% filter(Survival.future == 1), aes(x = Length..cm., y = Growth.future), size = 0.5) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  theme_jun1() +
  labs(
    x = "Plant length (cm)",
    y = "Annual growth (cm)",
    title = "Plant length - growth of focal species\n(individuals surviving in next year only)"
  )
plt_size_growth
ggsave(here("StageDependentComp","figures","census_analysis","focal_species","length_growth.png"), plt_size_growth, width = FIG_W, height = FIG_H)

# ===== Size - survival plot ====

plt_size_surv <- ggplot() +
  geom_point(data = sub_data %>% filter(Died.this.census.final == 0), aes(x = Length..cm., y = Survival.future), size = 1, alpha = 0.3) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  geom_smooth(data = sub_data, aes(x = Length..cm., y = Survival.future),
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE,
              linewidth = 0.5) +
  theme_jun1() +
  labs(
    x = "Plant length (cm)",
    y = "Survival",
    title = "Plant length - survival of focal species"
  )
plt_size_surv
ggsave(here("StageDependentComp","figures","census_analysis","focal_species","length_surv.png"), plt_size_surv, width = FIG_W, height = FIG_H)

