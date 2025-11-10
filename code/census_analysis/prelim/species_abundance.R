# This file contains code that summarises the abundance of different entitities within the census dataset.

# Load tidyverse package
pacman::p_load("tidyverse")

# Load census data
census_data <- read_csv("StageDependentComp/code/census_analysis/demography_2014-2022.csv")

# ===== Visualisation of raw counts =====

# Plot bar chart of different entities
plt_spprop <- ggplot(census_data, aes(x = Year, fill = Taxon)) +
  geom_bar(color = "black")

ggsave("StageDependentComp/figures/census_analysis/entity_count_prop.png", plt_spprop, width = 8, height = 6)

# Count number of entities in each year
ann_count <- census_data %>%
  count(Year) %>%
  rename(Total = n)

# Generate table with annual counts of different entities
ann_type_count <- census_data %>%
  group_by(Year) %>%
  count(Taxon) %>%
  pivot_wider(names_from = Taxon, values_from = n)

# Generate dataframe with proportions of different entity types
ann_type_prop <- ann_type_count
ann_type_prop$Total <- ann_count$Total
ann_type_prop <- ann_type_prop %>%
  mutate_at(
    vars(-c(Total, Year)),
    ~ . / Total
  ) %>%
  group_by(Year) %>%
  select(-Total) %>%
  pivot_longer(
    -Year,
    names_to = "sp",
    values_to = "proportion"
  )

ann_type_prop

# Obtain average proportions across years
avg_type_prop <- ann_type_prop %>%
  group_by(sp) %>%
  summarise(
    avg_proportion = mean(proportion, na.rm = TRUE)
  ) %>%
  mutate(avg_perc = avg_proportion * 100) %>%
  arrange(desc(avg_proportion))

# Print summary
print(avg_type_prop, n = 30)

# Plot proportions
plt_spprop_avg <- ggplot(avg_type_prop, aes(x = avg_perc, y = reorder(sp, avg_perc))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    x = "Average percentage of species across years (%)",
    y = "Species"
  )

ggsave("StageDependentComp/figures/census_analysis/entity_prop_avg.png", plt_spprop_avg, width = 8, height = 6)

# ===== Visualisation of 'cover'; count abundance scaled by length =====

# Generate table with annual counts of different entities
ann_type_cover <- census_data %>%
  group_by(Year, Taxon) %>%
  summarise(
    Total.length = sum(Length..cm.)
  )

print(ann_type_cover, n = 200)

# Calculate average 'cover' across the years for each species
avg_type_cover <- ann_type_cover %>%
  ungroup() %>%
  group_by(Taxon) %>%
  summarise(
    Avg.total.length = mean(Total.length)
  ) %>%
  arrange(desc(Avg.total.length))

print(avg_type_cover, n = 30)



