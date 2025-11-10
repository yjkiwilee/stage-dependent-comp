census_data <- read_csv(here("StageDependentComp",
                             "result_data",
                             "data_cleanup",
                             "demography_2014-2022_processed.csv"))


# Density of plants
temp_census_data <- census_data %>%
  filter(Died.this.census.final == 0, Year == 2014)
nrow(temp_census_data) / (length(unique(temp_census_data$Plot)) * 4)

# Proportion of reproductive small individuals

n_reprod <- census_data_id %>%
  filter(Length..cm. <= Stage.boundary & Is.flowering == 1) %>%
  group_by(Taxon) %>%
  summarise(
    n = n()
  )

prop_reprod <- census_data_id %>%
  group_by(Taxon) %>%
  filter(Length..cm. <= Stage.boundary) %>%
  summarise(
    n_tot = n()
  ) %>%
  left_join(n_reprod, by = join_by(Taxon)) %>%
  mutate(
    prop = n / n_tot
  )

prop_reprod %>% filter(!is.na(n))

# Abundance of different species

temp_census_data <- census_data %>%
  filter(Died.this.census.final == 0)

census_prop <- temp_census_data %>%
  count(Taxon) %>%
  arrange(desc(n)) %>%
  mutate(
    percentage = n / nrow(temp_census_data) * 100
  )

census_prop$Taxon <- factor(
  census_prop$Taxon,
  levels = unique(census_prop$Taxon)
)

tax_prop_plt <- ggplot(census_prop, aes(x = "", y = percentage, fill = Taxon)) +
  geom_bar(stat="identity", width = 1, color = "#333", linewidth = 0.4) +
  coord_polar(theta = "y") +
  theme_jun1() +
  theme(
    legend.position = "right",
    axis.title.y = element_blank(),
    panel.border = element_blank()
  ) +
  labs(
    title = "Taxonomic distribution of live individuals\nacross censuses",
    y = "Percentage"
  )

tax_prop_plt

ggsave(
  here("StageDependentComp","figures","census_analysis","prop_taxa.png"),
  tax_prop_plt,
  width = 8, height = 6
)

# Proportion of dormancy

temp_census_data <- census_data %>%
  filter(Died.this.census.final == 0 &
           Taxon %in% c(
             "Lupinus argenteus",
             "Ivesia gordonii",
             "Eriogonum umbellatum"
           )) %>%
  mutate(
    Is.dormant = Length..cm. == 0
  ) %>%
  group_by(Taxon, Tag) %>%
  summarise(
    Is.ever.dormant = TRUE %in% Is.dormant
  ) %>%
  group_by(Taxon) %>%
  summarise(
    Percentage.dormant = sum(Is.ever.dormant) / n()
  )

temp_census_data 

# Mean and standard deviation of fecundity

temp_census_data <- census_data %>%
  filter(Died.this.census.final == 0 & Is.flowering == 1)

temp_census_data %>%
  group_by(Taxon) %>%
  summarise(
    Mean.fecundity = mean(X..Capitulescences, na.rm = TRUE),
    Var.fecundity = var(X..Capitulescences, na.rm = TRUE)
  )
  

# Individuals with age data

temp_census_data <- census_data %>%
  filter(Died.this.census.final == 0 &
           Taxon %in% c(
             "Lupinus argenteus",
             "Ivesia gordonii",
             "Eriogonum umbellatum"
           )) %>%
  group_by(Taxon, Tag) %>%
  summarise(
    Age.available = !(NA %in% Age)
  ) %>%
  group_by(Taxon) %>%
  summarise(
    Percentage.age.available = sum(Age.available) / n() * 100
  )

temp_census_data 


# Age at first flowering

temp_census_data <- census_data %>%
  filter(Died.this.census.final == 0 &
           Taxon %in% c(
             "Lupinus argenteus",
             "Ivesia gordonii",
             "Eriogonum umbellatum"
           ) &
           Is.flowering == 1) %>%
  group_by(Taxon, Tag) %>%
  summarise(
    Age.first.flowering = min(Age)
  ) %>%
  group_by(Taxon) %>%
  count(Age.first.flowering)

temp_census_data

# Vital rates against size

temp_census_data <- census_data %>%
  filter(Taxon %in% c(
    "Lupinus argenteus",
    "Ivesia gordonii",
    "Eriogonum umbellatum"
  )) %>%
  mutate(
    Growth.future = Length..cm..next - Length..cm.
  ) %>%
  select(
    Taxon, Length..cm., Survival.future, Growth.future, X..Capitulescences
  ) %>%
  pivot_longer(
    c(Survival.future, Growth.future, X..Capitulescences),
    names_to = "Vital.rate",
    values_to = "Vital.rate.val"
  )

size_vr_plt <- ggplot(temp_census_data, aes(x = Length..cm., y = Vital.rate.val)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(Vital.rate ~ Taxon, scales = "free") +
  theme_jun1() +
  labs(
    title = "Vital rates against size"
  )

size_vr_plt

ggsave(
  here("StageDependentComp","figures","census_analysis","size_vr.png"),
  size_vr_plt,
  width = 8, height = 6
)



