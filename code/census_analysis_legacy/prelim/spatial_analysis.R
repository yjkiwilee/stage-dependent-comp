# Calculate the type-specific local densities of entities around each entity in the dataset

# Load tidyverse package
pacman::p_load("tidyverse", "ggforce")

# Load census data
census_data <- read_csv(here("StageDependentComp","data","demography_2014-2022.csv"))

count(census_data, Year)

spec(census_data)

# Get list of species
sp_list <- (census_data %>%
  distinct(Taxon) %>%
  arrange(Taxon))$Taxon

# Constants; define global figure size
FIG_W <- 6
FIG_H <- 4
# Radius defining neighbourhood in cm
R_C <- 20

# Import necessary functions
source(here("StageDependentComp","code","census_analysis","func","spatial_analysis.R"))

# ===== Plot location of each plant in the plots =====

plt_yr_plot_temp <- plt_yr_plot_sqrt(filter(census_data, Year == 2016 & Plot == 46))

plt_yr_plot_temp <- plt_yr_plot_temp +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent', color=NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color=NA), #transparent legend panel
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

ggsave(here("StageDependentComp","figures","census_analysis","prelim","test_plot_transparent.png"),
       plt_yr_plot_temp, width = FIG_W, height = FIG_H, scale = 0.5, bg = "transparent")

# ===== Calculate neighbourhood status between individuals  =====

# Select particular year & plot
test_data_sub <- filter(census_data, Year == 2016 & Plot == 1)

# Calculate neighbourhood matrices
test_dist_mat <- calc_dist_mat(test_data_sub)
test_const_nei <- calc_neighbour_mat_constr(test_dist_mat, R_C)
test_dyn_nei <- calc_neighbour_mat_dynr(test_dist_mat, test_data_sub)

# ===== Test out these neighbourhood algorithms by visualising =====

# Neighbourhood line plot for constant radius
test_plot_const_nei <- plt_yr_plot(test_data_sub %>% mutate(Length..cm. = R_C))
test_plot_const_nei <- add_nei_lines(test_plot_const_nei, test_data_sub, test_const_nei)
test_plot_const_nei
ggsave("StageDependentComp/figures/census_analysis/const_nei_2016_1.png", test_plot_const_nei, width = FIG_W, height = FIG_H)

# Neighbourhood line plot for dynamic radius
test_plot_dyn_nei <- plt_yr_plot(test_data_sub)
test_plot_dyn_nei <- add_nei_lines(test_plot_dyn_nei, test_data_sub, test_dyn_nei)
test_plot_dyn_nei
ggsave("StageDependentComp/figures/census_analysis/dyn_nei_2016_1.png", test_plot_dyn_nei, width = FIG_W, height = FIG_H)

# ===== Test summarised species-specific neighbourhood status =====

summ_nei_sp(test_dyn_nei, test_data_sub)

# ===== Apply neighbourhood summary across plots & years to obtain aggregate =====

# Dataframe to store the overall neighbourhood data
nei_data_constr <- tibble(
  year = numeric(),
  plot = character(),
  sp1 = character(),
  sp2 = character(),
  n = numeric(),
  sp1_n = numeric(),
  sp2_n = numeric()
)
nei_data_dynr <- tibble(
  year = numeric(),
  plot = character(),
  sp1 = character(),
  sp2 = character(),
  n = numeric(),
  sp1_n = numeric(),
  sp2_n = numeric()
)

# Get years
yrs <- (census_data %>%
  distinct(Year))$Year

# Iterate through years & run the neighbourhood summary
for(yr in yrs) {
  print(yr)
  
  # Subset data
  yr_data <- census_data %>%
    filter(Year == yr)
  # Get plots
  plots <- (yr_data %>%
    distinct(Plot))$Plot
  
  # Iterate through plots
  for(p in plots) {
    print(p)
    
    # Subset data
    plot_data <- yr_data %>%
      filter(Plot == p)
    
    # Generate summary with constant & dynamic radiii
    summ_tibble_constr <- summ_nei_sp_const(plot_data, R_C) %>%
      mutate(year = yr, plot = p) # Add column to indicate year & plot
    summ_tibble_dynr <- summ_nei_sp_dyn(plot_data) %>%
      mutate(year = yr, plot = p) # Add column to indicate year & plot
    
    # Append to global tibble
    nei_data_constr <- bind_rows(nei_data_constr, summ_tibble_constr)
    nei_data_dynr <- bind_rows(nei_data_dynr, summ_tibble_dynr)
  }
}

# Aggregate neighbourhood data across plots
nei_constr_agg <- nei_data_constr %>%
  mutate(
    sp1 = ifelse(sp1 > sp2, sp2, sp1),
    sp2 = ifelse(sp1 > sp2, sp1, sp2) # Swap based on alphabetical order
  ) %>%
  group_by(year, sp1, sp2) %>%
  summarise(
    n = sum(n),
    .groups = "drop"
  )

nei_dynr_agg <- nei_data_dynr %>%
  mutate(
    sp1 = ifelse(sp1 > sp2, sp2, sp1),
    sp2 = ifelse(sp1 > sp2, sp1, sp2) # Swap based on alphabetical order
  ) %>%
  group_by(year, sp1, sp2) %>%
  summarise(
    n = sum(n),
    .groups = "drop"
  )

# Aggregate across years
nei_constr_agg_grand <- nei_constr_agg %>%
  group_by(sp1, sp2) %>%
  summarise(
    n = sum(n),
    .groups = "drop"
  )
nei_dynr_agg_grand <- nei_dynr_agg %>%
  group_by(sp1, sp2) %>%
  summarise(
    n = sum(n),
    .groups = "drop"
  )

# Insert 0's to missing cells
for(i in 1:length(sp_list)) {
  for(j in i:length(sp_list)) {
    # If certain species pair doesn't exist
    if(nrow(nei_constr_agg_grand %>% filter(sp1 == sp_list[i] & sp2 == sp_list[j])) == 0) {
      # Insert 0
      nei_constr_agg_grand <- bind_rows(nei_constr_agg_grand,
                                        tibble(
                                          sp1 = sp_list[i],
                                          sp2 = sp_list[j],
                                          n = 0
                                        ))
    }
    
    # Same with dynamic r
    # If certain species pair doesn't exist
    if(nrow(nei_dynr_agg_grand %>% filter(sp1 == sp_list[i] & sp2 == sp_list[j])) == 0) {
      # Insert 0
      nei_dynr_agg_grand <- bind_rows(nei_dynr_agg_grand,
                                        tibble(
                                          sp1 = sp_list[i],
                                          sp2 = sp_list[j],
                                          n = 0
                                        ))
    }
  }
}

# Plot matrix
plt_constr_nei <- ggplot(nei_constr_agg_grand, aes(x = sp1, y = sp2, fill = n)) +
  geom_tile() +
  scale_x_discrete(limits = sp_list) +
  scale_fill_gradient(high = "#000000", low = "#ffffff", limits = c(0, 15000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "Species 1",
    y = "Species 2",
    fill = "Frequency of overlap",
    title = "Frequency of neighbourhood status with constant radius"
  )
plt_constr_nei
ggsave("StageDependentComp/figures/census_analysis/constr_nei_mat.png", plt_constr_nei, width = 10, height = 8)

plt_dynr_nei <- ggplot(nei_dynr_agg_grand, aes(x = sp1, y = sp2, fill = n)) +
  geom_tile() +
  scale_x_discrete(limits = sp_list) +
  scale_fill_gradient(high = "#000000", low = "#ffffff", limits = c(0, 15000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "Species 1",
    y = "Species 2",
    fill = "Frequency of overlap",
    title = "Frequency of neighbourhood status with dynamic radius"
  )
plt_dynr_nei
ggsave("StageDependentComp/figures/census_analysis/dynr_nei_mat.png", plt_dynr_nei, width = 10, height = 8)

# ===== Calculate & plot 'normalised' frquency of overlap between species ======

# Calculate 'normalised n'
nei_data_constr_norm <- nei_data_constr %>%
  mutate(
    n = ifelse(sp1 == sp2, 2*n / (sp1_n * (sp1_n - 1)), n / (sp1_n * sp2_n))
  )
nei_data_dynr_norm <- nei_data_dynr %>%
  mutate(
    n = ifelse(sp1 == sp2, 2*n / (sp1_n * (sp1_n - 1)), n / (sp1_n * sp2_n))
  )

# Take the average normalised n across all plots & years
nei_constr_norm_avg <- nei_data_constr_norm %>%
  mutate( # Swap species to fit alphabetical order
    sp1 = ifelse(sp1 > sp2, sp2, sp1),
    sp2 = ifelse(sp1 > sp2, sp1, sp2)
  ) %>%
  group_by(sp1, sp2) %>%
  summarise(
    n = mean(n),
    .groups = "drop"
  )
nei_dynr_norm_avg <- nei_data_dynr_norm %>%
  mutate( # Swap species to fit alphabetical order
    sp1 = ifelse(sp1 > sp2, sp2, sp1),
    sp2 = ifelse(sp1 > sp2, sp1, sp2)
  ) %>%
  group_by(sp1, sp2) %>%
  summarise(
    n = mean(n),
    .groups = "drop"
  )

# Insert 0's to missing cells
for(i in 1:length(sp_list)) {
  for(j in i:length(sp_list)) {
    # If certain species pair doesn't exist
    if(nrow(nei_constr_norm_avg %>% filter(sp1 == sp_list[i] & sp2 == sp_list[j])) == 0) {
      # Insert 0
      nei_constr_norm_avg <- bind_rows(nei_constr_norm_avg,
                                        tibble(
                                          sp1 = sp_list[i],
                                          sp2 = sp_list[j],
                                          n = 0
                                        ))
    }
    
    # Same with dynamic r
    # If certain species pair doesn't exist
    if(nrow(nei_dynr_norm_avg %>% filter(sp1 == sp_list[i] & sp2 == sp_list[j])) == 0) {
      # Insert 0
      nei_dynr_norm_avg <- bind_rows(nei_dynr_norm_avg,
                                      tibble(
                                        sp1 = sp_list[i],
                                        sp2 = sp_list[j],
                                        n = 0
                                      ))
    }
  }
}

# Plot matrix
plt_constr_nei <- ggplot(nei_constr_norm_avg, aes(x = sp1, y = sp2, fill = n)) +
  geom_tile() +
  scale_x_discrete(limits = sp_list) +
  scale_fill_gradient(high = "#000000", low = "#ffffff") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "Species 1",
    y = "Species 2",
    fill = "Frequency of overlap",
    title = "Mean normalised frequency of neighbourhood status with constant radius"
  )
plt_constr_nei
ggsave("StageDependentComp/figures/census_analysis/constr_nei_mat_norm.png", plt_constr_nei, width = 10, height = 8)

plt_dynr_nei <- ggplot(nei_dynr_norm_avg, aes(x = sp1, y = sp2, fill = n)) +
  geom_tile() +
  scale_x_discrete(limits = sp_list) +
  scale_fill_gradient(high = "#000000", low = "#ffffff") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "Species 1",
    y = "Species 2",
    fill = "Frequency of overlap",
    title = "Mean normalised frequency of neighbourhood status with dynamic radius"
  )
plt_dynr_nei
ggsave("StageDependentComp/figures/census_analysis/dynr_nei_mat_norm.png", plt_dynr_nei, width = 10, height = 8)


