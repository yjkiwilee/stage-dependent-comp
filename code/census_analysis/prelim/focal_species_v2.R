# Script to analyse the characteristics of focal species given the dataset

# Load tidyverse package
pacman::p_load("tidyverse", "ggforce", "here")

# Load census data
census_data <- read_csv(here("StageDependentComp", "result_data", "demography_2014-2022_with_future.csv"))

count(census_data, Year)

spec(census_data)

# Constants; define global figure size
FIG_W <- 6
FIG_H <- 4

# Two focal species
SP1 <- "Eriogonum umbellatum"
SP2 <- "Senecio crassulus"

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
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_distrib_v2.png", plt_size_distrib, width = FIG_W, height = FIG_H)

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
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_reprod_v2.png", plt_size_reprod, width = FIG_W, height = FIG_H)

# ===== Size - growwth plot ====

plt_size_growth <- ggplot() +
  geom_point(data = sub_data %>% filter(Survival.future == 1), aes(x = Length..cm., y = Growth.future), size = 0.5) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  theme_bw() +
  labs(
    x = "Plant length (cm)",
    y = "Annual growth (cm)",
    title = "Plant length - growth of focal species\n(individuals surviving in next year only)"
  )
plt_size_growth
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_growth_v2.png", plt_size_growth, width = FIG_W, height = FIG_H)

# ===== Size - survival plot ====

plt_size_surv <- ggplot() +
  geom_point(data = sub_data, aes(x = Length..cm..prev, y = Survival.future), size = 1, alpha = 0.3) +
  geom_vline(data = med_data, aes(xintercept = Length..cm.), color = "red", linetype = "dashed") +
  facet_wrap(~ Taxon, scales = "free") +
  geom_smooth(data = sub_data, aes(x = Length..cm..prev, y = Survival.future),
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
ggsave("StageDependentComp/figures/census_analysis/focal_species/length_surv_v2.png", plt_size_surv, width = FIG_W, height = FIG_H)

# ===== Size - seedling status plot ======

# ggplot() +
#   geom_boxplot(data = sub_data, aes(group = Is.seedling, y = Length..cm.)) +
#   facet_wrap(~ Taxon, scales = "free")

# ====== More plots to determine stage boundary =====

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
  filter(Length..cm..next > sp1_thres_len & Length..cm. <= sp1_thres_len)
sp2_retrog <- sp2_data %>%
  filter(Length..cm..next > sp2_thres_len & Length..cm. <= sp2_thres_len)

nrow(sp1_retrog) / nrow(sp1_data)
nrow(sp2_retrog) / nrow(sp2_data)

sp1_retrog %>%
  count(Year)
sp2_retrog %>%
  count(Year)

# Insert stage classification into the dataset
sp1_data$Stage <- ifelse(sp1_data$Length..cm. <= sp1_thres_len, "S", "L")
sp2_data$Stage <- ifelse(sp2_data$Length..cm. <= sp2_thres_len, "S", "L")

sp1_data$Stage.future <- ifelse(sp1_data$Length..cm..next <= sp1_thres_len, "S", "L")
sp2_data$Stage.future <- ifelse(sp2_data$Length..cm..next <= sp2_thres_len, "S", "L")

count(sp1_data, Stage)

# ===== Calculate stage-dependent vital rates =====

# Growth

sp1_g <- calc_stage_prog(sp1_data, sp1_thres_len)
sp2_g <- calc_stage_prog(sp2_data, sp2_thres_len)
# Incorporate into the data table
sp1_data$Progression <- sp1_g
sp2_data$Progression <- sp2_g

# Retrogression

sp1_r <- calc_stage_retrog(sp1_data, sp1_thres_len)
sp2_r <- calc_stage_retrog(sp2_data, sp2_thres_len)
# Incorporate into the data table
sp1_data$Retrogression <- sp1_r
sp2_data$Retrogression <- sp2_r

# Validate that these vital rates are NA for last survey time point
sp1_data %>%
  group_by(Stage, Stage.future, Progression, Retrogression) %>%
  count()
sp2_data %>%
  group_by(Stage, Stage.future, Progression, Retrogression) %>%
  count()

# Save data to csv files
write_csv(sp1_data, "StageDependentComp/code/census_analysis/sp1_vr.csv")
write_csv(sp2_data, "StageDependentComp/code/census_analysis/sp2_vr.csv")

# ===== Calculate annual recruitment & estimate recruitment factor =====


sp1_rec_fact <- calc_rec_fact(sp1_data)
sp2_rec_fact <- calc_rec_fact(sp2_data)

sp2_rec_fact

ggplot(sp1_rec_fact, aes(x = Total.capitulescences, y = N.new.tags)) +
  geom_point() +
  theme_bw()



# ====== Generate test annual matrices ======

# Function to generate annual matrices based on dataset with vital rates
calc_ann_mats <- function(sp_data) {
  # Get the range of years
  yrs <- unique(sp_data$Year)
  # List in which to store the matrices
  mats <- list()
  
  # Iterate through the years to generate an annual matrix
  for(yr in yrs) {
    # Select records from the given year
    records <- sp_data %>%
      filter(Year == yr)
    
    
  }
}

# Test year to generate annual matrix

sp1_data %>% filter(Died.this.census.final == 1)




