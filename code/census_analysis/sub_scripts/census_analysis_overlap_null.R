
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_packages.R"))
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_stat_packages.R"))
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_proc_data.R"))

source(here("StageDependentComp","code","census_analysis","func","spatial_analysis.R"))


simulate_dist_overlap <- function(n, rf, rn, prob_unit) {
  dist_overlap <- ifelse(
    runif(n) > (prob_unit * (rf + rn) ^ 2), # If there is no overlap
    0,
    (1 - sqrt(runif(n))) / (rf + rn)^2
  )
}

simulate_dist_overlap_vec <- Vectorize(simulate_dist_overlap)

simulate_tot_dist_overlap <- function(n, n_plant, rf, rn, prob_unit) {
  overlap_mat <- simulate_dist_overlap_vec(rep(n_plant, n), rf, rn, prob_unit)
  apply(overlap_mat, 2, sum)
}

# Function for simulating plot data
# simulate_plot <- function(n_foc, n_nei, foc_length, nei_length, plot_size = 1000) {
#   margin_width <- foc_length + nei_length
#   
  sim_foc_data <- tibble(
    Year = "0",
    Plot = "0",
    Tag = as.character(seq(1, n_foc)),
    X..cm. = runif(n_foc, min = margin_width, max = plot_size - margin_width),
    Y..cm. = runif(n_foc, min = margin_width, max = plot_size - margin_width),
    Length..cm. = foc_length
  ) %>%
    mutate(
      ID = paste(Plot, Year, Tag, sep = "/")
    )

  sim_nei_data <- tibble(
    Year = "0",
    Plot = "0",
    Tag = paste0(seq(1, n_nei), "n"),
    X..cm. = runif(n_nei, min = 0, max = plot_size),
    Y..cm. = runif(n_nei, min = 0, max = plot_size),
    Length..cm. = nei_length
  )
#   
#   sim_pairs <- calc_dist_overlap_pairs(sim_foc_data, sim_nei_data)
#   sim_data$Simulated.overlap <- calc_overlap_tot_frompair(sim_pairs, sim_nei_data)
#   
#   return(sim_data)
# }

# Simulate plot data

# sim_data <- simulate_plot(10000, 10000, 20, 20)

mean(simulate_dist_overlap(100000, 30, 20, 0.0001))

mean(sqrt(runif(10000)))

(1/3)*(0.0001 * (20 + 20)^2)

ggplot(tibble(x = simulate_tot_dist_overlap(10000, 500, 20, 20, 0.0001)), aes(x = x)) +
  geom_histogram()


test_x <- rnorm(1000, 0, 100)
test_y <- rnorm(1000, 10, 100)
test_z <- 2 * test_x * test_y
test_dat <- tibble(x = test_x, y = test_y, z = test_z)

summary(lm(z ~ x * y, test_dat))




