
pacman::p_load("here")

source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_packages.R"))
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_stat_packages.R"))
source(here("StageDependentComp","code","census_analysis","sub_scripts","census_analysis_load_proc_data.R"))

source(here("StageDependentComp","code","census_analysis","func","spatial_analysis.R"))
source(here("StageDependentComp","code","census_analysis","func","stage_func.R"))

# Create census dataset with unique ID
census_data_id <- census_data %>%
  mutate(
    ID = paste(Plot, Year, Tag, sep = "/")
  )

# Calculate stage boundaries
census_stage_bounds <- census_data %>%
  group_by(Taxon) %>%
  filter(Is.seedling == 1) %>%
  summarise(
    Stage.boundary = max(Length..cm., na.rm = TRUE)
  )

# Merge into dataset
census_data_id <- left_join(census_data_id, census_stage_bounds,
                            by = join_by(Taxon)) %>%
  mutate( # Insert column for life stage
    Stage = ifelse(Length..cm. <= Stage.boundary, "S", "L"),
    Stage.future = ifelse(Length..cm..next <= Stage.boundary, "S", "L")
  ) %>%
  mutate(
    Progression = ifelse(Stage == "S" & Stage.future == "L", 1, 0),
    Retrogression = ifelse(Stage == "L" & Stage.future == "S", 1, 0)
  )

# # Type of vital rate
# FOCAL_VR <- "Survival.future"
# FOCAL_VR_ABB <- "fec"
# 
# # Focal species
# SP_FOCAL <- "Lupinus argenteus"
# SP_NEI <- "Lupinus argenteus"
# SP_FOCAL_ABB <- "larg"
# SP_NEI_ABB <- "larg"
# 
# BOUNDARY_WIDTH <- 0
# PLOT_DIM <- 200


focal_data_unalt <- census_data_id %>%
  filter(
    Taxon == SP_FOCAL &
      !is.na(.data[[FOCAL_VR]]) &
      !is.na(Length..cm.) &
      !is.na(X..cm.) &
      !is.na(Y..cm.) &
      X..cm. >= BOUNDARY_WIDTH & X..cm. <= PLOT_DIM - BOUNDARY_WIDTH &
      Y..cm. >= BOUNDARY_WIDTH & Y..cm. <= PLOT_DIM - BOUNDARY_WIDTH
  )

focal_data <- focal_data_unalt

if(FOCAL_VR == "X..Capitulescences") {
  focal_data <- focal_data %>%
    filter( # Only flowering live individuals
      Is.flowering == 1 &
        Survival == 1
    )
} else if(FOCAL_VR %in% c("Growth.future")) {
  focal_data <- focal_data %>%
    filter(
      Survival.future == 1
    )
} else if(FOCAL_VR %in% c("Growth")) {
  focal_data <- focal_data %>%
    filter(
      Survival == 1
    )
} else if(FOCAL_VR == "Is.flowering") {
  focal_data <- focal_data %>%
    filter( # Exclude seedlings and dead individuals
      Is.seedling == 0 &
        Survival == 1
    )
}

nei_data <- census_data_id %>%
  filter(
    Taxon == SP_NEI &
      !is.na(Length..cm.) &
      Length..cm. > 0 &
      !is.na(X..cm.) &
      !is.na(Y..cm.)
  )

# Stage boundaries as maximum size of seedlings
# nei_stage_bound <- min(nei_data$Length..cm.[nei_data$X..Capitulescences > 0])
# nei_stage_bound <- median(nei_data$Length..cm., na.rm = TRUE)
nei_stage_bound <- max(nei_data$Length..cm.[nei_data$Is.seedling == 1], na.rm = TRUE)

nei_st <- list(
  s = nei_data %>% filter(Length..cm. <= nei_stage_bound),
  l = nei_data %>% filter(Length..cm. > nei_stage_bound)
)

other_data <- census_data_id %>%
  filter(
    Taxon != SP_NEI &
      !is.na(Length..cm.) &
      Length..cm. > 0 &
      !is.na(X..cm.) &
      !is.na(Y..cm.)
  )

# Get pairwise overlap & total overlap effects

nei_pairs <- lapply(nei_st, function(nei_dat) {
  calc_dens_overlap_pairs(focal_data, nei_dat)
})
other_pairs <- calc_dens_overlap_pairs(focal_data, other_data)

focal_data_stat <- focal_data

# Determine stage boundary
focal_stage_bound <- max(focal_data_unalt$Length..cm.[focal_data_unalt$Is.seedling == 1],
                         na.rm = TRUE)
focal_stage_bound

focal_data_stat <- focal_data_stat %>%
  mutate(
    Overlap.small = calc_overlap_tot_frompair(nei_pairs$s, focal_data),
    Overlap.large = calc_overlap_tot_frompair(nei_pairs$l, focal_data),
    Overlap.other = calc_overlap_tot_frompair(other_pairs, focal_data)
  ) %>%
  mutate(
    Overlap.focal = Overlap.small + Overlap.large,
    Overlap = Overlap.small + Overlap.large + Overlap.other,
    Stage = ifelse(Length..cm. <= focal_stage_bound, "S", "L")
  )

# Normalise length & overlap
norm_x <- function(data, var, x) {
  (x - mean(data[[var]])) / sd(data[[var]])
}
recover_x <- function(data, var, normx) {
  normx * sd(data[[var]]) + mean(data[[var]])
}

norm_vars <- c("Length..cm.",
               "Overlap.focal", "Overlap.small", "Overlap.large",
               "Overlap.other", "Overlap")

for(norm_var in norm_vars) {
  focal_data_stat[[paste0(norm_var, ".norm")]] <- norm_x(
    focal_data_stat, 
    norm_var, 
    focal_data_stat[[norm_var]]
  )
}

# Model equation

# Does the identity of the actor matter?
model_eq <- paste0(FOCAL_VR, " ~ Length..cm..norm +
  Overlap.small.norm +
  Overlap.large.norm +
  Overlap.other.norm + (1 | Year) + (1 | Plot)")

# Does the identity of the focal individual matter?
# model_eq <- paste0(FOCAL_VR, " ~ Length..cm..norm +
#   Overlap.focal.norm:Stage +
#   Overlap.other.norm + (1 | Year) + (1 | Plot)")

# model_eq <- paste0(FOCAL_VR, " ~ Length..cm..norm +
#   Overlap.small + Overlap.large +
#   Overlap.other.norm + Year")

# model_eq <- paste0(FOCAL_VR, " ~ Length..cm. +
#   Prop.overlap.focal.small + Prop.overlap.focal.small:Length..cm. +
#   Prop.overlap.focal.large + Prop.overlap.focal.large:Length..cm. +
#   Prop.overlap.other + Prop.overlap.other:Length..cm. +
#   (1 | Plot) + (1 | Year)")

# Reformulate dummy equation
model_eq <- as.formula(model_eq)

# ====== Frequentist test ======

# stat_model <- lmer(
#   model_eq,
#   data = focal_data
#   # control = lmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 7e4))
# )

# stat_model <- if(FOCAL_VR %in% c("Growth.future", "X..Capitulescences")) {
#   lmer(
#     model_eq,
#     data = focal_data
#   )
# } else if(FOCAL_VR %in% c("Survival.future", "Is.flowering")) {
#   glmer(
#     model_eq,
#     family = "binomial",
#     data = focal_data
#   )
# }

stat_model <- if(FOCAL_VR %in% c("Growth.future", "Growth", "X..Capitulescences")) {
  lmer(
    model_eq,
    data = focal_data_stat
  )
} else if(FOCAL_VR %in% c("Survival.future", "Survival", "Is.flowering")) {
  glmer(
    model_eq,
    family = "binomial",
    data = focal_data_stat,
    control = glmerControl(optimizer = "bobyqa")
  )
}

stat_model_all <- allFit(stat_model)

summary(stat_model_all)

linearHypothesis(stat_model, "0.029*Overlap.small.norm - 0.061*Overlap.large.norm = 0")
# linearHypothesis(stat_model,
#                  "Overlap.focal.norm:StageL - Overlap.focal.norm:StageS = 0")

car::Anova(stat_model, type=3)

summary(stat_model)

confint(stat_model)

ggplot(focal_data_stat, aes(x = Length..cm., y = Overlap.large,
                    color = .data[[FOCAL_VR]])) +
  geom_point(size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_jun1()

ggplot(focal_data_stat, aes(x = Overlap.large, y = .data[[FOCAL_VR]])) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_jun1()

ggplot(focal_data_stat, aes(x = Overlap.small, y = .data[[FOCAL_VR]])) +
  geom_point() +
  geom_smooth(method = "gam") +
  theme_jun1()

ggplot(focal_data_stat, aes(x = Length..cm., y = Overlap.small,
                            color = .data[[FOCAL_VR]])) +
  geom_point(size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_jun1()

ggplot(focal_data_stat, aes(x = Length..cm., y = Overlap.other,
                            color = .data[[FOCAL_VR]])) +
  geom_point(size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_jun1()

ggplot(focal_data_stat, aes(x = Overlap.large)) +
  geom_histogram()

# ====== Test model prediction ======

# Calculate model predictions
focal_data_pred <- focal_data
focal_data_pred$Vital.rate.predict <- predict(stat_model, focal_data)

# Plot estimates against real values
ggplot(focal_data_pred, aes(x = .data[[FOCAL_VR]], y = Vital.rate.predict)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black")

# ========= Plotting =========

# Get quantiles of variables
quants <- c(0, 0.05, 0.1, 0.5, 0.9, 0.95, 1)
quant_length <- quantile(focal_data$Length..cm., quants)
quant_os <- quantile(focal_data$Prop.overlap.focal.small, quants)
quant_ol <- quantile(focal_data$Prop.overlap.focal.large, quants)
quant_oo <- quantile(focal_data$Prop.overlap.other, quants)
quant_vr <- quantile(focal_data[[FOCAL_VR]], quants)

# Gradient of model output to generate heatmap
nstep <- 50
nstep_l <- nstep
nstep_os <- nstep
nstep_ol <- nstep

# Generate sequences
l_seq <- seq(quant_length[2], quant_length[6], length.out = nstep_l)
os_seq <- seq(quant_os[2], quant_os[6], length.out = nstep_os)
ol_seq <- seq(quant_ol[2], quant_ol[6], length.out = nstep_ol)


## ====== Plot modified data plots ======

# Function for calculating data points while only keeping the effect of some variables
only_effect_from <- function(in_data, fit_model, resp_var, eff_vars) {
  # Get raw response
  raw_resp <- in_data[[resp_var]]
  
  # Calculate effect from each of the variables
  var_effects <- list(
    "Plot" = coef(fit_model)$Plot[,"(Intercept)",drop = FALSE][in_data$Plot,],
    "Year" = coef(fit_model)$Year[,"(Intercept)",drop = FALSE][as.character(in_data$Year),],
    "(Intercept)" = rep(fixef(fit_model)["(Intercept)"], nrow(in_data)),
    "Length..norm" = fixef(fit_model)["Length..norm"] * in_data[["Length..norm"]],
    "Prop.overlap.focal.small.norm" = fixef(fit_model)["Prop.overlap.focal.small.norm"] *
      in_data[["Prop.overlap.focal.small.norm"]],
    "Prop.overlap.focal.large.norm" = fixef(fit_model)["Prop.overlap.focal.large.norm"] *
      in_data[["Prop.overlap.focal.large.norm"]],
    "Prop.overlap.other.norm" = fixef(fit_model)["Prop.overlap.other.norm"] *
      in_data[["Prop.overlap.other.norm"]],
    "Length..norm:Prop.overlap.focal.small.norm" = fixef(fit_model)["Length..norm:Prop.overlap.focal.small.norm"] *
      in_data[["Length..norm"]] * in_data[["Prop.overlap.focal.small.norm"]],
    "Length..norm:Prop.overlap.focal.large.norm" = fixef(fit_model)["Length..norm:Prop.overlap.focal.large.norm"] *
      in_data[["Length..norm"]] * in_data[["Prop.overlap.focal.large.norm"]],
    "Length..norm:Prop.overlap.other.norm" = fixef(fit_model)["Length..norm:Prop.overlap.other.norm"] *
      in_data[["Length..norm"]] * in_data[["Prop.overlap.other.norm"]]
  )
  
  # Remove desired variables
  var_effects <- var_effects[names(var_effects) %in% eff_vars == FALSE]
  
  # Calculate processed response
  proc_resp <- raw_resp - reduce(var_effects, `+`)
  names(proc_resp) <- NULL
  
  # Return processed response
  return(proc_resp)
}

summary(stat_model)

### ====== Small test ======

test_effect <- only_effect_from(focal_data, stat_model, FOCAL_VR,
                                c(
                                  "Length..norm",
                                  "Length..norm:Prop.overlap.focal.small.norm"
                                ))

test_df <- tibble(
  Length..cm. = focal_data$Length..cm.,
  Prop.overlap.focal.small = focal_data$Prop.overlap.focal.small,
  Vital.rate = test_effect
)

ggplot(test_df, aes(x = Length..cm., y = Prop.overlap.focal.small,
                    color = Vital.rate)) +
  geom_point(size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_jun1()

### ====== Large test ======

test_effect <- only_effect_from(focal_data, stat_model, FOCAL_VR,
                                c(
                                  "Length..norm",
                                  "Length..norm:Prop.overlap.focal.large.norm"
                                ))

test_df <- tibble(
  Length..cm. = focal_data$Length..cm.,
  Prop.overlap.focal.large = focal_data$Prop.overlap.focal.large,
  Vital.rate = test_effect
)

ggplot(test_df, aes(x = Length..cm., y = Prop.overlap.focal.large,
                    color = Vital.rate)) +
  geom_point(size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_jun1()

### ====== Other test ======

test_effect <- only_effect_from(focal_data, stat_model, FOCAL_VR,
                                c(
                                  "Length..norm",
                                  "Length..norm:Prop.overlap.other.norm"
                                ))

test_df <- tibble(
  Length..cm. = focal_data$Length..cm.,
  Prop.overlap.other = focal_data$Prop.overlap.other,
  Vital.rate = test_effect
)

ggplot(test_df, aes(x = Length..cm., y = Prop.overlap.other,
                    color = Vital.rate)) +
  geom_point(size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_jun1()


## ====== Plot expected overlap values at varying population-level density ======



# 
# ## ====== Plot gradient of vital rate against overlap =======
# 
# focal_small_gradients <- tibble(
#   Prop.overlap.focal.large = quant_ol[3],
#   Prop.overlap.focal.large.norm = recover_x(focal_data, "Prop.overlap.focal.large", quant_ol[3]),
#   Length..cm. = l_seq,
#   Length..norm = recover_x(focal_data, "Length..cm.", l_seq)
# ) %>%
#   mutate(
#     Small.gradient = ((fixef(stat_model)["Prop.overlap.focal.small.norm"] +
#       fixef(stat_model)["Length..norm:Prop.overlap.focal.small.norm"] * Length..norm) / 
#       sd(focal_data$Length..cm.)) / (pi * Length..cm.^2)
#   )
# 
# f_s_grad_plt <- ggplot(focal_small_gradients, aes(x = Length..cm., y = Small.gradient)) +
#   geom_line() +
#   theme_jun1()
# 
# f_s_grad_plt
# 
# focal_large_gradients <- tibble(
#   Prop.overlap.focal.small = quant_os[3],
#   Prop.overlap.focal.small.norm = recover_x(focal_data, "Prop.overlap.focal.small", quant_os[3]),
#   Length..cm. = l_seq,
#   Length..norm = recover_x(focal_data, "Length..cm.", l_seq)
# ) %>%
#   mutate(
#     Large.gradient = (fixef(stat_model)["Prop.overlap.focal.large.norm"] +
#                         fixef(stat_model)["Length..norm:Prop.overlap.focal.large.norm"] * Length..norm) / 
#       sd(focal_data$Length..cm.)
#   )
# 
# f_l_grad_plt <- ggplot(focal_large_gradients, aes(x = Length..cm., y = Large.gradient)) +
#   geom_line() +
#   theme_jun1()
# 
# f_l_grad_plt
# 
# ggplot(focal_data, aes(x = Length..cm., y = Prop.overlap.focal.large)) +
#   geom_point() +
#   geom_line(data = tibble(
#     Length..cm. = seq(min(focal_data$Length..cm.), max(focal_data$Length..cm.), length.out = 50)
#   ) %>% mutate(Prop.overlap.focal.large = pi * Length..cm.^2))
# 
# 
# ## ====== Plot heatmap of effect from focal neighbours =====
# 
# focal_plot <- focal_data$Plot[1]
# focal_year <- focal_data$Year[1]
# 
# plt_design <- "
#   12
#   13
# "
# 
# 
# 
# # Scale for heatmap colour
# color_lims <- if(FOCAL_VR == "Growth.future") {
#   c(
#     -max(abs(quant_vr[2]), abs(quant_vr[6])),
#     max(abs(quant_vr[2]), abs(quant_vr[6]))
#   )
# } else if(FOCAL_VR == "X..Capitulescences") {
#   c(0, quant_vr[7])
# } else if(FOCAL_VR == "Survival.future") {
#   c(0, 1)
# }
# color_lims
# 
# 
# 
# # Expand grid
# focal_data_grad <- expand.grid(Prop.overlap.focal.small = os_seq,
#                                Prop.overlap.focal.large = ol_seq) %>%
#   mutate( # Add arbitrary plot & year and median for non-focal overlap
#     Plot = focal_plot, Year = focal_year,
#     Length..cm. = quant_length[5], Prop.overlap.other = quant_oo[4]
#   )
# # Generate model predictions
# focal_data_grad[[FOCAL_VR]] <- predict(stat_model, focal_data_grad)
# 
# focal_data_grad
# 
# # Plot as heatmap
# focal_heatmap <- ggplot(focal_data_grad, aes(x = Prop.overlap.focal.small,
#                                              y = Prop.overlap.focal.large,
#                                              fill = .data[[FOCAL_VR]])) +
#   geom_tile() +
#   scale_fill_distiller(palette = "RdBu", limits = color_lims) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme_jun1()
# 
# focal_heatmap
# 
# # Plot line graphs to illustrate how vital rate responds to crowding
# focal_data_line <- tibble(
#   Plot = focal_plot,
#   Year = focal_year,
#   Length..cm. = quant_length[3],
#   Prop.overlap = o_seq
# )
# focal_data_line[[overlap_type]] <- o_seq
# focal_data_line[[other_type]] <- quant_other[4]
# focal_data_line[[FOCAL_VR]] <- predict(stat_model, focal_data_line)
# 
# focal_line_small <- ggplot(focal_data_line, aes(x = .data[[overlap_type]],
#                                                 y = .data[[FOCAL_VR]])) +
#   geom_line() +
#   coord_cartesian(ylim = c(quant_vr[2], quant_vr[6])) +
#   theme_jun1() +
#   labs(
#     caption = paste0(
#       "Length of plant affected (cm) = ", quant_length[3]
#     )
#   )
# 
# focal_data_line <- focal_data_line %>%
#   mutate(
#     Length..cm. = quant_length[5]
#   )
# focal_data_line[[FOCAL_VR]] <- predict(stat_model, focal_data_line)
# 
# focal_line_large <- ggplot(focal_data_line, aes(x = .data[[overlap_type]],
#                                                 y = .data[[FOCAL_VR]])) +
#   geom_line() +
#   coord_cartesian(ylim = c(quant_vr[2], quant_vr[6])) +
#   theme_jun1() +
#   labs(
#     caption = paste0(
#       "Length of plant affected (cm) = ", quant_length[5]
#     )
#   )
# 
# # Merge plots together
# focal_summ <- focal_heatmap + focal_line_small + focal_line_large +
#   plot_layout(design = plt_design) +
#   plot_annotation(
#     title = paste0(
#       "Effect of", if(overlap_type == "Prop.overlap.focal") { "" } else { " NOT" },
#       " ",
#       if(FOCAL_STAGE == "s") { "small" } else { "large" },
#       " individuals of ", SP_NEI,
#       " on vital rate of ", SP_FOCAL
#     ),
#     caption = paste0(
#       "Plot = ", focal_plot,
#       ", Year = ", focal_year,
#       ", Overlap from ",
#       if(overlap_type == "Prop.overlap.focal") { "other" } else { "focal" },
#       " individuals = ", quant_other[3]
#     )
#   )
# 
# return(focal_summ)
# 
# # Function for generating a summary plot
# gen_summ_plt <- function(overlap_type) {
#   other_type <- if(overlap_type == "Prop.overlap.focal") { "Prop.overlap.other" }
#     else { "Prop.overlap.focal" }
# 
#   # Get quantiles of variables
#   quant_length <- quantile(focal_data$Length..cm., c(0.025, 0.1, 0.5, 0.9, 0.95))
#   quant_o <- quantile(focal_data[[overlap_type]], c(0.05, 0.1, 0.5, 0.9, 0.95))
#   quant_other <- quantile(focal_data[[other_type]], c(0.05, 0.1, 0.5, 0.9, 0.95))
#   quant_vr <- quantile(focal_data[[FOCAL_VR]], c(0.05, 0.1, 0.5, 0.9, 0.95))
#   # Gradient of model output to generate heatmap
#   nstep <- 50
#   nstep_length <- nstep
#   nstep_o <- nstep
#   # Generate sequences
#   length_seq <- seq(quant_length[1], quant_length[5], length.out = nstep_length)
#   o_seq <- seq(quant_o[1], quant_o[5], length.out = nstep_o)
#   # Expand grid
#   focal_data_grad <- expand.grid(Length..cm. = length_seq,
#                                 Prop.overlap = o_seq) %>%
#     mutate( # Add arbitrary plot & year and median for non-focal overlap
#       Plot = focal_plot, Year = focal_year
#     )
#   focal_data_grad[[other_type]] <- quant_other[3]
#   focal_data_grad[[overlap_type]] <- focal_data_grad$Prop.overlap
#   # Generate model predictions
#   focal_data_grad[[FOCAL_VR]] <- predict(stat_model, focal_data_grad)
# 
#   # Plot as heatmap
#   focal_heatmap <- ggplot(focal_data_grad, aes(x = .data[[overlap_type]],
#                                                y = Length..cm.,
#                                                fill = .data[[FOCAL_VR]])) +
#     geom_tile() +
#     scale_fill_distiller(palette = "RdBu", limits = c(
#       -max(abs(quant_vr[1]), abs(quant_vr[5])),
#       max(abs(quant_vr[1]), abs(quant_vr[5]))
#     )) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme_jun1()
# 
#   # Plot line graphs to illustrate how vital rate responds to crowding
#   focal_data_line <- tibble(
#     Plot = focal_plot,
#     Year = focal_year,
#     Length..cm. = quant_length[2],
#     Prop.overlap = o_seq
#   )
#   focal_data_line[[overlap_type]] <- o_seq
#   focal_data_line[[other_type]] <- quant_other[3]
#   focal_data_line[[FOCAL_VR]] <- predict(stat_model, focal_data_line)
# 
#   focal_line_small <- ggplot(focal_data_line, aes(x = .data[[overlap_type]],
#                                                   y = .data[[FOCAL_VR]])) +
#     geom_line() +
#     coord_cartesian(ylim = c(quant_vr[1], quant_vr[5])) +
#     theme_jun1() +
#     labs(
#       caption = paste0(
#         "Length of plant affected (cm) = ", quant_length[2]
#       )
#     )
# 
#   focal_data_line <- focal_data_line %>%
#     mutate(
#       Length..cm. = quant_length[4]
#     )
#   focal_data_line[[FOCAL_VR]] <- predict(stat_model, focal_data_line)
# 
#   focal_line_large <- ggplot(focal_data_line, aes(x = .data[[overlap_type]],
#                                                   y = .data[[FOCAL_VR]])) +
#     geom_line() +
#     coord_cartesian(ylim = c(quant_vr[1], quant_vr[5])) +
#     theme_jun1() +
#     labs(
#       caption = paste0(
#         "Length of plant affected (cm) = ", quant_length[4]
#       )
#     )
# 
#   # Merge plots together
#   focal_summ <- focal_heatmap + focal_line_small + focal_line_large +
#     plot_layout(design = plt_design) +
#     plot_annotation(
#       title = paste0(
#         "Effect of", if(overlap_type == "Prop.overlap.focal") { "" } else { " NOT" },
#         " ",
#         if(FOCAL_STAGE == "s") { "small" } else { "large" },
#         " individuals of ", SP_NEI,
#         " on vital rate of ", SP_FOCAL
#       ),
#       caption = paste0(
#         "Plot = ", focal_plot,
#         ", Year = ", focal_year,
#         ", Overlap from ",
#         if(overlap_type == "Prop.overlap.focal") { "other" } else { "focal" },
#         " individuals = ", quant_other[3]
#       )
#     )
# 
#   return(focal_summ)
# }
# 
# gen_summ_plt("Prop.overlap.focal")
# 
# ggsave(here("StageDependentComp", "figures", "census_analysis", "overlap",
#             paste0(FOCAL_VR_ABB, "_", SP_FOCAL_ABB, "_by_",
#                    SP_NEI_ABB, "_", FOCAL_STAGE, ".png")),
#        width = 8, height = 6)
# 
# gen_summ_plt("Prop.overlap.other")
# 
# ggsave(here("StageDependentComp", "figures", "census_analysis", "overlap",
#             paste0(FOCAL_VR_ABB, "_", SP_FOCAL_ABB, "_by_",
#                    SP_NEI_ABB, "_", FOCAL_STAGE, "_other", ".png")),
#        width = 8, height = 6)

