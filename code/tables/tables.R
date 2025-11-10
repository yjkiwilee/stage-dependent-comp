# Script for generating tables in the text

# Packages
pacman::p_load("tidyverse", "here", "ggpubr", "patchwork", "cowplot", "flextable",
               "officer")

# Load ggplot themes
source(here("StageDependentComp", "code", "themes", "custom_themes.R"))

init_flextable_defaults()

set_flextable_defaults(
  font.family = "Arial",
  digits = 3,
  big.mark = ",",
  na_str = "<na>",
  border.color = "#222"
)

# ====== Table 1 ======

sp_names <- c("F2", "F1", "M", "S1", "S2")

realised_vr_df <- read_csv(
  here("StageDependentComp","result_data","wp3","realised_spp_vrs.csv")
)

spp_vrs_df <- read_csv(
  here("StageDependentComp","result_data","wp3","spp_vrs.csv")
)


spp_vrs_df_summ <- spp_vrs_df %>%
  mutate(
    sp_name = factor(sp_name, levels = sp_names)
  ) %>%
  group_by(sp_name) %>%
  summarise(
    sj = mean(sj),
    sa = mean(sa),
    g = mean(g),
    r = mean(r),
    min_f = min(f),
    max_f = max(f)
  ) %>%
  select(!sp_name)

spp_vrs_tab <- t(spp_vrs_df_summ)
colnames(spp_vrs_tab) <- sp_names
rownames(spp_vrs_tab) <- NULL
spp_vrs_tab <- as.data.frame(spp_vrs_tab)

spp_vrs_tab <- tibble(
    `Baseline vital rate` = c(
      "Juvenile survival",
      "Adult survival",
      "Progression",
      "Retrogression",
      "Minimum fertility",
      "Maximum fertility"
    )
  ) %>%
  bind_cols(spp_vrs_tab) %>%
  rename_with(
    function(vn) { paste0("Life history strategy_", vn) },
    !all_of(c("Baseline vital rate"))
  )

head_b <- fp_border(color="#222", width = 2)

spp_vrs_tab_rend <- flextable(spp_vrs_tab) %>%
  separate_header() %>%
  hline(i = 1, j = 2, border = head_b, part = "header") %>%
  align(i = 1, j = 1, align = "center", part = "header") %>%
  valign(i = 1, j = 1, valign = "center", part = "header") %>%
  bold(part = "header") %>%
  colformat_double() %>%
  autofit()

spp_vrs_tab_rend

save_as_image(
  spp_vrs_tab_rend,
  here(
    "StageDependentComp",
    "figures",
    "comm_sims",
    "spp_vrs_table.png"
  )
)







