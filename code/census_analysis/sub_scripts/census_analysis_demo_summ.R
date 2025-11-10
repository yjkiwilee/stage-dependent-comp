# Demographic summary of species in the system

census_data <- read_csv(here("StageDependentComp",
                             "result_data",
                             "data_cleanup",
                             "demography_2014-2022_processed.csv"))

census_spp_vr_summary <- census_data %>%
  group_by(Taxon) %>%
  filter(Died.this.census.final == 0) %>%
  group_modify(
    function(df, key) {
      df <- df %>%
        filter(
          !is.na(Survival.future) & !is.na(Stage.future)
        )
      
      return(tibble(
        Juvenile.survival = sum(df[df$Stage == "S",]$Survival.future) /
          nrow(df[df$Stage == "S",]),
        Adult.survival = sum(df[df$Stage == "L",]$Survival.future) /
          nrow(df[df$Stage == "L",]),
        Progression = nrow(df[df$Stage == "S" & df$Stage.future == "L" &
                                   df$Survival.future == 1,]) / nrow(df),
        Retrogression = nrow(df[df$Stage == "L" & df$Stage.future == "S" &
                                     df$Survival.future == 1,]) / nrow(df)
      ))
    }
  )

census_spp_fert_summary <- census_data %>%
  group_by(Taxon, Year) %>%
  group_modify(
    function(df, key) {
      return(tibble(
        N.Seedlings = nrow(df[df$Is.seedling == 1,]),
        N.Flowering = nrow(df[df$Is.flowering == 1,])
      ))
    }
  )

census_spp_fert_summary_future <- census_spp_fert_summary %>%
  mutate(Year = Year - 1)

census_spp_fert_summary <- census_spp_fert_summary %>%
  left_join(census_spp_fert_summary_future, by = join_by(Taxon, Year),
            suffix = c("", ".Future")) %>%
  mutate(
    Fertility = N.Seedlings.Future / N.Flowering
  ) %>%
  group_by(Taxon) %>%
  summarise(
    Fertility = mean(Fertility, na.rm = TRUE)
  )

census_spp_vr_summary <- census_spp_vr_summary %>%
  left_join(census_spp_fert_summary, by = join_by(Taxon))

census_spp_vr_summary <- census_spp_vr_summary %>%
  rename(
    sj = Juvenile.survival,
    sa = Adult.survival,
    g = Progression,
    r = Retrogression,
    f = Fertility
  )

census_spp_vr_summary_new <- NULL

for(row_i in 1:nrow(census_spp_vr_summary)) {
  spp_vrs_row <- census_spp_vr_summary[row_i,]
  
  vrs <- list(
    s = c(spp_vrs_row$sj[1], spp_vrs_row$sa[1]),
    g = c(spp_vrs_row$g[1]),
    r = c(spp_vrs_row$r[1]),
    f = c(spp_vrs_row$f[1])
  )
  
  sp_mpm <- vrs_to_mpm(2, vrs)
  
  if(!(NaN %in% sp_mpm$matA) && !(Inf %in% sp_mpm$matA)) {
    spp_vrs_row_new <- spp_vrs_row %>%
      mutate(
        m_1_1 = sp_mpm$matA[1,1],
        m_1_2 = sp_mpm$matA[1,2],
        m_2_1 = sp_mpm$matA[2,1],
        m_2_2 = sp_mpm$matA[2,2],
        gen_time = generation.time(sp_mpm$matA),
        net_rep_rate = net.reproductive.rate(sp_mpm$matA),
        life_expect = life_expect_mean(sp_mpm$matU),
        rep_val_j = reproductive.value(sp_mpm$matA)[1],
        rep_val_a = reproductive.value(sp_mpm$matA)[2],
        sens_dd_vr = d_lambda(dd_vr, vrs),
        juv_sens_dd_vr = d_prop_juv(sj, sa, g, r, f, dd_vr),
        max_real_lambda = max(
          f + sa + (sj - sa - f) * 1,
          f + sa + (sj - sa - f) * 0
        ),
        min_real_lambda = min(
          f + sa + (sj - sa - f) * 1,
          f + sa + (sj - sa - f) * 0
        ),
        diff_real_lambda = abs(
          sj - sa - f
        )
      ) %>%
      select(
        !(starts_with("sens_") & (!sens_dd_vr))
      )
  }
  
  if(is.null(census_spp_vr_summary_new)) {
    census_spp_vr_summary_new <- spp_vrs_row_new
  } else {
    census_spp_vr_summary_new <- bind_rows(census_spp_vr_summary_new, spp_vrs_row_new)
  }
}

census_spp_vr_summary <- census_spp_vr_summary_new

write_csv(census_spp_vr_summary,
          here("StageDependentComp",
               "result_data",
               "wp2",
               "spp_demo_summary.csv"))


