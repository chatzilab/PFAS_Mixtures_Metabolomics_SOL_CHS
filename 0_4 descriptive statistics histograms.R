library(ggridges)

# 0_4 descriptive statistics and histograms of outcomes 

exposure_outcome$chs$wave <- "MetaAir"

exposure_outcome_long <- exposure_outcome %>% 
  modify(~dplyr::select(., id, wave, all_of(exposures_continuous)) %>% 
           pivot_longer(., cols = all_of(exposures_continuous), 
                        names_to = "exposure",
                        values_to = "concentration") %>% 
           mutate(id = as.character(id))) %>%
  bind_rows() %>% 
  mutate(cohort = if_else(str_detect(id, "sol"), "SOLAR", "CHS"), 
         cohort_wave = str_c(cohort, wave, sep = ", ") %>%
           str_remove(", MetaAir") %>% 
           fct_relevel("SOLAR, first wave", "SOLAR, second wave"),
         cohort_wave_pfas = str_c(cohort, wave, exposure,sep = ", ")) %>% 
  group_by(cohort, exposure) %>% 
  mutate(lg2_concentration = log2(concentration)) %>% 
  ungroup()


# Rename PFAS
exposure_outcome_long <- exposure_outcome_long %>% 
  mutate(pfas = rename_pfas(exposure) %>% 
           fct_reorder(concentration, .desc = TRUE))


# Log transformed density plots
(sol_chs_densityplot_lg_var <- ggplot(exposure_outcome_long,
                                      aes(x = concentration, 
                                          y = fct_rev(cohort_wave))) + 
    # geom_density(alpha = .5, color = "black", ) + 
    geom_density_ridges(scale = 1.4) +
    theme(legend.position = "bottom", 
          legend.justification = "center") +
    facet_wrap(~pfas, scales = "free_x") +
    labs(fill='Study Wave') +
    scale_x_log10(expand=expansion(mult = c(0,.2))) +
    xlab("Concentration (µg/L)") +
    ylab(NULL))

# Save
ggsave(sol_chs_densityplot_lg_var, 
       filename = fs::path(dir_reports, "PFAS exposure density plot lg trns.jpg"),
       width = 9, height = 6, units = "in")

# Density plots (original scale)
(sol_chs_densityplot <- ggplot(exposure_outcome_long,
                               aes(x = concentration,
                                   y = fct_rev(cohort_wave))) +
    geom_density_ridges(scale = 1.4) +
    theme(legend.position = "bottom", 
          legend.justification = "center") +
    facet_wrap(~pfas, scales = "free_x") +
    labs(fill='Study Wave') +
    xlab("Concentration (µg/L)") +
    ylab(NULL))



ggsave(sol_chs_densityplot, 
       filename = fs::path(dir_reports, "PFAS exposure density plot raw scale.jpg"),
       width = 9, height = 6, units = "in")


# Test normality
normality_test_results_lg2 <- tapply(exposure_outcome_long$lg2_concentration, 
                                     exposure_outcome_long$cohort_wave_pfas, 
                                     shapiro.test)

normality_test_results_raw <- tapply(exposure_outcome_long$concentration, 
                                     exposure_outcome_long$cohort_wave_pfas, 
                                     shapiro.test)

# create df
normality_test_df <- modify2(
  normality_test_results_raw,
  normality_test_results_lg2, 
  ~data.frame(#statistic_raw = .x$statistic, 
    p_raw = .x$p, 
    # statistic_lg = .y$statistic, 
    p_lg = .y$p)) %>% 
  bind_rows(.id = "name")


normality_test_df$name <- names(normality_test_results_raw)

normality_test_df <- normality_test_df %>% 
  mutate(sig_raw = if_else(p_raw < 0.01, "sig", "not sig."), 
         sig_lg  = if_else(p_lg  < 0.01, "sig", "not sig."))




# Standard deviations
pfas_sd <- data.frame(sd=tapply(exposure_outcome_long$concentration, 
                                exposure_outcome_long$cohort_wave_pfas, 
                                sd))

