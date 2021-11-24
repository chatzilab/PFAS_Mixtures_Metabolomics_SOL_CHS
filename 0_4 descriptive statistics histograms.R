# 0_4 descriptive statistics and histograms of outcomes 
exposure_outcome$chs$wave <- "MetaAir"

exposure_outcome_long <- exposure_outcome %>% 
  modify(~dplyr::select(., id, wave, all_of(exposures_continuous)) %>% 
           pivot_longer(., cols = all_of(exposures_continuous), 
                        names_to = "exposure",
                        values_to = "concentration") %>% 
           mutate(id = as.character(id))) %>%
  bind_rows() %>% 
  mutate(cohort = if_else(str_detect(id, "sol"), "solar", "chs")) %>% 
  group_by(cohort, exposure) %>% 
  mutate(lg2_concentration = log2(concentration)) %>% 
  ungroup()

# Rename PFAS
exposure_outcome_long <- exposure_outcome_long %>% 
  mutate(pfas = rename_pfas(exposure) %>% 
           fct_reorder(concentration, .desc = TRUE))


# SOLAR 
(sol_hist <- ggplot(exposure_outcome_long %>% filter(cohort == "solar"), 
                    aes(x = concentration, center = TRUE, fill = wave)) + 
    geom_histogram(color= "grey20", bins = 30) + 
    theme(legend.position = "bottom", 
          legend.justification = "center") +
    facet_wrap(~pfas, scales = "free") +
    labs(fill='Study Wave') +
    xlab("Concentration (µg/L)") + 
    ylab("Count"))



# CHS 
(chs_hist <- ggplot(exposure_outcome_long %>% filter(cohort == "chs"), 
                    aes(x = concentration, center = TRUE)) + 
    geom_histogram(color= "grey20", fill = "darkgreen", bins = 30) + 
    facet_wrap(~pfas, scales = "free") +
    xlab("Concentration (µg/L)") + 
    ylab("Count") )


# Arrange plots 



sol_chs_hist <- cowplot::plot_grid(NULL, 
                   sol_hist, 
                   NULL,
                   chs_hist, 
                   ncol = 1, 
                   rel_heights = c(.05, 1, .05, .9),
                   labels = c("A) SOLAR", "",
                              "B) CHS", ""))

ggsave(sol_chs_hist, 
       filename = fs::path(dir_reports, "PFAS exposure histogram.jpg"),
       width = 7, height = 10, units = "in")
