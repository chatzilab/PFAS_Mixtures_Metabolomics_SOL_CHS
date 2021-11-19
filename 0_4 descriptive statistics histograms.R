# 0_4 descriptive statistics and histograms of outcomes 


exposure_outcome_long <- exposure_outcome %>% 
  modify(~dplyr::select(., id, all_of(exposures_continuous)) %>% 
           pivot_longer(., cols = all_of(exposures_continuous), 
                        names_to = "exposure",
                        values_to = "concentration") %>% 
           mutate(id = as.character(id))) %>%
  bind_rows() %>% 
  mutate(cohort = if_else(str_detect(id, "sol"), "solar", "chs")) %>% 
  group_by(cohort, exposure) %>% 
  mutate(lg2_concentration = log2(concentration)) %>% 
  ungroup()




# write.csv(xxx, "temp.csv")


xxx <- exposure_outcome_long %>% 

ggplot(exposure_outcome_long %>% filter(cohort == "solar"), 
       aes(x = scale(log2(concentration), center = TRUE))) + 
  geom_histogram(color= "grey20", fill = "grey90") + 
  facet_wrap(~exposure, scales = "free")
