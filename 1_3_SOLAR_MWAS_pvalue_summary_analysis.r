#read in temporary dataset
solar_c18_mwas_results_1 = read_rds(file = here::here('Temporary results', 
                                                      exposure_type,
                                                      '1_1_SOLAR_c18_mwas.rds')) 
solar_hilic_mwas_results_1 = read_rds(file =  here::here('Temporary results', 
                                                         exposure_type,
                                                         '1_1_SOLAR_hilic_mwas.rds')) 
# Get number of metabolites for HELIC and C18
solar_c18_n_metab <- table(solar_c18_mwas_results_1$exposure)[1] 
solar_hilic_n_metab <- table(solar_hilic_mwas_results_1$exposure)[1] 

# C18 ---------------------------------

## Calculate BH adjusted p-value, calculate summary statistics
solar_c18_mwas_results <-  solar_c18_mwas_results_1 %>% 
  group_by(exposure) %>%
  mutate(p.value = as.numeric(p.value),
         t.score = as.numeric(t.score),
         p.value_BH = p.adjust(p.value, 
                               method = 'BH', 
                               n = solar_c18_n_metab)) %>% 
  summarise(n = n(),
            pvalue_mean = mean(p.value),
            pvalue_BH_mean = mean(p.value_BH),
            perc_pval_0.05 = sum(p.value < 0.05)/length(p.value),
            perc_pval_0.05_BH = sum(p.value_BH <0.05) / length(p.value),
            n_pval_0.05_BH = sum(p.value_BH <0.05))

# HILIC -----------------------------------------
## Calculate BH adjusted p-value, calculate summary statistics
solar_hilic_mwas_results <- solar_hilic_mwas_results_1 %>% 
  group_by(exposure) %>%
  mutate(p.value = as.numeric(p.value),
         t.score = as.numeric(t.score),
         p.value_BH = p.adjust(p.value, 
                               n = solar_hilic_n_metab, 
                               method = 'BH')) %>% 
  summarise(n = n(),
            pvalue_mean = mean(p.value),
            pvalue_BH_mean = mean(p.value_BH),
            pec_pval_0.05 = sum(p.value < 0.05)/length(p.value),
            pec_pval_0.05_BH = sum(p.value_BH <0.05) / length(p.value),
            n_pval_0.05_BH = sum(p.value_BH <0.05))


rm(solar_c18_n_metab, solar_hilic_n_metab, 
   solar_c18_mwas_results_1, solar_hilic_mwas_results_1)
   # chs_c18_n_metab, chs_c18_mwas_results_1, chs_hilic_mwas_results_2,)

write_csv(solar_c18_mwas_results, here::here("Temporary results", 
                                             exposure_type,
                                             "solar_c18_mwas_results.csv"))

write_csv(solar_hilic_mwas_results, here::here("Temporary results", 
                                               exposure_type,
                                               "solar_hilic_mwas_results.csv"))