#read in temporary dataset
chs_c18_mwas_results_1 = read_rds(file = here::here('Temporary results',
                                                    exposure_type,
                                                    '1_1_chs_c18_mwas.rds')) 
chs_hilic_mwas_results_1 = read_rds(file =  here::here('Temporary results',
                                                       exposure_type,
                                                       '1_1_chs_hilic_mwas.rds')) 

chs_c18_n_metab <- table(chs_c18_mwas_results_1$exposure)[1] 
chs_hilic_n_metab <- table(chs_hilic_mwas_results_1$exposure)[1] 

#BH adjusted p-value - Check with Jesse
chs_c18_mwas_results <-  chs_c18_mwas_results_1 %>% 
  group_by(exposure) %>%
  arrange(exposure) %>%
  mutate(p.value = as.numeric(p.value)) %>%
  mutate(t.score = as.numeric(t.score)) %>%
  mutate(p.value_BH = p.adjust(p.value, 
                               method = 'BH', 
                               n = chs_c18_n_metab)) %>% 
  summarise(n = n(),
            pvalue_mean = mean(p.value),
            pvalue_BH_mean = mean(p.value_BH),
            perc_pval_0.05 = sum(p.value < 0.05)/length(p.value),
            perc_pval_0.05_BH = sum(p.value_BH <0.05) / length(p.value),
            n_pval_0.05_BH = sum(p.value_BH <0.05))

#For HILIC
chs_hilic_mwas_results <- chs_hilic_mwas_results_1 %>% 
  group_by(exposure) %>%
  mutate(p.value = as.numeric(p.value)) %>%
  mutate(t.score = as.numeric(t.score)) %>%
  mutate(p.value_BH = p.adjust(p.value, 
                               n = chs_hilic_n_metab, 
                               method = 'BH')) %>%
  summarise(n = n(),
            pvalue_mean = mean(p.value),
            pvalue_BH_mean = mean(p.value_BH),
            pec_pval_0.05 = sum(p.value < 0.05)/length(p.value),
            pec_pval_0.05_BH = sum(p.value_BH <0.05) / length(p.value),
            n_pval_0.05_BH = sum(p.value_BH <0.05))

# Clean workspace
rm(chs_c18_n_metab, chs_hilic_n_metab, 
   chs_c18_mwas_results_1, chs_hilic_mwas_results_1)


# Save data
write_csv(chs_c18_mwas_results, here::here("Temporary results", exposure_type,
                                           "chs_c18_mwas_results.csv"))

write_csv(chs_hilic_mwas_results, here::here("Temporary results", exposure_type,
                                             "chs_hilic_mwas_results.csv"))
