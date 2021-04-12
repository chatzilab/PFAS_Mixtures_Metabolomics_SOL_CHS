#read in temporary dataset
solar_c18_mwas_results_1 = read_rds(file = here::here('Temporary results',
                                                   '1.1.0_SOLAR_c18_mwas.rds')) 
solar_hilic_mwas_results_1 = read_rds(file =  here::here('Temporary results',
                                                      '1.1.0_SOLAR_hilic_mwas.rds')) 

solar_c18_n_metab <- table(solar_c18_mwas_results_1$exposure)[1] 
solar_hilic_n_metab <- table(solar_hilic_mwas_results_1$exposure)[1] 

#BH adjusted p-value - Check with Jesse
solar_c18_mwas_results_2 <-  solar_c18_mwas_results_1 %>% 
  dplyr::group_by(exposure) %>%
  dplyr::arrange(exposure) %>%
  dplyr::mutate(p.value = as.numeric(p.value)) %>%
  dplyr::mutate(t.score = as.numeric(t.score)) %>%
  dplyr::mutate(p.value_BH = p.adjust(p.value, 
                                      method = 'BH', 
                                      n = solar_c18_n_metab)) 
 
#some summary
solar_c18_mwas_results <- solar_c18_mwas_results_2 %>% 
  dplyr::group_by(exposure) %>%
  dplyr::summarise(n = n(),
                   pvalue_mean = mean(p.value),
                  pvalue_BH_mean = mean(p.value_BH),
                  perc_pval_0.05 = sum(p.value < 0.05)/length(p.value),
                  perc_pval_0.05_BH = sum(p.value_BH <0.05) / length(p.value),
                  n_pval_0.05_BH = sum(p.value_BH <0.05))
#For HILIC
solar_hilic_mwas_results_2 <- solar_hilic_mwas_results_1 %>% 
  dplyr::group_by(exposure) %>%
  dplyr::arrange(exposure) %>%
  dplyr::mutate(p.value = as.numeric(p.value)) %>%
  dplyr::mutate(t.score = as.numeric(t.score)) %>%
  dplyr::mutate(p.value_BH = p.adjust(p.value, 
                                      n = solar_hilic_n_metab, 
                                      method = 'BH')) # n = ### to adjust p value, mutate function is not correct


solar_hilic_mwas_results <- solar_hilic_mwas_results_2 %>% 
  dplyr::group_by(exposure) %>%
  dplyr::summarise(n = n(),
                   pvalue_mean = mean(p.value),
                   pvalue_BH_mean = mean(p.value_BH),
                   pec_pval_0.05 = sum(p.value < 0.05)/length(p.value),
                   pec_pval_0.05_BH = sum(p.value_BH <0.05) / length(p.value),
                   n_pval_0.05_BH = sum(p.value_BH <0.05))

rm(solar_c18_n_metab, solar_hilic_n_metab, 
   solar_c18_mwas_results_1, solar_hilic_mwas_results_1, 
   solar_c18_mwas_results_2, solar_hilic_mwas_results_2 )
   #chs_c18_n_metab, chs_c18_mwas_results_1, chs_hilic_mwas_results_2,)

write_csv(solar_c18_mwas_results, here::here("Temporary results", 
                                             "solar_c18_mwas_results.csv"))

write_csv(solar_hilic_mwas_results, here::here("Temporary results", 
                                             "solar_hilic_mwas_results.csv"))
