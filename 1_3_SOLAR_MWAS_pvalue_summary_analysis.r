#read in temporary dataset

# SOLAR
solar_mwas_results <- read_rds(
  file = here::here('Temporary results',
                    exposure_type,
                    '1_1_SOLAR_mwas_results_two_hour_met.rds')) %>% 
  modify(~mutate(., 
                 mz = str_split_fixed(name, "_", n = 2)[,1] %>% as.numeric(.),
                 rt = str_split_fixed(name, "_", n = 2)[,2] %>% as.numeric(.)))


solar_mwas_results_df <- solar_mwas_results %>%
  bind_rows(.id = "mode")

# CHS
chs_mwas_results <- read_rds(
  file = here::here('Temporary results',
                    exposure_type,
                    '1_1_CHS_mwas_results_two_hour_met.rds')) %>% 
  modify(~mutate(., 
                 mz = str_split_fixed(name, "_", n = 2)[,1] %>% as.numeric(.),
                 rt = str_split_fixed(name, "_", n = 2)[,2] %>% as.numeric(.)))

chs_mwas_results_df <- chs_mwas_results %>%
  bind_rows(.id = "mode")


# Volcano plots -----------------------
# # SOLAR 
# for(i in names(solar_mwas_results)){
#   figout <- ggplot(solar_mwas_results[[i]], 
#                    aes(x = estimate, 
#                        y = -1*log(p_value), 
#                        color = significancefdr)) +
#     geom_jitter(alpha = .4, size = .25) + 
#     facet_wrap(~exposure)
#   
#   ggsave(figout, filename = fs::path(here::here("figures", 
#                                                 "volcano plots",
#                                                 paste("solar", i, "volcano.jpeg", sep = "_"))))
# }
# 
# 
# # CHS 
# for(i in names(chs_mwas_results)){
#   figout <- ggplot(chs_mwas_results[[i]], 
#                    aes(x = estimate, 
#                        y = -1*log(p_value), 
#                        color = significancefdr)) +
#     geom_jitter(alpha = .4, size = .25) + 
#     facet_wrap(~exposure)
#   
#   ggsave(figout, filename = fs::path(here::here("figures", 
#                                                 "volcano plots",
#                                                 paste("chs", i, "volcano.jpeg", sep = "_"))))
# }




# Manhattan plots -----------------------
# SOLAR 
# # for(i in names(solar_mwas_results)){
#   ggplot(solar_mwas_results_df %>% filter(str_detect(exposure, "pfhxs")), 
#                    aes(x = mz, 
#                        y = -1*log(p_value), 
#                        color = significancefdr)) +
#     geom_point(alpha = .9, size = .75) +
#     facet_wrap(~mode, scales = "free_x") +
#   ylim(c(0,30))
#   
#   
#   # ggsave(figout, filename = fs::path(here::here("figures", 
#   #                                               "volcano plots",
#   #                                               paste("solar", i, "volcano.jpeg", sep = "_"))))
# # }
# 
# # CHS 
# for(i in names(chs_mwas_results)){
#   figout <- ggplot(chs_mwas_results[[i]], 
#                    aes(x = estimate, 
#                        y = -1*log(p_value), 
#                        color = significancefdr)) +
#     geom_jitter(alpha = .4, size = .25) + 
#     facet_wrap(~exposure)
#   
#   ggsave(figout, filename = fs::path(here::here("figures", 
#                                                 "volcano plots",
#                                                 paste("chs", i, "volcano.jpeg", sep = "_"))))
# }







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