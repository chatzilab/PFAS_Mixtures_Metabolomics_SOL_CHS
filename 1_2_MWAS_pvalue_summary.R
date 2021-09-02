#read in temporary dataset

# SOLAR
solar_mwas_results <- read_rds(
  file = here::here('Temporary results',
                    exposure_type,
                    '1_1_SOLAR_mwas_results_two_hour_met.rds')) %>% 
  modify(~mutate(., 
                 mz = str_split_fixed(name, "_", n = 2)[,1] %>% 
                   as.numeric(.),
                 rt = str_split_fixed(name, "_", n = 2)[,2] %>% 
                   as.numeric(.)))


solar_mwas_results_df <- solar_mwas_results %>%
  bind_rows(.id = "mode")

# CHS
chs_mwas_results <- read_rds(
  file = here::here('Temporary results',
                    exposure_type,
                    '1_1_CHS_mwas_results_two_hour_met.rds')) %>% 
  modify(~mutate(., 
                 mz = str_split_fixed(name, "_", n = 2)[,1] %>% 
                   as.numeric(.),
                 rt = str_split_fixed(name, "_", n = 2)[,2] %>% 
                   as.numeric(.)))

chs_mwas_results_df <- chs_mwas_results %>%
  bind_rows(.id = "mode")


# tables -----------------------

# SOLAR 
(solar_summary_mwas <- solar_mwas_results_df %>% 
  group_by(mode, exposure) %>% 
  summarise(n_features = length(name), 
            percent_significant_p05 = npct(significance, 
                                           "P<0.05", 
                                           n.digits = 2), 
            percent_significant_q05 = npct(significancefdr, 
                                           "P<0.05", 
                                           n.digits = 2)))

# chs 
(chs_summary_mwas <- chs_mwas_results_df %>% 
    group_by(mode, exposure) %>% 
    summarise(n_features = length(name), 
              percent_significant_p05 = npct(significance, 
                                             "P<0.05", 
                                             n.digits = 2), 
              percent_significant_q05 = npct(significancefdr, 
                                             "P<0.05", 
                                             n.digits = 2)))


# combine solar chs
mwas_results_summary <- full_join(solar_summary_mwas, chs_summary_mwas, 
                                  by = c("mode", "exposure", "n_features"), 
                                  suffix = c("_solar", "_chs"))

write_excel_csv(mwas_results_summary, 
          file = fs::path(dir_reports, "Summary of MWAS Results.csv"))

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