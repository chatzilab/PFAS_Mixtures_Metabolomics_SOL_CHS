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

write_excel_csv(mwas_results_summarysults_summary, 
          file = fs::path(dir_reports, "Summary of MWAS Results.csv"))






# Merge to get single data frame of all beta coefficients 
sol_mwas_betas <- solar_mwas_results %>% 
  modify(~select(.x, exposure, name, estimate) %>% 
           pivot_wider(id_cols = name, 
                       names_from = exposure, 
                       values_from = estimate)) %>% 
  bind_rows(.id = "mode")


chs_mwas_betas <- chs_mwas_results %>% 
  modify(~select(.x, exposure, name, estimate) %>% 
           pivot_wider(id_cols = name, 
                       names_from = exposure, 
                       values_from = estimate)) %>% 
  bind_rows(.id = "mode")

# Make list of MWAS results
mwas_results_long <- list(solar = solar_mwas_results %>% 
                            bind_rows(.id = "mode"), 
                          chs   = chs_mwas_results %>% 
                            bind_rows(.id = "mode"))

# 
write_rds(mwas_results_long, 
          here::here("Temporary results", 
                     "PFAS", 
                     "SOL CHS all MWAS results long.rds"))



# PCA analysis on the effect estimates
library(ggfortify)
temp <- chs_mwas_betas %>% select(-mode, -name)
pca_res <- prcomp(temp, scale. = TRUE)

autoplot(pca_res, data = chs_mwas_betas, colour = 'mode', alpha = .3)



# Save MWAS results 
mwas_beta_coefs <- list(solar = sol_mwas_betas, 
                        chs   = chs_mwas_betas)

write_rds(mwas_beta_coefs, 
          here::here("Temporary results", 
                     "PFAS", 
                     "SOL CHS all MWAS beta coef.rds"))
