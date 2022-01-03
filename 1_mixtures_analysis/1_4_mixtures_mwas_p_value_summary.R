# P-value Summary
library("jag2")
library("tidyverse")
source(here::here("!directories.r"))
# SOLAR --------------------------------------------------------
sol_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "sol_pfas_mixtures_results_w_09.csv")) %>% 
  as_tibble()


(solar_summary_mwas <- sol_mwas %>% 
   group_by(mode, exposure) %>% 
   summarise(n_features = length(feature), 
             percent_significant_p05 = jag2::npct(significance, 
                                            "p < 0.05", 
                                            n.digits = 2), 
             percent_significant_q05 = jag2::npct(significancefdr, 
                                            "q < 0.05", 
                                            n.digits = 2)))

# chs ----------------------------------------------------
chs_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "chs_pfas_mixtures_results_w_09.csv")) 

(chs_summary_mwas <- chs_mwas %>% 
    group_by(mode, exposure) %>% 
    summarise(n_features = length(feature), 
              percent_significant_p05 = npct(significance, 
                                             "p < 0.05", 
                                             n.digits = 2), 
              percent_significant_q05 = npct(significancefdr, 
                                             "q < 0.05", 
                                             n.digits = 2)))


# combine solar chs
mwas_results_summary <- full_join(solar_summary_mwas, 
                                  chs_summary_mwas, 
                                  by = c("mode", "exposure", "n_features"), 
                                  suffix = c("_solar", "_chs"))

write_excel_csv(mwas_results_summary, 
                file = fs::path(dir_reports, 
                                "Summary of Mixtures MWAS Results_w_09.csv"))


# Merge to get single data frame of all beta coefficients 
sol_mwas_betas <- sol_mwas   %>% 
  select(mode, exposure, feature, estimate_beta) %>% 
           pivot_wider(id_cols = c(mode, feature), 
                       names_from = exposure, 
                       values_from = estimate_beta) 


chs_mwas_betas <- chs_mwas  %>% 
  select(mode, exposure, feature, estimate_beta) %>% 
  pivot_wider(id_cols = c(mode, feature), 
              names_from = exposure, 
              values_from = estimate_beta) 



# Merge to get single data frame of all PIPs 
sol_mwas_pips <- sol_mwas   %>% 
  select(mode, exposure, feature, estimate_pip) %>% 
  pivot_wider(id_cols = c(mode, feature), 
              names_from = exposure, 
              values_from = estimate_pip ) 


chs_mwas_pips <- chs_mwas  %>% 
  select(mode, exposure, feature, estimate_pip ) %>% 
  pivot_wider(id_cols = c(mode, feature), 
              names_from = exposure, 
              values_from = estimate_pip ) 

# Make list of MWAS results
mwas_results_long <- list(solar = sol_mwas,
                          chs   = chs_mwas)

# 
write_rds(mwas_results_long, 
          fs::path(dir_results_mixtures, 
                   "SOL CHS all Mixtures MWAS results long w 09.rds"))

 
 
# # PCA analysis on the effect estimates
# library(ggfortify)
# temp <- chs_mwas_betas %>% select(-mode, -feature)
# pca_res <- prcomp(temp, scale. = TRUE)
# 
# autoplot(pca_res, data = chs_mwas_betas, colour = 'mode', alpha = .3)
# 


# Save MWAS beta coef and pip results 
mwas_beta_coefs <- list(solar = sol_mwas_betas, 
                        chs   = chs_mwas_betas, 
                        solar_pips = sol_mwas_pips,
                        chs_pips = chs_mwas_pips)

write_rds(mwas_beta_coefs, 
          fs::path(dir_results_mixtures, 
                   "SOL CHS all Mixtures MWAS beta coef w 09.rds"))