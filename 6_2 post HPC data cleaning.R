# Mixtures Analysis 
library(purrr)
library(tidyverse)
# Read in data that will be loaded on the HPC 

# Read and restructure results from mixtures analysis performed on HPC

source(here::here("!directories.r"))
source(here::here("!set_exposure_outcome_vars.r"))
source(here::here("!load_data.r"))
source(here::here("!functions.r"))

# List all files from HPC -----------------------------
temp = list.files(path = fs::path(dir_results_mixtures, 
                                  "from_hpc"), 
                  pattern="pw_pcs", 
                  full.names = TRUE)

# Get names
cohort_name <- case_when(str_detect(temp, "chs") ~ "CHS", 
                         str_detect(temp, "SOLAR") ~ "SOL")

# Read in all results --------------------------------------
results <- map(temp, 
               ~read_data_from_hpc(., n_col_in_df = 7)) 

# Rename lists 
names(results) <- cohort_name

# list to dataframe
results_df <- bind_rows(results, .id = "cohort") %>% 
  as_tibble()

# Clean results -------------------------------------
# Change "na" to missing 
results_df <- results_df %>% 
  mutate(p_value = na_if(p_value, "NA") %>% 
           as.numeric)



## Subset eta ------------------------------------
mixtures_eta <- filter(results_df,
                       term == "eta") %>% 
  rename(feature = metabolite)

## Pivot results wider on pips, ucl, lcl, etc  ---------------------------------
results_df_w <- results_df %>% 
  filter(term != "eta") %>%
  pivot_wider(id_cols = c("cohort", "metabolite", "exposure"), 
              names_from = term, 
              values_from = c(estimate, sd, lcl, ucl, p_value)) %>% 
  dplyr::select(cohort, metabolite, exposure, contains("beta"), contains("pip")) %>% 
  dplyr::select(-p_value_pip) %>% 
  dplyr::rename(p_value = p_value_beta)


# Calculate new p values
results_df_w <- results_df_w %>% 
  mutate(wald = (estimate_beta/sd_beta), 
         # p = 2*(1-pnorm(abs(wald),0,1)),
         p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
         p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
         neg_log_p = -log10(p_value)) %>% 
  group_by(cohort, exposure) %>% 
  mutate(q_value = p.adjust(p_value, method = "fdr"), 
         significance = if_else(p_value < 0.05, "p < 0.05", "Not Sig."), 
         significancefdr = if_else(q_value < 0.05, "q < 0.05", "Not Sig.")) %>% 
  ungroup()



# # Save final_results ---------------------------
# write_csv(results_final, 
#           file = fs::path(
#             dir_results_mixtures, ))
# 
# 
# 
# 
# # Get wide data frames of betas and pips --------------------------
# # Betas
# mwas_betas <- results_final %>% 
#   select(cohort, mode, exposure, feature, estimate_beta) %>% 
#   pivot_wider(id_cols = c(cohort, mode, feature), 
#               names_from = exposure, 
#               values_from = estimate_beta) 
# 
# # PIPs
# mwas_pips <- results_final   %>% 
#   select(cohort, mode, exposure, feature, estimate_pip) %>% 
#   pivot_wider(id_cols = c(cohort, mode, feature), 
#               names_from = exposure, 
#               values_from = estimate_pip) %>% 
#   select(-mixture)
# 
# # Make list of MWAS results
# mwas_results_long <- list(solar = results_final %>% filter(cohort == "SOL"),
#                           chs   = results_final %>% filter(cohort == "CHS"))
# 
# # Save MWAS results ------------------------------------------
# write_rds(mwas_results_long, 
#           fs::path(dir_results_mixtures, 
#                    ""))
# 
# 
# # Save MWAS beta coef and pip results -------------------------------
# mwas_beta_coefs <- list(solar = mwas_betas %>% filter(cohort == "SOL"),
#                         chs   = mwas_betas %>% filter(cohort == "CHS"), 
#                         solar_pips = mwas_pips %>% filter(cohort == "SOL"),
#                         chs_pips   = mwas_pips %>% filter(cohort == "CHS"))
# 
# write_rds(mwas_beta_coefs, 
#           fs::path(dir_results_mixtures, 
#                    ""))
# 
