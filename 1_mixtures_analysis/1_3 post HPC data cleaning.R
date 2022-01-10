# Mixtures Analysis 
library(purrr)
library(tidyverse)
# Read in data that will be loaded on the HPC 

# Read and restructure results from mixtures analysis performed on HPC

# The issue: Each node writes results to a single row, 
# so data is not rectangular (because each node does not run exactly
# the same number of models, as 23173 is not divisible by 128 nodes). 
source(here::here("!directories.r"))
source(here::here("!set_exposure_outcome_vars.r"))
source(here::here("!load_data.r"))
source(here::here("!functions.r"))

# List all files from HPC -----------------------------
temp = list.files(path = fs::path(dir_results_mixtures, 
                                  "from_hpc"), 
                  pattern="all_pfas", 
                  full.names = TRUE)

# Get names
cohort_name <- case_when(str_detect(temp, "chs_all") ~ "CHS", 
                         str_detect(temp, "sol_all") ~ "SOL")

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


# Remove PFAS from results 
results_df <- results_df %>% 
  tidylog::filter(
    metabolite != "412.9662665_238.1521714",    # PFOA,  C18 neg
    metabolite != "499.9333381_260.9808557",    # PFOS,  C18 neg
    metabolite != "398.9360546_243.5838959",    # PFHxS, C18 neg
    metabolite != "398.936398885251_33.3622678197535", # PFHxS, HELIC Neg
    metabolite != "462.962944_247.5910964",     # PFNA
    metabolite != "448.9331458_253.4631484",    # PFHpS
    metabolite != "412.966235310664_35.0213619776204") # PFOA, HELIC Neg 



## Subset eta ------------------------------------
mixtures_eta <- filter(results_df,
                       term == "eta") %>% 
  rename(feature = metabolite)

## Pivot results wider ---------------------------------
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

# Join with ft_metadata
results_final <- tidylog::inner_join(ft_metadata, results_df_w, 
                                     by = c("feature" = "metabolite")) %>% 
  as_tibble()

# Save final_results ---------------------------
write_csv(results_final, 
          file = fs::path(
            dir_results_mixtures, 
            "all_pfas_mixtures_results_hyper_g.csv"))




# Get wide data frames of betas and pips --------------------------
# Betas
mwas_betas <- results_final %>% 
  select(cohort, mode, exposure, feature, estimate_beta) %>% 
  pivot_wider(id_cols = c(cohort, mode, feature), 
              names_from = exposure, 
              values_from = estimate_beta) 

# PIPs
mwas_pips <- results_final   %>% 
  select(cohort, mode, exposure, feature, estimate_pip) %>% 
  pivot_wider(id_cols = c(cohort, mode, feature), 
              names_from = exposure, 
              values_from = estimate_pip) %>% 
  select(-mixture)

# Make list of MWAS results
mwas_results_long <- list(solar = results_final %>% filter(cohort == "SOL"),
                          chs   = results_final %>% filter(cohort == "CHS"))

# Save MWAS results ------------------------------------------
write_rds(mwas_results_long, 
          fs::path(dir_results_mixtures, 
                   "SOL CHS all Mixtures MWAS results long hyper_g.rds"))


# Save MWAS beta coef and pip results -------------------------------
mwas_beta_coefs <- list(solar = mwas_betas %>% filter(cohort == "SOL"),
                        chs   = mwas_betas %>% filter(cohort == "CHS"), 
                        solar_pips = mwas_pips %>% filter(cohort == "SOL"),
                        chs_pips   = mwas_pips %>% filter(cohort == "CHS"))

write_rds(mwas_beta_coefs, 
          fs::path(dir_results_mixtures, 
                   "SOL CHS all Mixtures MWAS beta coef hyper_g.rds"))

