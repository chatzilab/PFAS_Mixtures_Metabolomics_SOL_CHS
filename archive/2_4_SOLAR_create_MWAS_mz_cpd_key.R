#################### SOLAR!!
library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)
library(fs)
ggplot2::theme_set(cowplot::theme_cowplot())
i = 1
# Get folder Directory
folders <- list.files(here("Temporary results"))

folders <- folders[folders %in% exposures]
#Key for superpathways
key = readxl::read_xlsx(here::here("Supporting files",  
                                   "superpathway_key.xlsx"))


# 1) SOLAR: Read back in MWAS data to create table of MWAS results ----------------------------------------
# Set working directory for first folder, input data

mwas_results_c18 <- read_rds(here::here("Temporary results", "1_1_SOLAR_c18_mwas.rds"))
mwas_results_hilic <- read_rds(here::here("Temporary results", "1_1_SOLAR_hilic_mwas.rds"))
mwas_results <- bind_rows(mwas_results_c18, 
                          mwas_results_hilic) %>% 
  janitor::clean_names() %>% 
  rename(mz = m_z, rt = r_t) %>% 
  mutate(across(mz:t_score, as.numeric))

# Get lowest p value for each feature
mwas_min_p_vals <- mwas_results %>% 
  group_by(mode, mz, rt) %>% 
  filter(p_value == min(p_value)) %>% 
  ungroup() %>% 
  rename(min_p_value = p_value, 
         max_t_score = t_score)


# 2) SOLAR: Read in cpd_mz key ----------------------------------------------

solar_mzkey_w_superpathways <- read_rds(here::here("Temporary results",
                     "cpd_mz_matching", 
                     "solar_mult_mz_mult_cpd.RDS")) %>% 
  select(-mode)

# 3) SOLAR: Join significant features with mzkey from mummichog ------------------------------------------------
x1 <- left_join(mwas_min_p_vals %>% 
                      mutate(mz = as.character(mz)), 
                    solar_mzkey_w_superpathways %>% 
                      mutate(mz = as.character(mz))) %>% 
  filter(!is.na(cpd))

# Select only the most significant mz from each cpd 
one_cpd_mult_mz <- x1 %>% 
  group_by(cpd) %>% 
  filter(min_p_value == min(min_p_value))


# Create unique mz/cpd matches
solar_single_matching_only <- one_cpd_mult_mz %>% 
  group_by(mz) %>% 
  summarise(n_matches = length(cpd), 
            cpd = paste(cpd, collapse = "; "),
            across(everything(), ~.x[1])) %>% 
  select(mz, rt, cpd, pathways, main_super_pathway, main_super_pathway_categorized, n_matches, min_p_value) %>% 
  mutate(mode = "negative")





# Combine with MWAS results  ----------------------------------------------

solar_single_matches_w_mwas <- left_join(solar_single_matching_only, 
                                         mwas_results %>% 
                                           mutate(mz = as.character(mz))) %>% 
  mutate(across(c(mz, rt, p_value, t_score, beta, conf_low, conf_high), 
                as.numeric ))


solar_single_matches_w_mwas <- solar_single_matches_w_mwas %>% 
  droplevels() %>% 
  mutate(super_pathway = fct_lump(main_super_pathway, 7,
                                  other_level = "Multiple/other") %>% 
           replace_na("Multiple/other") %>% 
           fct_infreq()) %>% 
  arrange(p_value) %>%
  arrange(super_pathway) %>% 
  ungroup()

solar_single_matches_w_mwas <- solar_single_matches_w_mwas %>% 
  group_by(exposure) %>%
  mutate(q_value = p.adjust(p_value, method = "BY")) %>% 
  mutate(cpdnum = row_number(), 
         cpd_sig = if_else(q_value < 0.05, cpd, ""), 
         exposure2 = str_remove_all(exposure, 
                                    "_impute") %>% 
           str_remove_all("_ngml_detect") %>% 
           toupper() %>% 
           str_replace_all("_", " ") %>% 
           str_replace("HEXACHLOROBENZENE", "HCB"),
         exposure2 = if_else(str_detect(exposure, "ngml_detect"), 
                             paste0(exposure2, "*"), exposure2)) %>% 
  ungroup()


write_rds(solar_single_matches_w_mwas, here::here("Temporary results", 
                                                  "solar_mz_cpd_pathway_key_w_mwas.rds"))
