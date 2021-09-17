#################### CHS!!
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


# 1) CHS: Read back in MWAS data to create table of MWAS results ----------------------------------------
# Set working directory for first folder, input data

mwas_results_c18 <- read_rds(fs::path(dir_temp, "1_1_CHS_c18_mwas.rds")) #%>% 
  # left_join(chs_c18_rsd)

mwas_results_hilic <- read_rds(fs::path(dir_temp, "1_1_CHS_hilic_mwas.rds"))  #%>% 
  # left_join(chs_hilic_rsd)

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


# 2) CHS: Read in cpd_mz key ----------------------------------------------
chs_mzkey_w_superpathways <- read_rds(fs::path(dir_temp,
                                                 "cpd_mz_matching", 
                                                 "chs_mult_mz_mult_cpd.RDS")) %>% 
  select(-mode)

# 3) CHS: Join significant features with mzkey from mummichog ------------------------------------------------
x1 <- left_join(mwas_min_p_vals %>% 
                  mutate(mz = as.character(mz)), 
                chs_mzkey_w_superpathways %>% 
                  mutate(mz = as.character(mz))) %>% 
  filter(!is.na(cpd)) %>% 
  select(-contains("conf_"), -c(exposure, max_t_score,
                                min_p_value, beta, mode)) %>% 
  mutate(pathways_2 = fct_lump(main_sp_pthw, 10,
                               other_level = "Multiple/other"))



x2 <-  x1 %>% 
  group_by(metabolite, mz, rt) %>% 
  summarise(cpd = str_c(cpd, collapse = "; "), 
            across(c(aa_cho_lipid, pathways, main_sp_pthw, main_sp_pthw_categorized), 
                   ~str_c(unique(.x), collapse = "; "))) %>% 
  ungroup()

  

# Combine with MWAS results  ----------------------------------------------

chs_single_matches_w_mwas <- left_join(x2, 
                                       mwas_results %>% 
                                         mutate(mz = as.character(mz))) %>% 
  mutate(across(c(mz, rt, p_value, t_score, beta, conf_low, conf_high), 
                as.numeric ))


# Create AA_Lipid_cho variable
chs_single_matches_w_mwas <- chs_single_matches_w_mwas %>% 
  droplevels() %>% 
  mutate(super_pathway = fct_lump(main_sp_pthw, 7,
                                  other_level = "Multiple/other") %>% 
           replace_na("Multiple/other") %>% 
           fct_infreq(), 
         aa_cho_lipid = case_when(
           str_detect(main_sp_pthw, "Amino acid") & 
             str_detect(main_sp_pthw, "Lipid")    &
             str_detect(main_sp_pthw, "Carbohydrate") ~ "Amino Acid, Lipid, & CHO metab.", 
           str_detect(main_sp_pthw, "Amino acid") & 
             str_detect(main_sp_pthw, "Lipid")        ~ "Amino Acid & Lipid metab.", 
           str_detect(main_sp_pthw, "Amino acid") & 
             str_detect(main_sp_pthw, "Carbohydrate") ~ "Amino Acid & CHO metab.",
           str_detect(main_sp_pthw, "Amino acid")     ~ "Amino Acid metab.",
           str_detect(main_sp_pthw, "Lipid")          ~ "Lipid metab.",
           str_detect(main_sp_pthw, "Carbohydrate")   ~ "CHO metab.",
           str_detect(main_sp_pthw, "Nucleotide")   ~ "Nucleotide metab.",
           str_detect(main_sp_pthw, "Glycan")       ~ "Glycan metab.",
           str_detect(main_sp_pthw, "vitamins")     ~ "Metab. of cofactors and vit.",
           str_detect(main_sp_pthw, "Xenobiotics")  ~ "Xenobiotics biodegradation and metab.",
           str_detect(main_sp_pthw, "Steroid hormone")  ~ "Steroid hormone metab.",
           TRUE ~ "Other")) %>% 
  arrange(p_value) %>%
  arrange(super_pathway) %>% 
  ungroup()



chs_single_matches_w_mwas <- chs_single_matches_w_mwas %>% 
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


write_rds(chs_single_matches_w_mwas, fs::path(dir_temp, 
                                                "chs_mz_cpd_pathway_key_w_mwas.rds"))


