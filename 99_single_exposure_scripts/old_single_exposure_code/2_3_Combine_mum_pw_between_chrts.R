library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)
library(ggrepel)
library(fs)
ggplot2::theme_set(cowplot::theme_cowplot())

# Source setup scripts
source(here::here("0_0_1_format_vars_funs.R"))
source(here::here("!directories.R"))
source(here::here("!set_exposure_outcome_vars.R"))

# Temp Results folder architecture:
## Exposure Type > Exposure Name > Cohort > modes

#Key for superpathways
superpathwaykey <- readxl::read_xlsx(fs::path(dir_data, 
                                  "Supporting files",  
                                  "superpathway_key_sept_21.xlsx")) %>% 
  rename(path = pathway)

# Get list of all results folders ------------------------
dir_results_exposures <- fs::path(dir_results, exposure_type, exposures) 


dir_results_exposures_chrt_mode <- map(dir_results_exposures,
                                    ~fs::path(.x, cohort)) %>% 
  unlist() %>% 
  map(., ~fs::path(.x, modes)) %>% 
  unlist()

# Remove folder locations where mum. was not run.
dir_results_exposures_chrt_mode <- dir_results_exposures_chrt_mode[
  file_exists(fs::path(dir_results_exposures_chrt_mode, 
                       "mummichog_pathway_enrichment.csv"))]

# Get list of all mummichog files
mum_rds_files = dir(dir_results_exposures_chrt_mode, 
                    recursive=TRUE, 
                    full.names=TRUE,
                    pattern="\\.RDS$")

# 0) Load Mummichog RDS files --------------------------------------------------
# mum_results <- map(mum_rds_files, read_rds)
# # 
# # # Names of files 
# mum_rds_files_base <- basename(mum_rds_files) %>%
#   str_remove(., "_mummichog.RDS")
# # 
# # # Set Names
# names(mum_results) <-  mum_rds_files_base
# write_rds(mum_results, 
#           fs::path(dir_results,
#                    exposure_type, 
#                    "mum_pathway_results", 
#                    "raw_mum_results_files.RDS"))



# 1) Load Mummichog pathway results --------------------------------------------
mum_pw <- read_csv(fs::path(dir_results_exposures_chrt_mode, 
                            "mummichog_pathway_enrichment.csv"), 
                   id = "file_name") %>% 
  janitor::clean_names() %>% 
  rename(path = x1) %>% 
  mutate(file_name = str_replace_all(file_name, "/", "_") %>% 
           str_remove("_mummichog_pathway_enrichment.csv"))


# Get columns for PFAS, cohort, and mode
mum_pw1 <- mum_pw %>% 
  mutate( 
    temp = str_split(file_name,  '_PFAS_') %>% 
      map_chr(2), 
    pfas = str_split(temp,  '_') %>% 
      map_chr(1),
    cohort = str_split(temp,  '_') %>% 
      map_chr(2), 
    mode = str_split(temp,  '_') %>% 
      map_chr(3), 
    enrichment = hits_sig/hits_total, 
    neg_logp = -log(fet),
    name = str_c(cohort, pfas, mode, sep = "_") %>% 
      tolower(),
    pfasog = pfas, 
    pfas = rename_pfas(pfas,include_asterisk = TRUE), 
    path_2  = str_replace(path, "metabolism", "met.") %>% 
      str_replace("Metabolism", "met.") %>% 
      str_replace(., " pathway", "")) %>% 
  select(pfas, cohort, mode, name, everything(), -temp, -file_name, -pathway_number)


# Change mumichog pathway results to list for easier analysis down the line
mum_pw_lst <- mum_pw1 %>%
  filter(pathway_total > 3) %>% 
  split(., f = .$name)


# Collapse across modes, within PFAS
mum_pw2 <- mum_pw1 %>% 
  filter(pathway_total > 3) %>% 
  group_by(pfas, path, cohort) %>% 
  filter(fet == min(fet)) %>% 
  ungroup()


# Pivot wider
mum_pw_w1 <- pivot_wider(mum_pw2, 
                         id_cols = c(pfas,  path, path_2), 
                         names_from = cohort, 
                         values_from = c(mode, pathway_total:neg_logp))


# perform meta analysis of p values to get "combined" column
wgt_sol = sqrt(310)
wgt_chs = sqrt(136)
mum_pw_w1 <- mum_pw_w1 %>% 
  mutate(fet_solar2 = if_else(is.na(fet_solar), .99, fet_solar), 
         fet_chs2 = if_else(is.na(fet_chs), .99, fet_chs)) %>%
  rowwise() %>% 
  mutate(fet_meta = metap::sumz(p = c_across(fet_solar2:fet_chs2), 
                                    weights = c(wgt_sol, wgt_chs))$p, 
         enrichment_meta = 
           ((enrichment_solar*wgt_sol)+(enrichment_chs*wgt_chs))/(wgt_sol+wgt_chs), 
         neg_logp_meta = -log(fet_meta), 
         sig_meta = if_else(fet_meta < 0.01, "Sig.", "Not Sig.")) %>% 
  ungroup()


# select only pathways which were reported in both cohorts
mum_pw_final_w <- mum_pw_w1 %>% 
  filter(!is.na(fet_solar), !is.na(fet_chs)) %>%
  mutate(sig = case_when(fet_solar<0.05 & fet_chs < 0.05 ~ "Sig. Both Cohorts", 
                         fet_solar<0.05 ~ "Sig. SOLAR Only", 
                         fet_chs  <0.05 ~ "Sig. CHS Only", 
                         TRUE ~ "Not Significant")) %>% 
  left_join(superpathwaykey)


# Combine pathways with long data to get list to include
mum_pw_final <- tidylog::left_join(mum_pw2, 
                                   mum_pw_final_w %>% 
                                     select(pfas, path, path_2, sig,
                                            fet_meta, enrichment_meta), 
                                   by = c("pfas", "path", "path_2")) %>% 
  tidylog::left_join(superpathwaykey) %>% 
  filter(!is.na(sig))


# Clean Environment
rm(mum_pw, mum_pw1, mum_pw2, mum_pw_w1, wgt_chs, wgt_sol, mum_rds_files)

# Save Data 
write_rds(mum_pw_final_w, 
          fs::path(dir_results, 
                     exposure_type, 
                     "mum_pathway_results", 
                     "SOL CHS PFAS Mummichog wide sig PW.RDS"))


write_rds(mum_pw_final, 
          fs::path(dir_results, 
                     exposure_type, 
                     "mum_pathway_results", 
                     "SOL CHS PFAS Mummichog long sig PW.RDS"))
