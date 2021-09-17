library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)
library(ggrepel)
library(fs)
ggplot2::theme_set(cowplot::theme_cowplot())
source(here::here("0_0_1_format_vars_funs.R"))

# Set important vars
exposure_type = "PFAS"

exposures = c("netfosaa","nmefosaab","pfbs","pfda","pfdoa","pfds","pfhpa",
              "pfhps","pfhxa","pfhxs","pfna","pfns","pfoa","pfos","pfpes",
              "pfuda", "x82fts")
cohort = c("solar", "chs")
modes = c("c18pos","c18neg", "hilicpos", "hilicneg")

# Temp Results folder architecture:
## Exposure Type > Exposure Name > Cohort > modes

#Key for superpathways
key <- readxl::read_xlsx(fs::path(dir_data, 
                                  "Supporting files",  
                                  "superpathway_key_sept_21.xlsx")) %>% 
  rename(path = pathway)

# Get list of all results folders ------------------------
dir_temp_exposures <- fs::path(dir_temp, exposure_type, exposures) 

dir_temp_exposures_chrt_mode <- map(dir_temp_exposures,
                                    ~fs::path(.x, cohort)) %>% 
  unlist() %>% 
  map(., ~fs::path(.x, modes)) %>% 
  unlist()

# Remove folder locations where mum. was not run.
dir_temp_exposures_chrt_mode <- dir_temp_exposures_chrt_mode[
  file_exists(fs::path(dir_temp_exposures_chrt_mode, 
                       "mummichog_pathway_enrichment.csv"))]

# Get list of all mummichog files
mum_rds_files = dir(dir_temp_exposures_chrt_mode, 
                    recursive=TRUE, 
                    full.names=TRUE,
                    pattern="\\.RDS$")

# 0) Load Mummichog RDS files --------------------------------------------------
# mum_results <- map(mum_rds_files, read_rds)
# 
# # Names of files 
# mum_rds_files_base <- basename(mum_rds_files) %>% 
#   str_remove(., "_mummichog.RDS")
# 
# # Set Names
# names(mum_results) <-  mum_rds_files_base

# 1) Load Mummichog pathway results --------------------------------------------
mum_pw <- read_csv(fs::path(dir_temp_exposures_chrt_mode, 
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
wgt_sol =sqrt(312)
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
  left_join(key)


# Combine pathways with long data to get list to include
mum_pw_final <- tidylog::left_join(mum_pw2, 
                                   mum_pw_final_w %>% 
                                     select(pfas, path, path_2, sig,
                                            fet_meta, enrichment_meta), 
                                   by = c("pfas", "path", "path_2")) %>% 
  tidylog::left_join(key) %>% 
  filter(!is.na(sig))



# Save Data 
write_rds(mum_pw_final_w, 
          here::here("Temporary results", 
                     exposure_type, 
                     "mum_pathway_results", 
                     "SOL CHS PFAS Mummichog wide sig PW.RDS"))


write_rds(mum_pw_final, 
          here::here("Temporary results", 
                     exposure_type, 
                     "mum_pathway_results", 
                     "SOL CHS PFAS Mummichog long sig PW.RDS"))


# 2) Get Empirical compound to Pathway Key -------------------------------------

# First, you have to group by mode to get mode specific EC's'
# ec_pwk1 <- mum_pw_lst %>% 
#   modify(~filter(.x, pathway_total > 3) %>% 
#            )
  



# 3) Read in mz/rt data for all modes ------------------------------------------
mz_rt_key_files <- fs::path(dir_temp, exposure_type, exposures[1], cohort[1], modes, 
                            "mummichog_matched_compound_all.csv")

# Get mz and retention time key
mz_rt_key = lapply(mz_rt_key_files, read_csv) 
names(mz_rt_key) <- modes


# 4) Create mz_rt to empirical compound key -----------------------------------------------
# mz_rt_key to data frame
mzrtkey1 <- bind_rows(mz_rt_key, .id = "mode") %>%
  janitor::clean_names()

# Set primary ions
primary_ions <- c("M+H[1+]",
                  "M+Na[1+]",
                  "M-H2O+H[1+]",
                  "M-H[-]",
                  "M-2H[2-]",
                  "M-H2O-H[-]",
                  "M+H[1+]",
                  "M+Na[1+]",
                  "M-H2O+H[1+]",
                  "M-H[-]",
                  "M-2H[2-]",
                  "M-H2O-H[-]")

# Check to make sure primary ions are named correctly
length(unique(mzrtkey1$empirical_compound))


table(mzrtkey1$matched_form)


# Filter only compound matches which were
mzrtkey2 <- mzrtkey1 %>%
  filter(matched_form %in% primary_ions) %>%
  mutate(mode = as.factor(mode), 
         feature = str_c(query_mass, retention_time, sep = "_"))


# length(unique(mzrtkey2$empirical_compound))


# MXrt key list from mummichog
mzrtkey_lst <- split(mzrtkey2, f = mzrtkey2$mode)









# Read in MWAS Beta Coefficients --------------------
mwas_beta_coefs <- read_rds(
  here::here("Temporary results", 
             "PFAS", 
             "SOL CHS all MWAS beta coef.rds"))

#
#
#
#






# First: 

# Get all MWAS results

# Get file linking ECs and pathways


































# 
# # OLD: ------------------------------------
# 
# 
# # # get table of pathway p values
# solar_pfhxs_c18neg_mummichog <- mum_pw_lst$solar_pfhxs_c18neg
# output <- mum_pw_lst$solar_pfhxs_c18neg
# 
# # Get pathway names
# paths <- tibble(path = solar_pfhxs_c18neg_mummichog$path,
#                 path.hits = solar_pfhxs_c18neg_mummichog$ec_hits)
# 
# # Merge name of pathways and p values from output
# output1 <- inner_join(paths, output) %>%  
#   filter(pathway_total > 1)
# 
# # Merge compounds and pathways 
# output2 <- map2_dfr(output1$path.hits, 
#                     output1$path, 
#                     ~tibble(name = .y, 
#                             cpd = .x)) %>% 
#   mutate(hit = 1)
# 
# ###  Pivot compounds/pathways wider 
# hilic_cpds_to_pathways <- pivot_wider(output2, 
#                                       id_cols = cpd, 
#                                       values_from = hit, 
#                                       names_from = name)
# 
# # Merge compounds/pathways with superpathways
# output2_superpathways <- left_join(output2, 
#                                    key %>% dplyr::rename(name = pathway))
# 
# # Get number of pathways withing superpathways for each mz
# output3_superpathways <- output2_superpathways %>% 
#   group_by(cpd, super_pathway) %>% 
#   summarise(n_pathways_within_superpatwhay = sum(hit), 
#             pathways = paste(unique(name), collapse = "; ")) %>% 
#   ungroup()
# 
# # Concatenate superpathways 
# cpds_to_superpathways <- output3_superpathways %>% 
#   group_by(cpd) %>% 
#   summarise(n_superpathways = length(super_pathway), 
#             pathways = pathways[1],
#             main_sp_pthw = 
#               if_else(n_superpathways == 1, 
#                       super_pathway[1], 
#                       paste(super_pathway, collapse = "; ")), 
#             aa_cho_lipid = case_when(
#               str_detect(main_sp_pthw, "Amino acid") & 
#                 str_detect(main_sp_pthw, "Lipid")    &
#                 str_detect(main_sp_pthw, "Carbohydrate") ~ "Amino Acid, Lipid, & CHO metab.", 
#               str_detect(main_sp_pthw, "Amino acid") & 
#                 str_detect(main_sp_pthw, "Lipid")        ~ "Amino Acid & Lipid metab.", 
#               str_detect(main_sp_pthw, "Amino acid") & 
#                 str_detect(main_sp_pthw, "Carbohydrate") ~ "Amino Acid & CHO metab.",
#               str_detect(main_sp_pthw, "Amino acid")     ~ "Amino Acid metab.",
#               str_detect(main_sp_pthw, "Lipid")          ~ "Lipid metab.",
#               str_detect(main_sp_pthw, "Carbohydrate")   ~ "CHO metab.",
#               str_detect(main_sp_pthw, "Nucleotide")   ~ "Nucleotide metab.",
#               str_detect(main_sp_pthw, "Glycan")       ~ "Glycan metab.",
#               str_detect(main_sp_pthw, "vitamins")     ~ "Metab. of cofactors and vit.",
#               str_detect(main_sp_pthw, "Xenobiotics")  ~ "Xenobiotics biodegradation and metab.",
#               str_detect(main_sp_pthw, "Steroid hormone")  ~ "Steroid hormone metab.",
#               TRUE ~ "Other" )) %>% 
#   ungroup() %>% 
#   mutate(main_sp_pthw_categorized = 
#            fct_lump_min(main_sp_pthw, 
#                         10, 
#                         other_level = "Other Combination of pathways"))
# 
# ###  Join wide cpd/superpathway data with mzkey
# solar_mzkey_hilic_w_superpathways <- inner_join(mzkey_hilic, 
#                                                 cpds_to_superpathways)
# 
# 
# 
# 
# 
# 
# 
# 
# # 7) CHS: Load HILIC Mummichog data ------------------------------------------------
# CHS_hilic_set <- read_rds(here::here("Temporary results", exposure_type,
#                                      folders[2], 
#                                      "CHS",
#                                      "hilic", 
#                                      paste0(folders[2],
#                                             "_hilic_mumichog_results.rds")))
# # Get MZ and retention time key
# mzkey_hilic <- CHS_hilic_set$mz2cpd_dict
# 
# # Mutate from list to DF
# mzkey_hilic <- map2_dfr(mzkey_hilic, 
#                         names(mzkey_hilic),
#                         ~tibble(mz = .y, 
#                                 cpd = .x))
# 
# # 8) CHS: Map HILIC compounds to pathways -----------------------------------------------
# # # get table of pathway p values
# output <- as.data.frame(CHS_hilic_set$mummi.resmat) %>%
#   rownames_to_column(., var = "name") %>%
#   janitor::clean_names()
# 
# # Get name of paths
# paths <- tibble(name = CHS_hilic_set$path.nms,
#                 path.hits = CHS_hilic_set$path.hits)
# 
# # Merge name of pathways and p values from output
# output1 <- inner_join(paths, output) %>%  
#   filter(pathway_total > 1)
# 
# # Merge compounds and pathways 
# output2 <- map2_dfr(output1$path.hits, 
#                     output1$name, 
#                     ~tibble(name = .y, 
#                             cpd = .x)) %>% 
#   mutate(hit = 1)
# 
# ###  Pivot compounds/pathways wider 
# hilic_cpds_to_pathways <- pivot_wider(output2, 
#                                       id_cols = cpd, 
#                                       values_from = hit, 
#                                       names_from = name)
# 
# # Merge compounds/pathways with superpathways
# output2_superpathways <- left_join(output2, 
#                                    key %>% rename(name = pathway))
# 
# # Get number of pathways withing superpathways for each mz
# output3_superpathways <- output2_superpathways %>% 
#   group_by(cpd, super_pathway) %>% 
#   summarise(n_pathways_within_superpatwhay = sum(hit), 
#             pathways = paste(unique(name), collapse = "; ")) %>% 
#   ungroup()
# 
# # Concatenate superpathways 
# cpds_to_superpathways <- output3_superpathways %>% 
#   group_by(cpd) %>% 
#   summarise(n_superpathways = length(super_pathway),  
#             pathways = pathways[1],
#             main_sp_pthw = if_else(n_superpathways == 1, 
#                                    super_pathway[1], 
#                                    paste(super_pathway, collapse = "; "))) %>% 
#   ungroup() %>% 
#   mutate(main_sp_pthw_categorized = fct_lump_min(main_sp_pthw, 
#                                                  10, 
#                                                  other_level = "Other Combination of pathways"))
# 
# #  Join wide cpd/superpathway data with mzkey
# chs_mzkey_hilic_w_superpathways <- inner_join(mzkey_hilic, 
#                                               cpds_to_superpathways)
# 
# 
# 
# # 9) Combine SOLAR and CHS keys, match cpds and mzs -------------------
# 
# solar_mzkey_hilic_w_superpathways$mode <- "hilic"
# 
# solar_mzkey <- solar_mzkey_hilic_w_superpathways
# 
# # Merge mzkeys from C18 & HILIC, CHS
# # chs_mzkey_c18_w_superpathways$mode <- "c18"
# chs_mzkey_hilic_w_superpathways$mode <- "hilic"
# 
# chs_mzkey <- chs_mzkey_hilic_w_superpathways
# 
# 
# mz_cpd_df1 <- full_join(solar_mzkey, 
#                         chs_mzkey, 
#                         by = c("cpd", "mode"), 
#                         suffix = c("_solar", "_chs"))
# 
# mz_cpd_df1 <- mz_cpd_df1 %>% 
#   mutate(across(contains("mz_"), as.numeric))
# 
# 
# mz_cpd_df2 <- mz_cpd_df1 %>% 
#   mutate(delta_mz = abs(mz_solar-mz_chs)) %>% 
#   filter(delta_mz < 0.001) %>% 
#   select(cpd, mode, aa_cho_lipid, contains("solar"), contains("chs"))
# 
# 
# # 10) Save SOLAR and CHS results separately -------------------------------
# 
# solar_mult_mz_mult_cpd <- mz_cpd_df2 %>% 
#   select(cpd, mode, aa_cho_lipid, contains("solar")) %>% 
#   rename_all(~str_remove(.x, "_solar")) %>% 
#   write_rds(here::here("Temporary results", exposure_type,
#                        "cpd_mz_matching", 
#                        "solar_mult_mz_mult_cpd.RDS"))
# 
# 
# chs_mult_mz_mult_cpd <- mz_cpd_df2 %>% 
#   select(cpd, mode,aa_cho_lipid, contains("chs")) %>% 
#   rename_all(~str_remove(.x, "_chs")) %>% 
#   write_rds(here::here("Temporary results", exposure_type,
#                        "cpd_mz_matching", 
#                        "chs_mult_mz_mult_cpd.RDS"))
# 
