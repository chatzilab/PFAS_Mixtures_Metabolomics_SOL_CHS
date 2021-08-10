library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)
library(fs)
ggplot2::theme_set(cowplot::theme_cowplot())

i = 1
# Get folder Directory
folders <- list.files(here::here("Temporary results", exposure_type))

folders <- folders[folders %in% exposures]
#Key for superpathways
key = readxl::read_xlsx(here::here("Supporting files",  
                                   "superpathway_key.xlsx"))

# 3) SOLAR: Load HILIC Mummichog data ------------------------------------------------
SOLAR_hilic_set <- read_rds(here::here("Temporary results", 
                                       exposure_type,
                                       folders[i], 
                                       "SOLAR",
                                       "hilic", 
                                       paste0(folders[i], "_hilic_mumichog_results.rds")))
# Get MZ and retention time key
mzkey_hilic <- SOLAR_hilic_set$mz2cpd_dict

# Mutate from list to DF
mzkey_hilic <- map2_dfr(mzkey_hilic, 
                        names(mzkey_hilic),
                        ~tibble(mz = .y, 
                                cpd = .x))

# 4) SOLAR: Map HILIC compounds to pathways -----------------------------------------------
# # get table of pathway p values
output <- as.data.frame(SOLAR_hilic_set$mummi.resmat) %>%
  rownames_to_column(., var = "name") %>%
  janitor::clean_names()

# Get name of paths
paths <- tibble(name = SOLAR_hilic_set$path.nms,
                path.hits = SOLAR_hilic_set$path.hits)

# Merge name of pathways and p values from output
output1 <- inner_join(paths, output) %>%  
  filter(pathway_total > 1)

# Merge compounds and pathways 
output2 <- map2_dfr(output1$path.hits, 
                    output1$name, 
                    ~tibble(name = .y, 
                            cpd = .x)) %>% 
  mutate(hit = 1)

###  Pivot compounds/pathways wider 
hilic_cpds_to_pathways <- pivot_wider(output2, 
                                      id_cols = cpd, 
                                      values_from = hit, 
                                      names_from = name)

# Merge compounds/pathways with superpathways
output2_superpathways <- left_join(output2, 
                                   key %>% dplyr::rename(name = pathway))

# Get number of pathways withing superpathways for each mz
output3_superpathways <- output2_superpathways %>% 
  group_by(cpd, super_pathway) %>% 
  summarise(n_pathways_within_superpatwhay = sum(hit), 
            pathways = paste(unique(name), collapse = "; ")) %>% 
  ungroup()

# Concatenate superpathways 
cpds_to_superpathways <- output3_superpathways %>% 
  group_by(cpd) %>% 
  summarise(n_superpathways = length(super_pathway), 
            pathways = pathways[1],
            main_sp_pthw = 
              if_else(n_superpathways == 1, 
                      super_pathway[1], 
                      paste(super_pathway, collapse = "; ")), 
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
              TRUE ~ "Other" )) %>% 
  ungroup() %>% 
  mutate(main_sp_pthw_categorized = 
           fct_lump_min(main_sp_pthw, 
                        10, 
                        other_level = "Other Combination of pathways"))

###  Join wide cpd/superpathway data with mzkey
solar_mzkey_hilic_w_superpathways <- inner_join(mzkey_hilic, 
                                                cpds_to_superpathways)




# 7) CHS: Load HILIC Mummichog data ------------------------------------------------
CHS_hilic_set <- read_rds(here::here("Temporary results", exposure_type,
                                     folders[2], 
                                     "CHS",
                                     "hilic", 
                                     paste0(folders[2],
                                            "_hilic_mumichog_results.rds")))
# Get MZ and retention time key
mzkey_hilic <- CHS_hilic_set$mz2cpd_dict

# Mutate from list to DF
mzkey_hilic <- map2_dfr(mzkey_hilic, 
                        names(mzkey_hilic),
                        ~tibble(mz = .y, 
                                cpd = .x))

# 8) CHS: Map HILIC compounds to pathways -----------------------------------------------
# # get table of pathway p values
output <- as.data.frame(CHS_hilic_set$mummi.resmat) %>%
  rownames_to_column(., var = "name") %>%
  janitor::clean_names()

# Get name of paths
paths <- tibble(name = CHS_hilic_set$path.nms,
                path.hits = CHS_hilic_set$path.hits)

# Merge name of pathways and p values from output
output1 <- inner_join(paths, output) %>%  
  filter(pathway_total > 1)

# Merge compounds and pathways 
output2 <- map2_dfr(output1$path.hits, 
                    output1$name, 
                    ~tibble(name = .y, 
                            cpd = .x)) %>% 
  mutate(hit = 1)

###  Pivot compounds/pathways wider 
hilic_cpds_to_pathways <- pivot_wider(output2, 
                                      id_cols = cpd, 
                                      values_from = hit, 
                                      names_from = name)

# Merge compounds/pathways with superpathways
output2_superpathways <- left_join(output2, 
                                   key %>% rename(name = pathway))

# Get number of pathways withing superpathways for each mz
output3_superpathways <- output2_superpathways %>% 
  group_by(cpd, super_pathway) %>% 
  summarise(n_pathways_within_superpatwhay = sum(hit), 
            pathways = paste(unique(name), collapse = "; ")) %>% 
  ungroup()

# Concatenate superpathways 
cpds_to_superpathways <- output3_superpathways %>% 
  group_by(cpd) %>% 
  summarise(n_superpathways = length(super_pathway),  
            pathways = pathways[1],
            main_sp_pthw = if_else(n_superpathways == 1, 
                                   super_pathway[1], 
                                   paste(super_pathway, collapse = "; "))) %>% 
  ungroup() %>% 
  mutate(main_sp_pthw_categorized = fct_lump_min(main_sp_pthw, 
                                                 10, 
                                                 other_level = "Other Combination of pathways"))

#  Join wide cpd/superpathway data with mzkey
chs_mzkey_hilic_w_superpathways <- inner_join(mzkey_hilic, 
                                              cpds_to_superpathways)



# 9) Combine SOLAR and CHS keys, match cpds and mzs -------------------

solar_mzkey_hilic_w_superpathways$mode <- "hilic"

solar_mzkey <- solar_mzkey_hilic_w_superpathways

# Merge mzkeys from C18 & HILIC, CHS
# chs_mzkey_c18_w_superpathways$mode <- "c18"
chs_mzkey_hilic_w_superpathways$mode <- "hilic"

chs_mzkey <- chs_mzkey_hilic_w_superpathways


mz_cpd_df1 <- full_join(solar_mzkey, 
                        chs_mzkey, 
                        by = c("cpd", "mode"), 
                        suffix = c("_solar", "_chs"))

mz_cpd_df1 <- mz_cpd_df1 %>% 
  mutate(across(contains("mz_"), as.numeric))


mz_cpd_df2 <- mz_cpd_df1 %>% 
  mutate(delta_mz = abs(mz_solar-mz_chs)) %>% 
  filter(delta_mz < 0.001) %>% 
  select(cpd, mode, aa_cho_lipid, contains("solar"), contains("chs"))


# 10) Save SOLAR and CHS results separately -------------------------------

solar_mult_mz_mult_cpd <- mz_cpd_df2 %>% 
  select(cpd, mode, aa_cho_lipid, contains("solar")) %>% 
  rename_all(~str_remove(.x, "_solar")) %>% 
  write_rds(here::here("Temporary results", exposure_type,
                       "cpd_mz_matching", 
                       "solar_mult_mz_mult_cpd.RDS"))


chs_mult_mz_mult_cpd <- mz_cpd_df2 %>% 
  select(cpd, mode,aa_cho_lipid, contains("chs")) %>% 
  rename_all(~str_remove(.x, "_chs")) %>% 
  write_rds(here::here("Temporary results", exposure_type,
                       "cpd_mz_matching", 
                       "chs_mult_mz_mult_cpd.RDS"))

