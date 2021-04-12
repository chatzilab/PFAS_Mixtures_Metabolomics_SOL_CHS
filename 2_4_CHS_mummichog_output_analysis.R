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

# Clean MWAS data  --------------------------------------------------------
#  Set working directory for first folder, input data
mwas_results_c18 <- read_rds(here::here("Temporary results", "1.1.0_chs_c18_mwas.rds"))
mwas_results_hilic <- read_rds(here::here("Temporary results", "1.1.0_chs_hilic_mwas.rds"))
mwas_results <- bind_rows(mwas_results_c18, 
                          mwas_results_hilic) %>% 
  janitor::clean_names() %>% 
  rename(mz = m_z, rt = r_t)

# Get min p value for each feature
mwas_min_p_vals <- mwas_results %>% 
  group_by(mode, mz, rt) %>% 
  filter(p_value == min(p_value)) %>% 
  ungroup() %>% 
  rename(min_p_value = p_value, 
         max_t_score = t_score)

# Get significant features with any exposure C18
c18_mwas_min_p_vals <- mwas_min_p_vals %>% 
  filter(mode == "negative") %>% 
  mutate(exposure = fct_infreq(exposure))


# Get significant features with any exposure hilic
hilic_mwas_min_p_vals <- mwas_min_p_vals %>% 
  filter(mode == "positive") %>% 
  mutate(exposure = fct_infreq(exposure))

# C18: Load Mummichog data ------------------------------------------------
chs_c18_set <- read_rds(here::here("Temporary results",
                                   folders[i], "chs",
                                   "c18", 
                                   paste0(folders[i],
                                          "_c18_mumichog_results.rds")))
# Get MZ and retention time key
mzkey_c18 <- chs_c18_set$mz2cpd_dict

# Mutate from list to DF
mzkey_c18 <- map2_dfr(mzkey_c18, 
                      names(mzkey_c18),
                      ~tibble(mz = .y, 
                              cpd = .x))


# C18: Map compounds to pathways -----------------------------------------------
# # get table of pathway p values
output <- as.data.frame(chs_c18_set$mummi.resmat) %>%
  rownames_to_column(., var = "name") %>%
  janitor::clean_names()

# Get name of paths
paths <- tibble(name = chs_c18_set$path.nms,
                path.hits = chs_c18_set$path.hits)

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
c18_cpds_to_pathways <- pivot_wider(output2, 
                                    id_cols = cpd, 
                                    values_from = hit, 
                                    names_from = name)

###  Join wide cpd/pathway data with mzkey
# mzkey_c18_w_pathways <- inner_join(mzkey_c18, 
#                                    cpds_to_pathways)

# Merge compounds/pathways with superpathways
output2_superpathways <- left_join(output2, 
                                   key %>% rename(name = pathway))

# Get number of pathways withing superpathways for each mz
output3_superpathways <- output2_superpathways %>% 
  select(-name) %>% 
  group_by(cpd, super_pathway) %>% 
  summarise(n_pathways_within_superpatwhay = sum(hit)) %>% 
  ungroup()

# Concatenate superpathways 
cpds_to_superpathways <- output3_superpathways %>% 
  group_by(cpd) %>% 
  summarise(n_superpathways = length(super_pathway),  
            main_super_pathway = if_else(n_superpathways == 1, 
                                         super_pathway[1], 
                                         paste(super_pathway, collapse = ", "))) %>% 
  ungroup() %>% 
  mutate(main_super_pathway_categorized = fct_lump_min(main_super_pathway, 
                                                       10, 
                                                       other_level = "Other Combination of pathways"))

###  Join wide cpd/superpathway data with mzkey
mzkey_c18_w_superpathways <- inner_join(mzkey_c18, 
                                        cpds_to_superpathways)


# C18: Join significant features with mzkey from mummichog ------------------------------------------------
x1_c18 <- left_join(c18_mwas_min_p_vals %>% 
                      mutate(mz = as.character(mz)), 
                    mzkey_c18_w_superpathways) %>% 
  filter(!is.na(cpd))


# Select only the most significant mz from each cpd 
c18_1cpd_mult_mz <- x1_c18 %>% 
  group_by(cpd) %>% 
  filter(min_p_value == min(min_p_value))

# Filter only single mz/cpd matches
c18_single_matching_only <- c18_1cpd_mult_mz %>% 
  group_by(mz) %>% 
  summarise(n_matches = length(cpd), 
            across(everything(), ~.x[1])) %>% 
  filter(n_matches == 1) %>% 
  select(mz, rt, cpd, main_super_pathway, main_super_pathway_categorized) %>% 
  mutate(mode = "negative")


write_csv(c18_single_matching_only, 
          here::here("Temporary results", 
                     "chs_C18_mz_cpd_superpathway_key.csv"))

# HILIC: Load Mummichog data ------------------------------------------------
chs_hilic_set <- read_rds(here::here("Temporary results",
                                     folders[i], "chs",
                                     "hilic", 
                                     paste0(folders[i],
                                            "_hilic_mumichog_results.rds")))
# Get MZ and retention time key
mzkey_hilic <- chs_hilic_set$mz2cpd_dict

# Mutate from list to DF
mzkey_hilic <- map2_dfr(mzkey_hilic, 
                        names(mzkey_hilic),
                        ~tibble(mz = .y, 
                                cpd = .x))

# Map compounds to pathways -----------------------------------------------
# # get table of pathway p values
output <- as.data.frame(chs_hilic_set$mummi.resmat) %>%
  rownames_to_column(., var = "name") %>%
  janitor::clean_names()

# Get name of paths
paths <- tibble(name = chs_hilic_set$path.nms,
                path.hits = chs_hilic_set$path.hits)

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

###  Join wide cpd/pathway data with mzkey
# mzkey_hilic_w_pathways <- inner_join(mzkey_hilic, 
#                                    cpds_to_pathways)


# Merge compounds/pathways with superpathways
output2_superpathways <- left_join(output2, 
                                   key %>% rename(name = pathway))

# Get number of pathways withing superpathways for each mz
output3_superpathways <- output2_superpathways %>% 
  select(-name) %>% 
  group_by(cpd, super_pathway) %>% 
  summarise(n_pathways_within_superpatwhay = sum(hit)) %>% 
  ungroup()

# Concatenate superpathways 
cpds_to_superpathways <- output3_superpathways %>% 
  group_by(cpd) %>% 
  summarise(n_superpathways = length(super_pathway),  
            main_super_pathway = if_else(n_superpathways == 1, 
                                         super_pathway[1], 
                                         paste(super_pathway, collapse = ", "))) %>% 
  ungroup() %>% 
  mutate(main_super_pathway_categorized = fct_lump_min(main_super_pathway, 
                                                       10, 
                                                       other_level = "Other Combination of pathways"))

###  Join wide cpd/superpathway data with mzkey
mzkey_hilic_w_superpathways <- inner_join(mzkey_hilic, 
                                          cpds_to_superpathways)


# Join significant features hilic with mzkey from mummichog ------------------------------------------------
x1_hilic <- left_join(hilic_mwas_min_p_vals %>% 
                        mutate(mz = as.character(mz)), 
                      mzkey_hilic_w_superpathways) %>% 
  filter(!is.na(cpd))


# Select only the most significant mz from each cpd 
hilic_1cpd_mult_mz <- x1_hilic %>% 
  group_by(cpd) %>% 
  filter(min_p_value == min(min_p_value))

# Filter only single mz/cpd matches
hilic_single_matching_only <- hilic_1cpd_mult_mz %>% 
  group_by(mz) %>% 
  summarise(n_matches = length(cpd), 
            across(everything(), ~.x[1])) %>% 
  filter(n_matches == 1) %>% 
  select(mz, rt, cpd, main_super_pathway, main_super_pathway_categorized) %>% 
  mutate(mode = "positive")


# 
write_csv(hilic_single_matching_only, 
          here::here("Temporary results", 
                     "chs_hilic_mz_cpd_superpathway_key.csv"))



# Combine with MWAS results  ----------------------------------------------
chs_single_matches <- bind_rows(hilic_single_matching_only, 
                                c18_single_matching_only)


chs_single_matches_w_mwas <- left_join(chs_single_matches, 
                                       mwas_results %>% 
                                         mutate(mz = as.character(mz))) %>% 
  mutate(across(c(mz, rt, p_value, t_score, beta), as.numeric ))


chs_single_matches_w_mwas <- chs_single_matches_w_mwas %>% 
  droplevels() %>% 
  mutate(super_pathway = fct_lump(main_super_pathway, 7,
                                  other_level = "Multiple/other") %>% 
           replace_na("Multiple/other") %>% 
           fct_infreq()) %>% 
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
                             paste0(exposure2, "*"), exposure2), 
  ) %>% 
  ungroup()

#########################
ggplot(single_matches_w_mwas,
       aes(x = cpdnum,
           y = -log(p_value),
           color = super_pathway)) +
  geom_hline(yintercept = -log(0.05/532),
             color = "grey80",
             size = .5, linetype =2) +
  geom_point(size = .75, shape = 21) +
  # geom_text(aes(label = cpd_sig)) + 
  facet_wrap(~exposure2) + 
  xlab("Compound Number (ordered by p-value within superpathway)") + 
  ylab("-log P")



ggplot(chs_single_matches_w_mwas,
       aes(x = beta,
           y = -log(p_value),
           color = super_pathway)) +
  # geom_hline(yintercept = -log(0.05),
  #            color = "grey80",
  #            size = .5, linetype =2) +
  geom_point(size = .75, shape = 21) +
  # geom_text(aes(label = cpd_sig)) + 
  facet_wrap(~exposure2)
xlab("Compound Number (ordered by p-value within superpathway)") + 
  ylab("-log P")




hist(chs_single_matches_w_mwas$beta, breaks = 100)