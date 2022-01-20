# Get most significant pathway from each PFAS
library(ggExtra)
library(tidylog)
# 2) Plot Mummichog Pathway Results
library(colorspace)
library(janitor)

# Format functions
source(here::here("0_0_1_format_vars_funs.R"))

# 1) read in data  ------------------------
## Mummichog Data ------------------------------------
mum_pw_wide <- read_rds(
  fs::path(dir_results_mum_mixtures,
           "Mixture effect hyper_g", 
           "SOL CHS PFAS Mummichog wide sig PW.RDS")) %>% 
  clean_names() %>% 
  mutate(q_meta = p.adjust(pval_meta), 
         sig_fdr = if_else(q_meta < 0.2, 
                           "Sig. FDR < 0.2",
                           "Not. Sig"), 
         mean_num_sig = (hits_sig_solar + hits_sig_chs)/2)

# Filter only sig pw in both cohorts
mum_pw_wide <- mum_pw_wide %>% 
  filter(combined_pvals_solar < 0.05, 
         combined_pvals_chs < 0.05)

## Read in mzrt key--------------------------------------
mzrt_key <- read_rds(
  fs::path(dir_data_local,
           "Supporting Files", 
           "mummichog_pw_ec_feature_key_with_cpd_names.rds")) %>% 
  filter(pathway %in% mum_pw_wide$path) %>% 
  select(pathway, met_name, everything())

## MWAS Data ------------------------------------
# Create full ft for both cohorts
ftdata2 <- ftdata %>%
  map(~purrr::reduce(., .f = left_join))

# select only mz/rt from one of the pathways of interest
ftdata_selected <- ftdata2 %>% 
  modify(~select(., id, all_of(mzrt_key$name)))


# 2) Get dataframe with column df's for all mzrt for a given PW --------
solar_data <- mzrt_key %>% 
  select(-cohort) %>% 
  distinct() %>%
  group_by(pathway) %>% 
  nest() %>% 
  mutate(ft = map(data, 
                  ~select(ftdata_selected$solar, 
                          id, all_of(.$name)) %>% 
                    column_to_rownames("id"))) %>% 
  rename(mzrt_key = data)


chs_data <- mzrt_key %>% 
  select(-cohort) %>% 
  distinct() %>%
  group_by(pathway) %>% 
  nest() %>% 
  mutate(ft = map(data, 
                  ~select(ftdata_selected$chs, 
                          id, all_of(.$name)) %>% 
                    column_to_rownames("id"))) %>% 
  rename(mzrt_key = data)


# 3) Perform PCA for each pathway --------------------------------------
# SOLAR
solar_data <- solar_data %>% 
  mutate(pc_res = map(ft, 
                      ~prcomp(x = ., center = TRUE, scale = TRUE)), 
         pc_12 = map(pc_res, 
                     ~as_tibble(.$x, rownames = "id") %>% 
                       select(id, PC1, PC2)), 
         pc_ft = map2(ft, pc_12, 
                      bind_cols),
         pc_ft = map(pc_ft, ~select(., id, contains("PC"))))

# CHS
chs_data <- chs_data %>% 
  mutate(pc_res = map(ft, 
                      ~prcomp(x = ., center = TRUE, scale = TRUE)), 
         pc_12 = map(pc_res, 
                     ~as_tibble(.$x, rownames = "id") %>% 
                       select(id, PC1, PC2)), 
         pc_ft = map2(ft, pc_12, 
                      bind_cols),
         pc_ft = map(pc_ft, ~select(., id, contains("PC"))))


# 4) Get first and second PC scores for each pathway --------------------------
#solar
solar_pathway_pcs <- solar_data %>% 
  select(pathway, pc_ft) %>% 
  unnest(pc_ft) %>% 
  pivot_wider(id_cols = "id", 
              names_from = "pathway",
              values_from = "PC1":"PC2") %>% 
  rename_all(tolower) %>% 
  clean_names()

#chs
chs_pathway_pcs <- chs_data %>% 
  select(pathway, pc_ft) %>% 
  unnest(pc_ft) %>% 
  pivot_wider(id_cols = "id", 
              names_from = "pathway",
              values_from = "PC1":"PC2") %>% 
  rename_all(tolower) %>% 
  clean_names()

# Join solar and chs into single list
pca_dat_per_path <- list(solar = solar_pathway_pcs, 
                         chs = chs_pathway_pcs)


# Save data 
write_rds(pca_dat_per_path, 
            fs::path(dir_results_mum_mixtures,
                     "SOL CHS PC Scores of sig pathways.RDS"))