# Get cpds from mum

library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)
library(ggrepel)
library(fs)
ggplot2::theme_set(cowplot::theme_cowplot())


# Source setup scripts
source(here::here("!directories.R"))
source(here::here("!set_exposure_outcome_vars.R"))
# source(here::here("2_3_Combine_mum_pw_between_chrts_mixtures.R"))


#1) Get Empirical compound to Pathway Key -------------------------------------

# First, you have to group by mode to get mode specific EC's'
# ec_pwk1 <- mum_pw_lst %>% 
#   modify(~filter(.x, pathway_total > 3) %>% 
#            )


#2) Read in mz/rt data for all modes ------------------------------------------
mz_rt_key_files <- fs::path(dir_results_mum_mixtures,
                            "Mixture effect w09", cohort[1], 
                            "mummichog_matched_compound_all.csv")

# Get mz and retention time key
mz_rt_key = lapply(mz_rt_key_files, read_csv) 
names(mz_rt_key) <- modes


# 3) Create mz_rt to empirical compound key ------------------------------------
# mz_rt_key to data frame
mzrtkey1 <- bind_rows(mz_rt_key, .id = "mode") %>%
  janitor::clean_names() %>% 
  mutate(ec_mode = str_c(empirical_compound, mode, sep = "_"))

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

# Filter only compound matches from prinary ions
mzrtkey_primary_ions <- mzrtkey1 %>%
  filter(matched_form %in% primary_ions) %>%
  mutate(mode = as.factor(mode), 
         feature = str_c(query_mass, retention_time, sep = "_"))

length(unique(mzrtkey_primary_ions$matched_compound))
# mzrt key list from mummichog
# mzrtkey_lst <- split(mzrtkey2, f = mzrtkey2$mode)

rm(mz_rt_key, mz_rt_key_files, mzrtkey1)

# Create pw to ec dataset -------------------------------------
# Read in raw mum_ results files
raw_mum_rslts <- read_rds(fs::path(dir_results_mum_mixtures,
                                   "mum_pathway_results_w_09", 
                                   "raw_mum_results_files.RDS"))

# #renamefxn
name_ecpathwayfxn <- function(x, nms){
  names(x) <- nms
  return(x)
}

# Select PW results from results
mum_pws <- raw_mum_rslts %>% 
  map(~.x$pathways) 

# Get pathway names
mum_pw_names <- mum_pws %>% 
  map(~.x$name) 

# Get list of cpds
mum_pw_ecs <- mum_pws %>% 
  map(~.x$emp_cpds) 

# Merge and unlist nested dataframes
ecd_pw_key <- map2(mum_pw_ecs, mum_pw_names, 
                   ~name_ecpathwayfxn(.x, .y)) %>% 
  map_depth(.depth = 2, 
            ~enframe (.)) %>% 
  map(~bind_rows(.x, .id = "pathway") %>% 
        rename(empcpd_num = name)) %>% 
  bind_rows(.id = "chrt_pfas_mode")

# get cohort, pfas, and mode
ecd_pw_key_2 <- ecd_pw_key %>% 
  mutate(cohort = str_split_fixed(chrt_pfas_mode,"_", n = 3)[,1], 
         effect = str_split_fixed(chrt_pfas_mode,"_", n = 3)[,2]) %>%
  select(cohort, effect, pathway, empcpd_num, value) %>% 
  rename(empirical_compound = value) 

# Remove duplicates
ecd_pw_key_final <- ecd_pw_key_2 %>% 
  select(pathway, empirical_compound) %>%
  tidylog::distinct()

# Combine mzrt key and ecd_pw_key -----------------------------------------
final_key <- tidylog::full_join(ecd_pw_key_final, 
                                mzrtkey_primary_ions, 
                                by = c("empirical_compound")) %>% 
  tidylog::filter(!is.na(feature)) %>% 
  rename(name = feature) %>% 
  distinct()

# Get Compound Names ----------------------------------------------------------
# This data has 1,212 unique matched compounds:
length(unique(final_key$matched_compound))

# Save data to CSV to create key linking compound ids to compound names
# write_csv(data.frame(compound_id = unique(final_key$matched_compound)), 
#           fs::path(dir_data_local,
#                    "Supporting Files", 
#                    "Cpd id to name keys",
#                    "Mummichog cpd id to cpd name key raw.csv"))
# 
# write_rds(final_key, fs::path(dir_data_local,
#                               "Supporting Files", 
#                               "mummichog_pw_ec_feature_key_cpd_id_only.rds"))