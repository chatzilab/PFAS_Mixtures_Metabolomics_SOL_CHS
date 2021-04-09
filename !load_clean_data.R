library(tidyverse)
library(fs)
setwd(here::here())
source(here::here("!directories.R"))

# load formats for variables
source(here::here("0.0.1_format_vars_funs.R"))

# load separate data sets

##  SOLAR
solar_c18 <- read_csv(fs::path(dir_solar_data_secure_temp, 
                         "metabolomics_data_cleaned",
                         "solar_c18_final_ft.csv"))
solar_hilic <- read_csv(fs::path(dir_solar_data_secure_temp, 
                           "metabolomics_data_cleaned",
                           "solar_hilic_final_ft.csv")) 

##  CHS
chs_c18 <- read_csv(fs::path(dir_data_chs, 
                         "metabolomics_data_cleaned",
                         "metaair_c18_final_ft.csv"))
chs_hilic <- read_csv(fs::path(dir_data_chs, 
                           "metabolomics_data_cleaned",
                           "metaair_hilic_final_ft.csv"))



source(here::here("0.1.2_SOLAR_exposure_outcome_data_cleaning.R"))
#
 


# # combine data into single dataframe for analysis
# df_selection_messy <- 
#   
#   # liver outcome
#   liv_enz_data$liver_enz %>%
#   
#   # exposome covariates
#   left_join(
#     exposome_covs$table1_data_select,
#     by = c("helixid", "h_cohort", "alt_n", "ast_n", "ggt_n", "ck18_n", 
#            "i_q90_alt", "i_q90_ast", "i_q90_ggt", "i_q90_ck18", "i_q90_any", "i_alt_sex")
#   ) %>% 
#   
#   # SNP exposure
#   left_join(g_snp_genome$g_snps, by = "helixid") %>% 
#   
#   # transcriptome intermediates # almost 300 participants with transcriptomics measurements not in liver enzyme data...is it an error?
#   #  no, there were extra samples taken. From QC report:
#   #  "308 extra samples from three HELIX cohorts (extra HELIX samples)"
#   left_join(lst_transcriptome$helix_trnscrptm, by = "helixid") %>%
#   
#   # proteome intermediates
#   left_join(lst_proteome$proteome, by = "helixid") %>% 
#   
#   # serum metabolome intermediates
#   left_join(lst_metabolome_serum$metab_serum, by = "helixid") %>% 
#   
#   # urine metabolome intermediates
#   left_join(lst_metabolome_urine$metab_urine, by = "helixid") %>% 
#   
#   dplyr::arrange(helixid)
# 
# 
# # df_selection_messy %>% glimpse
