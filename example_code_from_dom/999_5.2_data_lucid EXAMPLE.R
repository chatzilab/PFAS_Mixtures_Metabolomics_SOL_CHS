library(tidyverse)
library(Biobase)
library(BiocGenerics)


# join data: outcome, exposure, omics intermediates
df_lucid_messy <- liv_enz_data$liver_enz %>%          # liver outcome
  left_join(g_snp_genome$g_snps, by = "helixid") %>%  # SNP exposure
  left_join(lst_proteome$proteome, by = "helixid")    # proteome intermediate (okay to have missing Z)
# df_lucid_messy %>% glimpse

# covariates for outcome: sociodemographics
covs_sociodem <- exposome_covs$table1_data_select %>%
  dplyr::select(helixid, 
                exposome_covs$names_covs$cov_child, #  child sociodems
                exposome_covs$names_covs$cov_parent # parent sociodems
  ) %>% rename_with( ~ paste("cov_y", .x, sep = "_")) %>% dplyr::rename(helixid = cov_y_helixid) %>% 
  # convert factors to numeric
  dplyr::mutate(
    cov_y_e3_sex_None = as.numeric(cov_y_e3_sex_None) - 1,
    cov_y_h_ethnicity_c_None = as.numeric(cov_y_h_ethnicity_c_None) - 1,
    cov_y_e3_bwc_None = as.numeric(cov_y_e3_bwc_None) - 1,
    cov_y_h_edufc_None = as.numeric(cov_y_h_edufc_None), 
    cov_y_h_edumc_None = as.numeric(cov_y_h_edumc_None),
    cov_y_h_mbmic_None = as.numeric(cov_y_h_mbmic_None)
  )

# covariates for technical vars for latent cluster on genetic exposure
covs_techvar <- lst_proteome$proteome %>% 
  dplyr::select(helixid, all_of(lst_proteome$names_proteome$ptm_cov)) %>% 
  rename_with( ~ paste("cov_y", .x, sep = "_")) %>% dplyr::rename(helixid = cov_y_helixid) 

# add covariates to LUCID data
df_lucid_messy <- df_lucid_messy %>% 
  left_join(covs_sociodem, by = "helixid") %>% 
  left_join(covs_techvar, by = "helixid")

# clean data for LUCID
df_lucid <- df_lucid_messy %>%
  # remove participants with missing genetic data
  dplyr::filter_at(vars(starts_with("rs")), .vars_predicate = (function(x) !is.na(x))) %>% 
  # remove participants with missing proteomics data
  dplyr::filter_at(vars(starts_with("pro_")), .vars_predicate = (function(x) !is.na(x))) %>% 
  dplyr::filter_at(vars(starts_with("cov_y")), .vars_predicate = (function(x) !is.na(x)))

# how many participants have complete outcome, exposure and intermediate data?
df_lucid %>% nrow
# 674 observations with complete cases
  # only 985 of HELIX participants with Liver Outcomes have complete SNPs data
  # only 723 of HELIX participants with Liver Outcomes have complete proteomics data
                