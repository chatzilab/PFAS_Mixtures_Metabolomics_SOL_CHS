library(janitor)
library(jag2)
# library(tidylog)

# load formats for variables
source(here::here("0_0_1_format_vars_funs.R"))
source(here::here("!directories.R"))

# Load Metabolomics Feature Tables --------------------------------------
ftdata <- read_rds(fs::path(dir_data, 
                            "feature_tables_batch_corrected.rds"))

## Remove qc samples and get list of ids 
metab_samples <- ftdata$feature_tables$c18neg$id
# IDs length 5 are CHS IDs, length 13 are SOLAR ids
metab_ids <- metab_samples[nchar(metab_samples) == 5 | 
                             nchar(metab_samples) == 13]

# Get Feature Tables by cohort: SOLAR
fts_sol <-ftdata$feature_tables %>% 
  modify(~filter(., id %in% metab_ids, 
                 str_detect(id, "sol")) %>% 
           select(-class))

# Get Feature Tables by cohort: CHS
fts_chs <-ftdata$feature_tables %>% 
  modify(~filter(., id %in% metab_ids, 
                 str_detect(id, "sol", negate = TRUE)) %>% 
           select(-class))


met <- list(solar = fts_sol, 
            chs = fts_chs)

# Load Sample Metadata
samp_metadata <- read_rds(fs::path(dir_data,
                                   "sample_metadata_with_summaries.rds"))

# Load Exposure Outcome Data from drive  ------------------------
load(fs::path(dir_metabolomics, 
              "All Final Datasets with HRE PFAS.Rdata"))

# Remove unnecessary datasets
rm(sol_longitudinal, sol_ogtt, chs_ogtt)

# Rename Datasets 
solar_exposure_outcome <- sol_baseline; rm(sol_baseline)
chs_exposure_outcome <- chs; rm(chs)


# Calculate additional exposure variables: OC Chemichals --------------------------
solar_exposure_outcome <- solar_exposure_outcome  %>%
  mutate(across(contains("detect"), 
                ~if_else(str_detect(.,"non"), 0, 1)), 
         across(all_of(exposures_continuous), 
                ~log2(.),
                .names = "lg2_{col}")) %>% 
  mutate(pcb_num_detect = pcb_180_ngml_detect+pcb_153_ngml_detect+
           pcb_138_ngml_detect+pcb_118_ngml_detect, 
         pbde_num_detect = pbde_154_ngml_detect + 
           pbde_47_ngml_detect + 
           pbde_100_ngml_detect + 
           pbde_153_ngml_detect+
           pbde_85_ngml_detect, 
         ocs = dde_impute + hexachlorobenzene_impute) %>% 
  filter(id %in% metab_ids) 


#CHS 
chs_exposure_outcome <- chs_exposure_outcome %>% 
  mutate(across(contains("detect"), 
                ~if_else(str_detect(.,"non"), 0, 1)), 
         across(all_of(exposures_continuous), 
                ~log2(.),
                .names = "lg2_{col}")) %>% 
  mutate(pcb_num_detect = pcb_180_ngml_detect+pcb_153_ngml_detect+ pcb_138_ngml_detect+pcb_118_ngml_detect, 
         pbde_num_detect = pbde_154_ngml_detect + pbde_47_ngml_detect + 
           pbde_100_ngml_detect + pbde_153_ngml_detect+pbde_85_ngml_detect, 
         ocs = dde_impute + hexachlorobenzene_impute) 


#bind to list
exposure_outcome = list(solar = solar_exposure_outcome,
                        chs = chs_exposure_outcome)

# Clean Environment
rm(chs_exposure_outcome, solar_exposure_outcome, ftdata, fts_chs, fts_sol)




# Annotations 
common_metabolite_annotation <- read_rds(fs::path(
  dir_data, 
  "4_Common_Metabolites_Annotation", 
  "Common_Metabolites_SOLAR_CHS_V1.RDS"))
