library(janitor)
library(jag2)
# library(tidylog)

# load formats for variables
source(here::here("0_0_1_format_vars_funs.R"))


# Load Metabolomics Feature Tables --------------------------------------
ftdata <- read_rds(fs::path(dir_data, 
                            "sol_chs_batch_cor_scaled_untargeted_fts.rds"))

# Load Sample Metadata
samp_metadata <- read_rds(fs::path(dir_data,
                                   "sample_metadata_with_summaries.rds"))

# Load Exposure Outcome Data from drive  ------------------------
load(fs::path(dir_data, 
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
  filter(id %in% ftdata$solar$c18neg$id) 


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
rm(chs_exposure_outcome, solar_exposure_outcome, ftdata)


  
# Load annotated data ----------------------------------------------
annotated_fts <- read_rds(
  fs::path(dir_data, 
           "sol_chs_batch_cor_scaled_annotated_fts.rds"))
