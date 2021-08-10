# 0.1.2 exposure_outcome Data Cleaning

# Load Exposure Outcome Data from Server 
load(file = fs::path(dir_solar_data_secure, 
                 "11.Jesse EDCs T2D 2021",
                 "Datasets", "Final data",
                 "All Final Datasets with HRE PFAS.Rdata", 
                 sep = ""))

# Remove unnecessary datasets
rm(sol_longitudinal, sol_ogtt, chs_ogtt)

# Rename Datasets 
solar_exposure_outcome <- sol_baseline; rm(sol_baseline)
chs_exposure_outcome <- chs; rm(chs)


# Calculate additional exposure variables: OC Chemichals --------------------------
solar_exposure_outcome <- solar_exposure_outcome  %>%
  mutate(across(contains("detect"), 
                ~if_else(str_detect(.,"non"), 0, 1))) %>% 
  mutate(pcb_num_detect = pcb_180_ngml_detect+pcb_153_ngml_detect+
           pcb_138_ngml_detect+pcb_118_ngml_detect, 
         pbde_num_detect = pbde_154_ngml_detect + 
           pbde_47_ngml_detect + 
           pbde_100_ngml_detect + 
           pbde_153_ngml_detect+
           pbde_85_ngml_detect, 
         ocs = dde_impute + hexachlorobenzene_impute) %>% 
  filter(id %in% solar_c18$id)

#CHS 
chs_exposure_outcome <- chs_exposure_outcome %>% 
  mutate(across(contains("detect"), 
                ~if_else(. == "detect", 1, 0))) %>% 
  mutate(pcb_num_detect = pcb_180_ngml_detect+pcb_153_ngml_detect+ pcb_138_ngml_detect+pcb_118_ngml_detect, 
         pbde_num_detect = pbde_154_ngml_detect + pbde_47_ngml_detect + 
           pbde_100_ngml_detect + pbde_153_ngml_detect+pbde_85_ngml_detect, 
         ocs = dde_impute + hexachlorobenzene_impute) 



##  Remove unnecessary variables


  
# Filter out non-physiological responses
chs_exposure_outcome <- chs_exposure_outcome %>% filter(og_glu_5 < 240)

# rm(cuberoot, pfas_dicot, transpose_ft, untranspose_ft, vzdt)

