# script for MWAS analysis (regression)
# C8 and hilic metabolites separated
# Jesse Goodrich 081721

library(broom)
source(here::here("0_0_3_analysis_funs.R"))

# Get vector of lcms modes for analysis
modes <- names(met$chs)

# chs mwas ------------------------------
# Initialize list of results
chs_mwas_results <- vector("list", length(modes))
names(chs_mwas_results) <- modes

## Run MWAS for All modes ---------------
for(l in modes){
  print(paste("Start", l, Sys.time()))
  # Run MWAS for all exposures
  mode_mwas <- mwas_all_exposures(cohort = "chs", 
                                  metab_dat = met,
                                  lcms_mode = modes[1],
                                  exp_cov_dat = exposure_outcome, 
                                  exposures = exposures) %>% 
    list()
  
  # Add to results list
  chs_mwas_results[l] = mode_mwas
}

#version 1: new LC-MS Data
saveRDS(chs_mwas_results,
        file = here::here('Temporary results',
                          exposure_type,
                          '1_1_CHS_mwas_results_two_hour_met.rds'))

