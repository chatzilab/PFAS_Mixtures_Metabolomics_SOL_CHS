# script for MWAS analysis (regression)
# Jesse Goodrich 081721

library(broom)
source(here::here("0_0_3_analysis_funs.R"))

# Get vector of lcms modes for analysis
modes <- names(ftdata$solar)

# Solar mwas ------------------------------
# Initialize list of results
solar_mwas_results <- vector("list", length(modes))
names(solar_mwas_results) <- modes

## Run MWAS for All modes ---------------
# Get smaller working dataset
# temp_met <- list(solar = met$solar %>% modify(~select(., 1:10)),
#                  chs = met$chs %>% modify(~select(., 1:10)) )

# For Loop: For all modes, run the analysis
for(l in modes){
  print(paste("Start", l, Sys.time()))
  # Run MWAS for all exposures
  mode_mwas <- mwas_all_exposures(
    cohort = "solar", 
    metab_dat = met,
    lcms_mode = l,
    exp_cov_dat = exposure_outcome, 
    analysis_exposures = exposures_for_analysis) %>% 
    list()
  
  # Add to results list
  solar_mwas_results[l] = mode_mwas
}

# Save results
saveRDS(solar_mwas_results,
        file = fs::path(dir_results,
                        exposure_type,
                        '1_1_SOLAR_mwas_results_two_hour_met.rds'))


# CHS MWAS -------------------------------
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
                                  lcms_mode = l,
                                  exp_cov_dat = exposure_outcome, 
                                  analysis_exposures = exposures_for_analysis) %>% 
    list()
  
  # Add to results list
  chs_mwas_results[l] = mode_mwas
}

#version 1: new LC-MS Data
saveRDS(chs_mwas_results,
        file = fs::path(dir_results,
                        exposure_type,
                        '1_1_CHS_mwas_results_two_hour_met.rds'))

