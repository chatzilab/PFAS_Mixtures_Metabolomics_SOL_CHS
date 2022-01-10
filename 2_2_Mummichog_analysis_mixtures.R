# RUN MUMMICHOG METABOANALYST code
library(tidyverse)
library(MetaboAnalystR)
# library(igraph)

source(here::here("!directories.R"))
source(file = fs::path(dir_home,
                       "1_Code",
                       "troubleshooting files",
                       "internal mum variables.R"))

# Set important vars
exposure_type = "PFAS"
exposures = c("Mixture effect hyper_g")
cohort = c("solar", "chs")
modes = c("mixed")


# Get Dataframe of all conbinations
# all_mum_comb <- expand.grid(list(modes = modes, 
#                                  cohort = cohort, 
#                                  exposures = exposures)) %>% 
#   as_tibble() %>% 
#   mutate(across(.cols = everything(), ~as.character(.)), 
#          row = row_number(exposures)) %>% 
#   dplyr::select(row, exposures, cohort, modes) 


all_mum_comb <- data.frame(row = 1:2, 
           exposures = rep("Mixture effect hyper_g",2),
           cohort = c("solar", "chs"),
           modes = rep("mixed", 2))

 
rm(modes, cohort, exposures)

# Temp Results folder architecture:
## Exposure Name > Cohort > modes
i = 1
# Run all analyses -------------------------
# for(i in 1:nrow(all_mum_comb)){
  # Get names of variables
  exposures = all_mum_comb$exposures[i]
  chrt = all_mum_comb$cohort[i]
  mode = all_mum_comb$modes[i]
  
  print(paste("Begin",
              exposures,
              chrt,
              mode, 
              "-- Analysis", i, "out of", nrow(all_mum_comb)))
  t1 <- Sys.time()
  #  Set working directory for first folder
  setwd(fs::path(dir_results_mum_mixtures,
                 exposures, chrt))
  
  
  ## Set up for mummichog -------------------------
  ##  Create objects for storing data from the MS peaks to pathways module
  mSet <- InitDataObjects(data.type = "mass_all",
                          anal.type = "mummichog", 
                          paired = FALSE)
  
  ##  Set peak format - contains m/z features, p-values, t-scores, and r.t.
  mSet <- SetPeakFormat(mSetObj = mSet, "mptr")
  mSet<-SetRTincluded(mSet, "seconds")
  
  ## Update metaboanalyst parameters: 
  mSet<-UpdateInstrumentParameters(mSet,
                                   force_primary_ion = TRUE,
                                   instrumentOpt = 5.0, # Mass accuracy
                                   msModeOpt = "mixed",  # mode
                                   rt_tol = 20) 
  # Read in MWAS results
  mSet <- Read.PeakListData(
    mSetObj = mSet, 
    filename =  paste(chrt, 
                      "mixture_effect_mixed_mode_MWAS.csv",
                      sep = "_"))
  
  
  # Perform Sanity check
  mSet<-SanityCheckMummichogData(mSet)
  # Set peak enrichment algorithm options
  mSet <- SetPeakEnrichMethod(mSet, algOpt = "integ", version = "v2")
  ## Set mummichog p value:  # THIS MATTERS A LOT!
  mSet <- SetMummichogPval(mSet, 0.2)
  ##  Perform mummichog algorithm 
  Set <- PerformPSEA(mSet, 
                     "hsa_mfn", 
                     "current", 
                     minLib = 3, 
                     permNum =  10000)
  ## Save mummichog output
  saveRDS(Set, 
          file = paste(chrt, 
                       "mixture_effect_mixed_mode_mummichog.RDS", 
                       sep = "_"))
  
  
  ## Get pathway data from mum results -------------
  pathways <- Set$mummi.resmat
  
  pathways <- as_tibble(pathways, 
                        rownames = "pathway") %>% 
    janitor::clean_names() %>% 
    mutate(mode = mode, 
           exposure = exposures) %>% 
    select(exposure, mode, everything())
  
  # Write pathway results 
  # write_csv(pathways, 
  #           fs::path(dir_results_mum_mixtures, 
  #                    "mum_pathway_results_hyper_g",
  #                    paste(chrt, 
  #                          exposures,
  #                          mode,
  #                          "mum_pathways.csv", 
  #                          sep = "_")))
  
  print(Sys.time()-t1)
  print(paste(all_mum_comb$exposures[i], chrt, mode, "Complete"))
  rm(mSet)
} 
