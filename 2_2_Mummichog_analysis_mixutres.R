# RUN MUMMICHOG METABOANALYST code
library(tidyverse)
library(fs)
library(MetaboAnalystR)


source(here::here("!directories.R"))
source(file = fs::path(dir_home,
                       "1_Code", 
                       "troubleshooting files",
                       "internal mum variables.R"))

# Set important vars
exposure_type = "PFAS"
exposures = c("Mixture effect",
              "pfda",
              "pfhps",
              "pfhxs",
              "pfna",
              "pfoa",
              "pfos")

cohort = c("solar", "chs")
modes = c("c18pos","c18neg", "hilicpos", "hilicneg")


# Get Dataframe of all conbinations
all_mum_comb <- expand.grid(list(modes = modes, 
                                 cohort = cohort, 
                                 exposures = exposures)) %>% 
  as_tibble() %>% 
  mutate(across(.cols = everything(), ~as.character(.)), 
         row = row_number(exposures)) %>% 
  dplyr::select(row, exposures, cohort, modes) 

rm(modes, cohort, exposures)

# # Error: "Begin pfhps solar c18neg" (row 18)

# Temp Results folder architecture:
## Exposure Name > Cohort > modes

i = 18
all_mum_comb$exposures[i]
# Run all analyses -------------------------
for(i in 1:nrow(all_mum_comb)){
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
  setwd(fs::path(dir_results_mum_mixtures,exposures, chrt, mode))
  
  xxx<-read_csv("solar_pfhps_c18neg_MWAS.csv")
  table(xxx$p.value<0.05)
  
  # Get mode for mummichog analysis
  pos_neg =if_else(str_detect(mode, "pos"), "positive", "negative")
  
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
                                   msModeOpt = pos_neg,  # mode
                                   rt_tol = 20) 
  # Read in MWAS results
  mSet <- Read.PeakListData(
    mSetObj = mSet, 
    filename =  paste(chrt, 
                      exposures,
                      mode, 
                      "MWAS.csv", 
                      sep = "_"))
  
  
  # Perform Sanity check
  mSet<-SanityCheckMummichogData(mSet)
  # Set peak enrichment algorithm options
  mSet <- SetPeakEnrichMethod(mSet, algOpt = "integ", version = "v2")
  ## Set mummichog p value:  # THIS MATTERS A LOT!
  mSet <- SetMummichogPval(mSet, 0.05)
  ##  Perform mummichog algorithm 
  Set <- PerformPSEA(mSet, 
                     "hsa_mfn", 
                     "current", 
                     minLib = 3, 
                     permNum =  10000)
  ## Save mummichog output
  saveRDS(Set, 
          file = paste(chrt, 
                       exposures, 
                       mode, 
                       "mummichog.RDS", 
                       sep = "_"))
  
  
  ## Get pathway data from mum results -------------
  pathways <- Set$mummi.resmat
  
  pathways <- as_tibble(pathways, rownames = "pathway") %>% 
    janitor::clean_names() %>% 
    mutate(mode = mode, 
           exposure = exposures) %>% 
    select(exposure, mode, everything())
  
  # Write pathway results 
  write_csv(pathways, 
            fs::path(dir_results_mum_mixtures, 
                     "mum_pathway_results",
                     chrt,
                     paste(exposures,
                           mode,
                           "mum_pathways.csv", 
                           sep = "_")))
  
  print(Sys.time()-t1)
  print(paste(all_mum_comb$exposures[i], chrt, mode, "Complete"))
  rm(mSet)
} 


