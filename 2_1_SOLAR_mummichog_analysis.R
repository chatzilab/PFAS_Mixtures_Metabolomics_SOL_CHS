# MUMMICHOG METABOANALYST code
library(MetaboAnalystR)
library(tidyverse)
library(fs)
library(here)

##  Set Base Working Directory
# rm(list = ls())

##  Get names of folders in Results directory
folders <- list.files(here("Temporary results"))

folders <- folders[folders %in% exposures]


for(i in 1:length(folders)){
  # Clean MWAS data  --------------------------------------------------------
  #  Set working directory for first folder, input data
  setwd(here::here("Temporary results",
                   folders[i],  "SOLAR"))
  
  dir.create(file.path(here::here("Temporary results",
                                  folders[i],  "SOLAR"), 
                       "c18"), showWarnings = FALSE)
  
  dir.create(file.path(here::here("Temporary results",
                                  folders[i],  "SOLAR"), 
                       "hilic"), showWarnings = FALSE)
  
  
  ##  Read .csv
  mwas_results <- read_csv(list.files()[str_detect(list.files(), ".csv")])
  
  ##  Subset hilic & c18
  hilic <- mwas_results %>% filter(mode == "positive") %>% select(m.z, p.value, t.score, r.t)
  c18 <- mwas_results %>% filter(mode == "negative") %>% select(m.z, p.value, t.score, r.t)
  
  ##  Save hilic and c18
  write_csv(hilic, fs::path("hilic", paste(folders[i], "MWAS_hilic.csv", sep = "_")))
  write_csv(c18, fs::path("c18", paste(folders[i], "MWAS_c18.csv", sep = "_")))
  
  ## Clean Working Directory
  rm(mwas_results, hilic, c18)
  
  # Run mummichog on HILIC positive  ----------------------------------------
  setwd(here::here("Temporary results", folders[i], "SOLAR", "hilic"))
  
  ##  Create objects for storing processed data from the MS peaks to pathways module
  mSet<-InitDataObjects(data.type = "mass_all",
                        anal.type = "mummichog", 
                        paired = FALSE)
  ##  Set peak format - contains m/z features, p-values and t-scores
  SetPeakFormat("mptr")
  ## Update metaboanalyst parameters: 
  mSet<-UpdateInstrumentParameters(mSet, 
                                   force_primary_ion = TRUE,
                                   instrumentOpt = 5.0, # Mass accuracy
                                   msModeOpt = "positive",  # mode
                                   rt_tol = 20) 
  # Read in MWAS results
  mSet <- Read.PeakListData(mSetObj = mSet, 
                            filename = paste(folders[i], "MWAS_hilic.csv", sep = "_")) #change this
  # Perform Sanity check
  mSet<-SanityCheckMummichogData(mSet)
  # Set peak enrichment algorithm options
  mSet <- SetPeakEnrichMethod(mSet, algOpt = "mum", version = "v2")
  ## Set mummichog p value:  # THIS MATTERS A LOT!
  mSet <- SetMummichogPval(mSet, 0.05)
  ##  Perform mummichog algorithm 
  Set <- PerformPSEA(mSet, "hsa_mfn", "current", minLib = 2, permNum =  10000)
  ## Save output
  saveRDS(Set, file = paste0(folders[i], "_hilic_mumichog_results.RDS"))
  
  ## Get pathway data 
  pathways <- Set$mummi.resmat
  
  pathways <- as_tibble(pathways, rownames = "pathway") %>% 
    janitor::clean_names() %>% 
    mutate(mode = "hilic", 
           exposure = folders[i]) %>% 
    select(exposure, mode, everything())
  
  # Write pathway results 
  write_csv(pathways, here::here("Temporary results", 
                                 "mum_pathway_results",
                                 "SOLAR",
                                 paste(folders[i], "mum_pathway_hilic.csv", sep = "_")))
  
  # Run mummichog on C18 negative  ----------------------------------------
  rm(Set, mSet)
  # Set wd for c18 analysis
  setwd(here::here("Temporary results", folders[i],"SOLAR", "c18"))
  
  ##  Create objects for storing processed data from the MS peaks to pathways module
  mSet<-InitDataObjects(data.type = "mass_all",
                        anal.type = "mummichog", 
                        paired = FALSE)
  ##  Set peak format - contains m/z features, p-values and t-scores
  SetPeakFormat("mptr")
  ## Update metaboanalyst parameters: 
  mSet<-UpdateInstrumentParameters(mSet, 
                                   force_primary_ion = TRUE,
                                   instrumentOpt = 5.0, # Mass accuracy
                                   msModeOpt = "positive",  # mode
                                   rt_tol = 20) 
  # Read in MWAS results
  mSet <- Read.PeakListData(mSetObj = mSet, 
                            filename = paste(folders[i], "MWAS_c18.csv", sep = "_")) #change this
  # Perform Sanity check
  mSet<-SanityCheckMummichogData(mSet)
  # Set peak enrichment algorithm options
  mSet <- SetPeakEnrichMethod(mSet, algOpt = "mum", version = "v2")
  ## Set mummichog p value:  # THIS MATTERS A LOT!
  mSet <- SetMummichogPval(mSet, 0.05)
  ##  Perform mummichog algorithm with 10,000 permutations
  Set <- PerformPSEA(mSet, "hsa_mfn", "current", minLib = 2, permNum =  10000)
  ## Save output
  saveRDS(Set, file = paste0(folders[i], "_c18_mumichog_results.RDS"))
  
  ## Get pathway data 
  pathways <- Set$mummi.resmat
  
  pathways <- as_tibble(pathways, rownames = "pathway") %>% 
    janitor::clean_names() %>% 
    mutate(mode = "c18", 
           exposure = folders[i]) %>% 
    select(exposure, mode, everything())
  
  # Write pathway results 
  write_csv(pathways, here::here("Temporary results", 
                                 "mum_pathway_results",
                                 "SOLAR",
                                 paste(folders[i], "mum_pathway_c18.csv", sep = "_")))
}
