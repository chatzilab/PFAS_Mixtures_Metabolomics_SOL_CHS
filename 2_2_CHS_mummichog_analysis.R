# RUN MUMMICHOG METABOANALYST code
library(MetaboAnalystR)
library(tidyverse)
library(fs)
library(here)

##  Set Base Working Directory
# rm(list = ls())
exposure_type = "PFAS"

##  Get names of folders in Results directory
folders <- list.files(here::here("Temporary results", exposure_type))

folders <- folders[folders %in% exposures]

# Create loop to automaticially run all analysis
for(i in 1:length(folders)){
  # Clean MWAS data  --------------------------------------------------------
  print(paste("Begin", folders[i]))
  # Create chs Folder
  dir.create(file.path(here::here("Temporary results", 
                                  exposure_type, 
                                  folders[i], 
                                  "chs")), 
             showWarnings = TRUE)
  
  #  Set working directory for first folder, input data
  setwd(here::here("Temporary results",
                   exposure_type,
                   folders[i], 
                   "chs"))
  
  # # Create C18 folder
  # dir.create(file.path(here::here("Temporary results",
  #                                 exposure_type, 
  #                                 folders[i], 
  #                                 "chs", 
  #                                 "c18")), 
  #            showWarnings = FALSE)
  
  # Create HILIC folder
  dir.create(file.path(here::here("Temporary results",
                                  exposure_type, 
                                  folders[i], 
                                  "chs", 
                                  "hilic")),
             showWarnings = FALSE)
  
  
  ##  Read in mwas results
  hilic <- read_rds(fs::path(dir_home, 
                             "1_Code", 
                             "Temporary results",
                             exposure_type, 
                             "1_1_chs_hilic_mwas.rds")) %>% 
    filter(exposure == folders[i]) %>% 
    select(m.z, p.value, t.score, r.t)
  
  
  
  if(sum(is.na(hilic$p.value)) == 0){
    
    ##  Save hilic data
    write_csv(hilic, fs::path("hilic", paste(folders[i], "MWAS_hilic.csv", sep = "_")))
    
    ## Clean Working Directory
    rm(hilic)
    
    # Run mummichog on HILIC positive  ----------------------------------------
    setwd(here::here("Temporary results", exposure_type, folders[i], "chs", "hilic"))
    
    ##  Create objects for storing processed data from the MS peaks to pathways module
    mSet <- InitDataObjects(data.type = "mass_all",
                            anal.type = "mummichog", 
                            paired = FALSE)
    
    ##  Set peak format - contains m/z features, p-values, t-scores, and r.t.
    mSet <- SetPeakFormat(mSetObj = mSet, "mptr")
    ## Update metaboanalyst parameters: 
    mSet<-UpdateInstrumentParameters(mSet, 
                                     force_primary_ion = TRUE,
                                     instrumentOpt = 5.0, # Mass accuracy
                                     msModeOpt = "positive",  # mode
                                     rt_tol = 20) 
    # Read in MWAS results
    mSet <- Read.PeakListData(mSetObj = mSet, 
                              filename = paste(folders[i], 
                                               "MWAS_hilic.csv", 
                                               sep = "_")) #change this
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
    write_csv(pathways, here::here("Temporary results", exposure_type,
                                   "mum_pathway_results",
                                   "chs",
                                   paste(folders[i], "mum_pathway_hilic.csv", sep = "_")))
    
    print(paste(folders[i], "Complete"))
    
    rm(mSet)
  } else{ 
    print(paste(folders[i], ": mummichog not performed (no non-na values)"))
  }
  
  
  # # Run mummichog on C18 negative  ----------------------------------------
  # 
  ##  Read in mwas results
  # c18 <- read_rds(fs::path(dir_home, 
  #                      "1_Code", 
  #                      "Temporary results",
  #                      exposure_type, 
  #                      "1_1_chs_c18_mwas.rds")) %>% 
  #   filter(exposure == folders[i]) %>% 
  #   select(m.z, p.value, t.score, r.t)
  # 
  ##  Save c18 data
  # write_csv(c18, 
  #           fs::path("c18", paste(folders[i], "MWAS_c18.csv", sep = "_")))
  # 
  ## Clean Working Directory
  # rm(c18)
  # rm(Set, mSet)
  # # Set wd for c18 analysis
  # setwd(here::here("Temporary results", folders[i],"chs", "c18"))
  # 
  # ##  Create objects for storing processed data from the MS peaks to pathways module
  # mSet<-InitDataObjects(data.type = "mass_all",
  #                       anal.type = "mummichog", 
  #                       paired = FALSE)
  # ##  Set peak format - contains m/z features, p-values and t-scores
  # SetPeakFormat("mptr")
  # ## Update metaboanalyst parameters: 
  # mSet<-UpdateInstrumentParameters(mSet, 
  #                                  force_primary_ion = TRUE,
  #                                  instrumentOpt = 5.0, # Mass accuracy
  #                                  msModeOpt = "positive",  # mode
  #                                  rt_tol = 20) 
  # # Read in MWAS results
  # mSet <- Read.PeakListData(mSetObj = mSet, 
  #                           filename = paste(folders[i], "MWAS_c18.csv", sep = "_")) #change this
  # # Perform Sanity check
  # mSet<-SanityCheckMummichogData(mSet)
  # # Set peak enrichment algorithm options
  # mSet <- SetPeakEnrichMethod(mSet, algOpt = "mum", version = "v2")
  # ## Set mummichog p value:  # THIS MATTERS A LOT!
  # mSet <- SetMummichogPval(mSet, 0.05)
  # ##  Perform mummichog algorithm with 10,000 permutations
  # Set <- PerformPSEA(mSet, "hsa_mfn", "current", minLib = 2, permNum =  10000)
  # ## Save output
  # saveRDS(Set, file = paste0(folders[i], "_c18_mumichog_results.RDS"))
  # 
  # ## Get pathway data 
  # pathways <- Set$mummi.resmat
  # 
  # pathways <- as_tibble(pathways, rownames = "pathway") %>% 
  #   janitor::clean_names() %>% 
  #   mutate(mode = "c18", 
  #          exposure = folders[i]) %>% 
  #   select(exposure, mode, everything())
  # 
  # # Write pathway results 
  # write_csv(pathways, here::here("Temporary results", 
  #                                "mum_pathway_results",
  #                                "chs",
  #                                paste(folders[i], "mum_pathway_c18.csv", sep = "_")))
}
      