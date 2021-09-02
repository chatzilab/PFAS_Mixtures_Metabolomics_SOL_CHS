# RUN MUMMICHOG METABOANALYST code
library(tidyverse)
library(fs)
library(MetaboAnalystR)


source(here::here("!directories.R"))
##  Set Base Working Directory
# rm(list = ls())
exposure_type = "PFAS"


exposures = c("netfosaa","nmefosaab","pfbs","pfda","pfdoa","pfds","pfhpa",
              "pfhps","pfhxa","pfhxs","pfna","pfns","pfoa","pfos","pfpes",
              "pfuda", "x82fts")


##  Get names of folders in Results directory
folders <- list.files(here::here("Temporary results", exposure_type))

folders <- exposures

# Create loop to automaticially run all analysis
for(i in 2:length(folders)){
  # Clean MWAS data  --------------------------------------------------------
  print(paste("Begin", folders[i]))
  
  # Create Exposure Folder
  dir.create(file.path(here::here("Temporary results", 
                                  exposure_type, 
                                  folders[i])), 
             showWarnings = TRUE)
  
  # Create SOLAR Folder
  dir.create(file.path(here::here("Temporary results", 
                                  exposure_type, 
                                  folders[i], 
                                  "SOLAR")), 
             showWarnings = TRUE)
  
  #  Set working directory for first folder, input data
  setwd(here::here("Temporary results",exposure_type,folders[i], 
                   "SOLAR"))
  
  # Create C18 folder
  dir.create(file.path(here::here("Temporary results",exposure_type,folders[i], 
                                  "SOLAR", "c18pos")), showWarnings = FALSE)
  
  dir.create(file.path(here::here("Temporary results",exposure_type,folders[i], 
                                  "SOLAR", "c18neg")), showWarnings = FALSE)
  
  
  # Create HILIC folder
  dir.create(file.path(here::here("Temporary results",exposure_type,folders[i], 
                                  "SOLAR", "hilicpos")), showWarnings = FALSE)
  
  dir.create(file.path(here::here("Temporary results",exposure_type,folders[i], 
                                  "SOLAR", "hilicneg")), showWarnings = FALSE)
  
  
  
  
  ##  Read in mwas results (list of 4), select "HILIC Pos" from list
  hilicpos <- read_rds(
    fs::path(dir_home, 
             "1_Code", 
             "Temporary results",
             exposure_type, 
             "1_1_SOLAR_mwas_results_two_hour_met.rds"))[["hilicpos"]] %>% 
    dplyr::filter(exposure == folders[i]) %>% 
    dplyr::mutate(m.z = str_split_fixed(name, pattern = "_", n = 2)[,1] %>% 
                    as.numeric(), 
           r.t = str_split_fixed(name, pattern = "_", n = 2)[,2] %>% 
             as.numeric()) %>%
    dplyr::rename(p.value = p_value, 
           t.score = statistic) %>%
    dplyr::select(m.z, p.value, t.score, r.t)
  
  
  
  if(sum(is.na(hilicpos$p.value)) == 0){
    
    ##  Save hilic data
    write_csv(hilicpos, fs::path("hilicpos", 
                                 paste(folders[i], "MWAS_hilicpos.csv", sep = "_")))
    
    ## Clean Working Directory
    rm(hilicpos)
    
    # Run mummichog on HILIC positive  ----------------------------------------
    setwd(here::here("Temporary results", exposure_type, folders[i], 
                     "SOLAR", "hilicpos"))
    
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
                                               "MWAS_hilicpos.csv", 
                                               sep = "_")) #change this
    # Perform Sanity check
    mSet<-SanityCheckMummichogData(mSet)
    # Set peak enrichment algorithm options
    mSet <- SetPeakEnrichMethod(mSet, algOpt = "integ", version = "v2")
    ## Set mummichog p value:  # THIS MATTERS A LOT!
    mSet <- SetMummichogPval(mSet, 0.05)
    ##  Perform mummichog algorithm 
    Set <- PerformPSEA(mSet, "hsa_mfn", "current", minLib = 2, permNum =  10000)
    ## Save output
    saveRDS(Set, file = paste0(folders[i], "_hilicpos_mumichog_results.RDS"))
    
    ## Get pathway data 
    pathways <- Set$mummi.resmat
    
    pathways <- as_tibble(pathways, rownames = "pathway") %>% 
      janitor::clean_names() %>% 
      mutate(mode = "hilicpos", 
             exposure = folders[i]) %>% 
      select(exposure, mode, everything())
    
    # Write pathway results 
    write_csv(pathways, here::here("Temporary results", exposure_type,
                                   "mum_pathway_results",
                                   "SOLAR",
                                   paste(folders[i], "mum_pathway_hilicpos.csv", sep = "_")))
    
    print(paste(folders[i], "Complete"))
    
    rm(mSet)
  } else{ 
    print(paste(folders[i], ": mummichog not performed (no non-na values)"))
  }
}
