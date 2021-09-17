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

exposures = c("netfosaa","nmefosaab","pfbs","pfda","pfdoa","pfds","pfhpa",
              "pfhps","pfhxa","pfhxs","pfna","pfns","pfoa","pfos","pfpes",
              "pfuda", "x82fts")
cohort = c("solar", "chs")
modes = c("c18pos","c18neg", "hilicpos", "hilicneg")

# Temp Results folder architecture:
## Exposure Type > Exposure Name > Cohort > modes

# exposure_name = "netfosaa"
# chrt = cohort[1]
# mode = modes[1]
# rm(exposure_name, chrt, mode)

# Run all analyses -------------------------
for(exposure_name in exposures){
  for(chrt in cohort){
    for(mode in modes){
      print(paste("Begin",exposure_name,chrt, mode))
      t1 <- Sys.time()
      #  Set working directory for first folder
      setwd(fs::path(dir_temp,exposure_type,exposure_name, chrt, mode))

      # Check to see if mwas results are in the file:
      if(file.info(paste(chrt, 
                        exposure_name,
                        mode, 
                        "MWAS.csv", sep="_"))$size > 48){
        
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
                            exposure_name, 
                            mode, 
                            "MWAS.csv", 
                            sep = "_")) #change this
        
        
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
                             exposure_name, 
                             mode, 
                             "mummichog.RDS", 
                             sep = "_"))
        
        
        ## Get pathway data from mum results -------------
        pathways <- Set$mummi.resmat
        
        pathways <- as_tibble(pathways, rownames = "pathway") %>% 
          janitor::clean_names() %>% 
          mutate(mode = mode, 
                 exposure = exposure_name) %>% 
          select(exposure, mode, everything())
        
        # Write pathway results 
        write_csv(pathways, 
                  fs::path(dir_temp, 
                             exposure_type,
                             "mum_pathway_results",
                             chrt,
                             paste(exposure_name,
                                   mode,
                                   "mum_pathways.csv", 
                                   sep = "_")))
        
        print(Sys.time()-t1)
        print(paste(exposure_name,chrt, mode,  "Complete"))
        rm(mSet)
      } else{ 
        print(paste(exposure_name,chrt, mode, 
                    ": mummichog not performed (no non-na values)"))
        
      }
    }
  }
}
  
  
  
