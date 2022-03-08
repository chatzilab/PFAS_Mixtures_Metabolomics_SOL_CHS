# Set up datasets for Mummichog

source(here::here("!directories.R"))
##  Set Base Working Directory
# rm(list = ls())
exposure_type = "PFAS"


exposures = c("netfosaa","nmefosaab","pfbs","pfda","pfdoa","pfds","pfhpa",
              "pfhps","pfhxa","pfhxs","pfna","pfns","pfoa","pfos","pfpes",
              "pfuda", "x82fts")
cohort = c("solar", "chs")
modes = c("c18pos","c18neg", "hilicpos", "hilicneg")

# set up folder structure ------------------------------
for(exposure_name in exposures){
  # Create Exposure Folder: Level 1 (Exposure)
  dir.create(file.path(fs::path(dir_results, exposure_type, 
                                exposure_name)), 
             showWarnings = TRUE)
  # Level 2 (Cohort): Create SOLAR Folder
  dir.create(file.path(fs::path(dir_results,exposure_type, exposure_name, 
                                "solar")), 
             showWarnings = TRUE)
  # Level 2 (Cohort): Create CHS Folder
  dir.create(file.path(fs::path(dir_results,exposure_type, exposure_name, 
                                "chs")), 
             showWarnings = TRUE)
  # Level 3 (Cohort): Create Modes
  for(mode in modes){
    #solar
    dir.create(file.path(fs::path(dir_results,exposure_type, exposure_name, 
                                  "solar", mode)),
               showWarnings = TRUE)
    #chs
    dir.create(file.path(fs::path(dir_results,exposure_type, exposure_name, 
                                  "chs", mode)),
               showWarnings = TRUE)
  }
  rm(mode, exposure_name)
}


# save mwas in analysis folders ------------------------------
mwas_results_solar <- read_rds(fs::path(dir_results, exposure_type, 
                                        "1_1_SOLAR_mwas_results_two_hour_met.rds"))

mwas_results_chs <- read_rds(fs::path(dir_results, exposure_type, 
                                      "1_1_CHS_mwas_results_two_hour_met.rds"))

# Join mwas for each cohort into one list
mwas_results <- list(solar = mwas_results_solar, 
                     chs = mwas_results_chs)
# Clean environment
rm(mwas_results_solar, mwas_results_chs)

# Loop over exposures, cohort, and mode to save mwas data
for(exposure_name in exposures){
  print(paste("Begin", exposure_name))
  t1 <- Sys.time()
  for(chrt in cohort){
    for(mode in modes){
      # Subset by mode and cohort
      temp_mwas <- mwas_results[[chrt]][[mode]]
      
      temp_mwas2 <- temp_mwas %>% 
        dplyr::filter(str_detect(exposure, exposure_name)) %>% 
        dplyr::mutate(m.z = str_split_fixed(name, pattern = "_", n = 2)[,1] %>% 
                        as.numeric(), 
                      r.t = str_split_fixed(name, pattern = "_", n = 2)[,2] %>% 
                        as.numeric()) %>%
        dplyr::rename(p.value = p_value, 
                      t.score = statistic) %>%
        dplyr::select(m.z, p.value, t.score, r.t)
      
      # Save mwas to correct folder
      if(sum(is.na(temp_mwas2$p.value)) == 0){
        write_csv(temp_mwas2, 
                  fs::path(dir_results, exposure_type, exposure_name, chrt, mode,
                           paste(chrt, 
                                 exposure_name, 
                                 mode, 
                                 "MWAS.csv", 
                                 sep = "_")))
      }
    }
  }
  print(Sys.time()-t1)
}


