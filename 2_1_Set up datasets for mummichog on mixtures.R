# Set up datasets for Mummichog

source(here::here("!directories.R"))
##  Set Base Working Directory
# rm(list = ls())
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

dir_results_mixtures <- fs::path(dir_results, "PFAS_mixtures", "mummichog")

# set up folder structure ------------------------------
for(exposure_name in exposures){
  # Create Exposure Folder: Level 1 (Exposure)
  dir.create(file.path(fs::path(dir_results_mixtures, 
                                exposure_name)), 
             showWarnings = TRUE)
  # Level 2 (Cohort): Create SOLAR Folder
  dir.create(file.path(fs::path(dir_results_mixtures,exposure_name, 
                                "solar")), 
             showWarnings = TRUE)
  # Level 2 (Cohort): Create CHS Folder
  dir.create(file.path(fs::path(dir_results_mixtures,exposure_name, 
                                "chs")), 
             showWarnings = TRUE)
  # Level 3 (Cohort): Create Modes
  for(mode in modes){
    #solar
    dir.create(file.path(fs::path(dir_results_mixtures,exposure_name, 
                                  "solar", mode)),
               showWarnings = TRUE)
    #chs
    dir.create(file.path(fs::path(dir_results_mixtures,exposure_name, 
                                  "chs", mode)),
               showWarnings = TRUE)
  }
  rm(mode, exposure_name)
}


# save mwas in analysis folders ------------------------------
sol_mwas  <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "sol_pfas_mixtures_results_final_v3.csv")) %>% 
  as_tibble()

chs_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "chs_pfas_mixtures_results_final_v3.csv")) %>% 
  as_tibble()


# Join mwas for each cohort into one list
mwas_results <- list(solar = sol_mwas, 
                     chs = chs_mwas)
# Clean environment
rm(sol_mwas, chs_mwas)


# Loop over exposures, cohort, and mode to save mwas data
for(exposure_name in exposures){
  print(paste("Begin", exposure_name))
  t1 <- Sys.time()
  for(chrt in cohort){
    for(mde in modes){
      # Subset by mode and cohort
      temp_mwas <- mwas_results[[chrt]] %>% filter(mode == mde)
      
      temp_mwas2 <- temp_mwas %>% 
        dplyr::filter(str_detect(exposure, exposure_name)) %>% 
        dplyr::mutate(m.z = str_split_fixed(feature, pattern = "_", n = 2)[,1] %>% 
                        as.numeric(), 
                      r.t = str_split_fixed(feature, pattern = "_", n = 2)[,2] %>% 
                        as.numeric()) %>%
        dplyr::rename(p.value = p_value, 
                      t.score = wald) %>%
        dplyr::select(m.z, p.value, t.score, r.t)
      
      # Save mwas to correct folder
      if(sum(is.na(temp_mwas2$p.value)) == 0){
        write_csv(temp_mwas2, 
                  fs::path(dir_results_mixtures, exposure_name, chrt, mde,
                           paste(chrt, 
                                 exposure_name, 
                                 mde, 
                                 "MWAS.csv", 
                                 sep = "_")))
      }
    }
  }
  print(Sys.time()-t1)
}