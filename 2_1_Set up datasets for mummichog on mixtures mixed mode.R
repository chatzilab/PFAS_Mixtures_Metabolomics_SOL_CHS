# Set up datasets for Mummichog. 
# This was originally written to run for 

source(here::here("!directories.R"))
##  Set Base Working Directory
# rm(list = ls())
exposure_type = "PFAS"
exposures = c("Mixture effect w 09")
cohort = c("solar", "chs")
modes = c("c18pos","c18neg", "hilicpos", "hilicneg")


dir_results_mixtures <- fs::path(dir_results, "PFAS_mixtures", "mummichog")

# set up folder structure ------------------------------
# Create Exposure Folder: Level 1 (Exposure)
dir.create(file.path(fs::path(dir_results_mixtures, 
                              exposures)), 
           showWarnings = TRUE)
# Level 2 (Cohort): Create SOLAR Folder
dir.create(file.path(fs::path(dir_results_mixtures,exposures, 
                              "solar")), 
           showWarnings = TRUE)
# Level 2 (Cohort): Create CHS Folder
dir.create(file.path(fs::path(dir_results_mixtures,exposures, 
                              "chs")), 
           showWarnings = TRUE)

# save mwas in analysis folders ------------------------------
sol_mwas  <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "sol_pfas_mixtures_results_w_09.csv")) %>% 
  as_tibble()

chs_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "chs_pfas_mixtures_results_w_09.csv")) %>% 
  as_tibble()

xxx <- sol_mwas %>% filter(exposure == "mixture")
sum(xxx$p_value<0.2)#/length(xxx$p_value)


# Join mwas for each cohort into one list
mwas_results <- list(solar = sol_mwas, 
                     chs = chs_mwas)
# Clean environment
rm(sol_mwas, chs_mwas)


# Format data for mummichog (data structure: m (mz), p, t, r (retention time), mode (pos/neg))
temp_mwas <- mwas_results %>% 
  modify(. %>% 
           dplyr::mutate(mode = if_else(str_detect(mode, "neg"), 
                                        "negative",
                                        "positive"), 
                         m.z = str_split_fixed(feature, 
                                               pattern = "_", 
                                               n = 2)[,1] %>% 
                           as.numeric(), 
                         r.t = str_split_fixed(feature, 
                                               pattern = "_",
                                               n = 2)[,2] %>% 
                           as.numeric()) %>%
           dplyr::rename(p.value = p_value, 
                         t.score = wald) %>%
           dplyr::filter(exposure == "mixture") %>% 
           dplyr::select(m.z, p.value, t.score, r.t, mode))


# Save SOLAR data 
write_csv(temp_mwas[["solar"]], 
          fs::path(dir_results_mixtures, exposure_name, "solar",
                   "solar_mixture_effect_mixed_mode_MWAS.csv"))

# Save CHS data 
write_csv(temp_mwas[["chs"]], 
          fs::path(dir_results_mixtures, exposure_name, "chs",
                   "chs_mixture_effect_mixed_mode_MWAS.csv"))