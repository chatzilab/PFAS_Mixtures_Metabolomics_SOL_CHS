# Set up datasets for Mummichog. 
# This was originally written to run for 

source(here::here("!directories.R"))
##  Set Base Working Directory
# rm(list = ls())
exposure_type = "PFAS"
exposures = c("Mixture effect hyper_g")
cohort = c("solar", "chs")
modes = c("c18pos","c18neg", "hilicpos", "hilicneg")


dir_results_mum_mixtures 

# set up folder structure ------------------------------
# Create Exposure Folder: Level 1 (Exposure)
dir.create(file.path(fs::path(dir_results_mum_mixtures, 
                              exposures)), 
           showWarnings = TRUE)
# Level 2 (Cohort): Create SOLAR Folder
dir.create(file.path(fs::path(dir_results_mum_mixtures,
                              exposures, 
                              "solar")), 
           showWarnings = TRUE)
# Level 2 (Cohort): Create CHS Folder
dir.create(file.path(fs::path(dir_results_mum_mixtures,
                              exposures, 
                              "chs")), 
           showWarnings = TRUE)

# save mwas in analysis folders ------------------------------
mwas_results_cohort_list  <- read_rds(
  file = fs::path(dir_results_mixtures, 
                  "SOL CHS all Mixtures MWAS results long hyper_g_v2.rds"))


# Format data for mummichog (data structure: m (mz), p, t, r (retention time), mode (pos/neg))
temp_mwas <- mwas_results_cohort_list %>% 
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
write_csv(temp_mwas[["SOL_all pfas"]], 
          fs::path(dir_results_mum_mixtures, 
                   exposures, 
                   "solar",
                   "solar_all_pfas_mixture_effect_mixed_mode_MWAS.csv"))

write_csv(temp_mwas[["SOL_pfcas"]], 
          fs::path(dir_results_mum_mixtures, 
                   exposures, 
                   "solar",
                   "solar_pfcas_mixture_effect_mixed_mode_MWAS.csv"))

write_csv(temp_mwas[["SOL_pfsas"]], 
          fs::path(dir_results_mum_mixtures, 
                   exposures, 
                   "solar",
                   "solar_pfsas_mixture_effect_mixed_mode_MWAS.csv"))


# Save CHS data 
write_csv(temp_mwas[["CHS_all pfas"]],
          fs::path(dir_results_mum_mixtures,
                   exposures,
                   "chs",
                   "chs_all_pfas_mixture_effect_mixed_mode_MWAS.csv"))

write_csv(temp_mwas[["CHS_pfcas"]],
          fs::path(dir_results_mum_mixtures,
                   exposures,
                   "chs",
                   "chs_pfcas_mixture_effect_mixed_mode_MWAS.csv"))

write_csv(temp_mwas[["CHS_pfsas"]],
          fs::path(dir_results_mum_mixtures,
                   exposures,
                   "chs",
                   "chs_pfsas_mixture_effect_mixed_mode_MWAS.csv"))