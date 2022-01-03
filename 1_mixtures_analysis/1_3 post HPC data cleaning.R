# Mixtures Analysis 
library(purrr)
library(tidyverse)
# Read in data that will be loaded on the HPC 

# Read and restructure results from mixtures analysis performed on HPC

# The issue: Each node writes results to a single row, so data is not rectangular 
# (because each node does not run exactly the same number of models, as 23173 is
# not divisible by 128 nodes). 
source(here::here("!directories.r"))
source(here::here("!set_exposure_outcome_vars.r"))
source(here::here("!load_data.r"))

# SOLAR -------------------------------------------------------------------
# Read in data, without headers.
sol_og_for_colnames <- read.table(
  fs::path(dir_results, 
           "PFAS_mixtures",
           "from_hpc",
           "sol_pfas_mixtures_mwas_w_09.csv"),
  sep = ",", na.strings = "",
  as.is = TRUE, fill = TRUE, header = FALSE)

# Get temp col names. Col names are not indexed correctly though- they 
# Need to be shifted to the right by 1, and we need to 
row1 <- sol_og_for_colnames[1,] %>% as.character(.)
col_names <- c("beta",
               row1[-length(row1)], 
               "var.names.9999", 
               "Mean.9999", 
               "SD.9999",         
               "X2.5..9999",      
               "X97.5..9999",    
               "p.val.9999",      
               "metabolite.9999")

# Read in final data with headers
sol_og <- read.table(fs::path(dir_results, 
                              "PFAS_mixtures",
                              "from_hpc",
                              "sol_pfas_mixtures_mwas_w_09.csv"),
                     sep = ",",col.names = col_names,
                     na.strings = "",
                     as.is = TRUE, fill = TRUE, header = TRUE) %>% 
  dplyr::select(-beta)

# Clean Col names
sol_1 <- sol_og %>% 
  janitor::clean_names() %>% 
  rename_all(~str_replace(., "x2_5", "lcl") %>% 
               str_replace(., "x97_5", "ucl") %>% 
               str_replace(., "p_val", "p") %>% 
               str_replace(., "var_names", "var"))

# Get new column names in a dataframe
new_colnames <- tibble(cnms = colnames(sol_1)) %>% 
  mutate(variable = str_split_fixed(cnms, "_", n = 2)[,1], 
         group = str_split_fixed(cnms, "_", n = 2)[,2] %>% 
           if_else(. == "", "0", .) %>% 
           as.numeric)

# create a list of colnames by column group
cnms_bygroup <- split(new_colnames, new_colnames$group)


# create a list of sol result by column groups
sol_2 <- cnms_bygroup %>% 
  modify(~dplyr::select(sol_1, all_of(.$cnms)))

# rename all cols to match names, then cbind
sol_3 <- sol_2 %>% 
  modify(~rename_all(., ~str_split_fixed(., "_", 2)[,1])) %>% 
  bind_rows() %>% 
  filter(!is.na(metabolite))

# Separate beta and pips
sol_4 <- sol_3 %>% 
  mutate(effect = case_when(str_detect(var, "beta") ~ "beta", 
                            str_detect(var, "gamma") ~ "pip",
                            str_detect(var, "psi") ~ "beta", 
                            str_detect(var, "eta") ~ "eta"), 
         var = str_remove(var, ".beta") %>% 
           str_remove(".gamma") %>% 
           str_remove("eta.") %>%
           str_replace("psi", "mixture")) %>% 
  rename(exposure = var, 
         estimate = mean, 
         p_value = p) %>% 
  dplyr::select(exposure, effect, everything())


# Select Profiles data (for later? )
eta_profiles <- sol_4 %>% 
  filter(exposure == "high" | exposure == "low") %>% 
  dplyr::select(-p_value)

# Pivot wider
sol_5 <- sol_4 %>% 
  filter(exposure != "high", exposure != "low") %>%
  pivot_wider(id_cols = c("metabolite", "exposure"), names_from = effect, 
              values_from = c(estimate, sd, lcl, ucl, p_value)) %>% 
  dplyr::select(metabolite, exposure, contains("beta"), contains("pip")) %>% 
  dplyr::select(-p_value_pip) %>% 
  dplyr::rename(p_value = p_value_beta)


# Calculate new p values
sol_6 <- sol_5 %>% 
  mutate(wald = (estimate_beta/sd_beta), 
         # p = 2*(1-pnorm(abs(wald),0,1)),
         p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
         p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
         neg_log_p = -log10(p_value)) %>% 
  group_by(exposure) %>% 
  mutate(q_value = p.adjust(p_value, method = "fdr"), 
         significance = if_else(p_value < 0.05, "p < 0.05", "Not Sig."), 
         significancefdr = if_else(q_value < 0.05, "q < 0.05", "Not Sig.")) %>% 
  ungroup()

# Join with ft_metadata
sol_final <- left_join(ft_metadata, sol_6, by = c("feature" = "metabolite"))


# Remove PFAS from data: 
sol_final <- sol_final %>% 
  tidylog::filter(feature != "412.9662665_238.1521714",    # PFOA,  C18 neg
                  feature != "499.9333381_260.9808557",    # PFOS,  C18 neg
                  feature != "398.9360546_243.5838959",    # PFHxS, C18 neg
                  feature != "398.936398885251_33.3622678197535", # PFHxS, HELIC Neg
                  feature != "462.962944_247.5910964",     # PFNA
                  feature != "448.9331458_253.4631484",    # PFHpS
                  feature != "412.966235310664_35.0213619776204") # PFOA, HELIC Neg 



# Save 
write_csv(sol_final, 
          file = fs::path(dir_results, 
                          'PFAS_Mixtures', 
                          "sol_pfas_mixtures_results_w_09.csv"))

rm(cnms_bygroup, new_colnames, row1, col_names,
   sol_1, sol_2,sol_3, sol_4,sol_5, sol_6, sol_final, sol_og, sol_og_for_colnames)

# CHS -------------------------------------------------------------------
# Read in data, without headers.
chs_og_for_colnames <- read.table(fs::path(dir_results, 
                                           "PFAS_mixtures",
                                           "from_hpc",
                                           "chs_pfas_mixtures_mwas_w_09.csv"),
                                  sep = ",", na.strings = "",
                                  as.is = TRUE, fill = TRUE, header = FALSE)

# Get temp col names. Col names are not indexed correctly though- they 
# Need to be shifted to the right by 1, and we need to 
row1 <- chs_og_for_colnames[1,] %>% as.character(.)
col_names <- c("beta",
               row1[-length(row1)], 
               "var.names.9999", 
               "Mean.9999", 
               "SD.9999",         
               "X2.5..9999",      
               "X97.5..9999",    
               "p.val.9999",      
               "metabolite.9999")

# Read in final data with headers
chs_og <- read.table(fs::path(dir_results, 
                              "PFAS_mixtures",
                              "from_hpc",
                              "chs_pfas_mixtures_mwas_w_09.csv"),
                     sep = ",",
                     col.names = col_names,
                     na.strings = "",
                     as.is = TRUE, fill = TRUE, header = TRUE) %>% 
  dplyr::select(-beta)


# Clean Col names
chs_1 <- chs_og %>% 
  janitor::clean_names() %>% 
  rename_all(~str_replace(., "x2_5", "lcl") %>% 
               str_replace(., "x97_5", "ucl") %>% 
               str_replace(., "p_val", "p") %>% 
               str_replace(., "var_names", "var"))

# Get new column names in a dataframe
new_colnames <- tibble(cnms = colnames(chs_1)) %>% 
  mutate(variable = str_split_fixed(cnms, "_", n = 2)[,1], 
         group = str_split_fixed(cnms, "_", n = 2)[,2] %>% 
           if_else(. == "", "0", .) %>% 
           as.numeric)

# create a list of colnames by column group
cnms_bygroup <- split(new_colnames, new_colnames$group)


# create a list of chs result by column groups
chs_2 <- cnms_bygroup %>% 
  modify(~dplyr::select(chs_1, all_of(.$cnms)))

# rename all cols to match names, then cbind
chs_3 <- chs_2 %>% 
  modify(~rename_all(., ~str_split_fixed(., "_", 2)[,1])) %>% 
  bind_rows() %>% 
  filter(!is.na(metabolite))

# Separate beta and pips
chs_4 <- chs_3 %>% 
  mutate(effect = case_when(str_detect(var, "beta") ~ "beta", 
                            str_detect(var, "gamma") ~ "pip",
                            str_detect(var, "psi") ~ "beta", 
                            str_detect(var, "eta") ~ "eta"), 
         var = str_remove(var, ".beta") %>% 
           str_remove(".gamma") %>% 
           str_remove("eta.") %>%
           str_replace("psi", "mixture")) %>% 
  rename(exposure = var, 
         estimate = mean, 
         p_value = p) %>% 
  dplyr::select(exposure, effect, everything())

# Select Profiles data (for later? )
eta_profiles <- chs_4 %>% 
  filter(exposure == "high" | exposure == "low") %>% 
  dplyr::select(-p_value)

# Pivot wider
chs_5 <- chs_4 %>% 
  filter(exposure != "high", exposure != "low") %>%
  pivot_wider(id_cols = c("metabolite", "exposure"), names_from = effect, 
              values_from = c(estimate, sd, lcl, ucl, p_value)) %>% 
  dplyr::select(metabolite, exposure, contains("beta"), contains("pip")) %>% 
  dplyr::select(-p_value_pip) %>% 
  dplyr::rename(p_value = p_value_beta)


# Calculate new p values
chs_6 <- chs_5 %>% 
  mutate(wald = (estimate_beta/sd_beta), 
         # p = 2*(1-pnorm(abs(wald),0,1)),
         p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
         p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
         neg_log_p = -log10(p_value)) %>% 
  group_by(exposure) %>% 
  mutate(q_value = p.adjust(p_value, method = "fdr"), 
         significance = if_else(p_value < 0.05, "p < 0.05", "Not Sig."), 
         significancefdr = if_else(q_value < 0.05, "q < 0.05", "Not Sig.")) %>% 
  ungroup()

# Join with ft_metadata
chs_final <- left_join(ft_metadata, chs_6, by = c("feature" = "metabolite"))


# Remove PFAS from data: 
chs_final <- chs_final %>% 
  tidylog::filter(feature != "412.9662665_238.1521714",    # PFOA,  C18 neg
                  feature != "499.9333381_260.9808557",    # PFOS,  C18 neg
                  feature != "398.9360546_243.5838959",    # PFHxS, C18 neg
                  feature != "398.936398885251_33.3622678197535", # PFHxS, HELIC Neg
                  feature != "462.962944_247.5910964",     # PFNA
                  feature != "448.9331458_253.4631484",    # PFHpS
                  feature != "412.966235310664_35.0213619776204") # PFOA, HELIC Neg 


# Save 
write_csv(chs_final, 
          file = fs::path(dir_results, 
                          'PFAS_Mixtures', 
                          "chs_pfas_mixtures_results_w_09.csv"))


rm(cnms_bygroup, new_colnames, row1, col_names,
   chs_1, chs_2,chs_3, chs_4, chs_final, chs_og, chs_og_for_colnames)
