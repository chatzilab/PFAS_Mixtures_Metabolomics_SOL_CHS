# 0.1 SOLAR Metabolomics Data Cleaning
setwd(here::here())
rm(list = ls())
source("!directories.R")
source("0.0.2_data_cleaning_funs.r")

## Finalized Batch Correction Analysis
library(statTarget)
library(tidyverse)
library(gridExtra)
library(janitor)
library(readxl)
library(data.table)

# Read in Raw C18 and hilic LCMS Data
solar_hilic <- read_csv(fs::path(dir_solar_data_secure_temp, 
                             "metabolomics_pre_batch_corrected", 
                             "solar_hilic_raw_ft.csv")) %>% 
  janitor::clean_names() %>% 
  rename_all(~str_replace_all(., "_qc", "_QC"))

solar_c18 <- read_csv(fs::path(dir_solar_data_secure_temp,
                           "metabolomics_pre_batch_corrected",
                           "solar_C18_raw_ft.csv"))

# Read in Sample Metadata
solar_hilic_metadata <- read_csv(fs::path(dir_solar_data_secure_temp,  
                                      "metabolomics_pre_batch_corrected", 
                                      "hilic batch correction",
                                      "solar_hilic_sample_metadata_for_qcrf.csv"))

solar_c18_metadata <- read_csv(fs::path(dir_solar_data_secure_temp,  
                                    "metabolomics_pre_batch_corrected", 
                                    "c18 batch correction",
                                    "solar_c18_sample_metadata_for_qcrf.csv"))


# prep hilic data for qcrf-sc  ---------------------------------------------
solar_hilic <- solar_hilic %>% 
  mutate(name = str_c(mz, time, sep = "_")) %>% 
  select(name, everything(), -mz, -time)
## Transpose hilic data and mutate id's to match metadata
solar_hilic_pre = transpose_ft(solar_hilic, pct_metab_na_allowed = 0.2, mzrt = FALSE) 


##  Merge with metadata
solar_hilic_pre <- left_join(solar_hilic_metadata, solar_hilic_pre, by = "sample")

# Scale QC Samples to the mean batch intensity 
##  Select QC samples from original data
solar_hilic_qc <- solar_hilic_pre %>%
  filter(str_detect(sample, "q3")) %>% 
  select(sample, batch, everything(), -class, -order)

##  Get mz_rt for later
solar_hilic_qc_og_names <- colnames(solar_hilic_qc)

##  Remove sample var and change to data.table
solar_hilic_qc_mean1 <- solar_hilic_qc %>% 
  select(-sample) %>% 
  as.data.table()

##  Calculate mean intensity for each feature in QC samples by batch
solar_hilic_qc_mean2 <- 
  solar_hilic_qc_mean1[ , lapply(.SD,mean), by=batch] 

##  Add important Vars
solar_hilic_qc_mean2 <- solar_hilic_qc_mean2 %>% 
  mutate(samples = "QC")

# Now, do the same with samples:
##  With Samples: Remove QC Samples and change to data.table
solar_hilic_pre1 <- solar_hilic_pre %>%
  filter(str_detect(sample, "q3", negate = TRUE)) %>% 
  select(batch, everything(),-sample, -class, -order) %>% 
  as.data.table()

##  With Samples: Calculate batch mean intensity for each feature in all samples (except QC) 
solar_hilic_features_mean <- solar_hilic_pre1[, lapply(.SD,mean), by=batch] 

##  With Samples: Add important vars
solar_hilic_features_mean <- solar_hilic_features_mean %>%
  mutate(samples = "All Samples") %>% 
  select(batch, samples, everything())

# Merge qc batch means and sample batch means  
solar_hilic_batch_means <- bind_rows(solar_hilic_features_mean, solar_hilic_qc_mean2)

##  Calculate difference between qc batch means and batch means for all samples
solar_hilic_batch_deltas <- solar_hilic_batch_means %>%
  group_by(batch) %>%
  arrange(samples) %>%
  summarise(across(where(is.numeric), function(x){x[1] - x[2]}))

## Select qc sample names and batch
solar_hilic_qc_ids <- solar_hilic_qc %>% 
  select(sample, batch) 

##  Join batch deltas with ids (expands table)
solar_hilic_batch_delta_for_id <- left_join(solar_hilic_qc_ids, solar_hilic_batch_deltas)

## Bind QC sample values with Batch Deltas to be able to take sum later
solar_hilic_temp_delta_data <- bind_rows(solar_hilic_qc, 
                                         solar_hilic_batch_delta_for_id, 
                                         by = NULL) 

##  Change 0's to na before calculating batch means
solar_hilic_temp_delta_data <- map_dfc(solar_hilic_temp_delta_data, 
                                       ~na_if(., 0))

##  Correct all QC samples by average intensity of batch
solar_hilic_qc_corrected <- solar_hilic_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

##  Get minimum values (excluding 0) for each feature from the corrected QC data 
##  Gives a warning: no non-missing arguments to min; returning Inf: THIS IS OK! 
min_values_qc <- solar_hilic_qc_corrected %>%
  select(-c(sample, batch)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then get min value
  mutate(value = "qc_min") 

## Get minimum values (excluding 0) for each feature from the original data 
min_values <- solar_hilic_pre %>% 
  select(-c(sample, batch, class, order)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then get min value
  mutate(value = "og_min") 

##  Calculate difference between min qc and min original
min_values1 <- bind_rows(min_values_qc, min_values) %>% 
  # select(value, everything()) %>% 
  summarise(across(where(is.numeric), ~(.[2]-.[1])))

##  Cross join to get full data for second sum
min_values2 <- left_join(solar_hilic_qc_ids, min_values1, by = character())

##  Add difference between minimum QC value and original minimum value (prevents negative qc values)
solar_hilic_qc_corrected <- solar_hilic_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

solar_hilic_qc_corrected_by_minimum <- bind_rows(solar_hilic_qc_corrected, 
                                                 min_values2, 
                                                 by = NULL)  %>% 
  group_by(sample, batch) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 


# Merge everything back together
##  Get QC sample names
sample <- solar_hilic_qc_corrected_by_minimum$sample

## Transpose Corrected QC samples
solar_hilic_qc_corrected_t <- solar_hilic_qc_corrected_by_minimum %>%
  select(-sample, -batch) %>%
  t() %>%
  as.data.frame()

##  Rename Columns with correct sample names
colnames(solar_hilic_qc_corrected_t) <- sample

##  Rownames to variable
solar_hilic_qc_corrected_t$name <- rownames(solar_hilic_qc_corrected_t)
rownames(solar_hilic_qc_corrected_t) <- NULL
solar_hilic_qc_corrected_t1 <- solar_hilic_qc_corrected_t %>% select(name, everything())

# Merge qc_corrected with Original Data
##  Also, remove QC saples (corrected samples are merged in next step)
solar_hilic_pre_samples <- solar_hilic %>%
  select(!contains("q3"))

##  Merge Samples and Batch Corrected QC's (left join also drops features which were only selected based on missingness )
solar_hilic_pre_batch_correction <- left_join(solar_hilic_qc_corrected_t1, 
                                              solar_hilic_pre_samples)

# Reorder
solar_hilic_pre_batch_correction <- solar_hilic_pre_batch_correction %>%
  select(name, all_of(solar_hilic_metadata$sample)) %>% 
  map_dfc(~replace_na(., 0))



# Run solar hilic QC-RFSC -------------------------------------------------------
##  Set working directory
setwd(fs::path(dir_solar_data_secure_temp,
           "metabolomics_pre_batch_corrected", 
           "hilic batch correction"))

##  Save QC batch corrected LC-MS Data 
write_csv(solar_hilic_pre_batch_correction, "solar_hilic_raw_ft_for_qcsc.csv")

rm(list = ls())
##  Run Batch Correction
shiftCor(samPeno = "solar_hilic_sample_metadata_for_qcrf.csv",
         samFile = "solar_hilic_raw_ft_for_qcsc.csv",
         Frule = 0.8, # Removes features with >20% non-detects
         coCV = 30, # Removes features with QC sample CV >30
         MLmethod = 'QCRFSC',
         QCspan = 0,
         imputeM = "minHalf",
         plot = FALSE)




# prep c18 data for qcrf-sc  ---------------------------------------------

## Transpose c18 data and mutate id's to match metadata
solar_c18_pre = transpose_ft(solar_c18, pct_metab_na_allowed = 0.2, mzrt = FALSE) 

##  Merge with metadata
solar_c18_pre <- left_join(solar_c18_metadata, solar_c18_pre, by = "sample")

# Scale QC Samples to the mean batch intensity 
##  Select QC samples from original data
solar_c18_qc <- solar_c18_pre %>%
  filter(str_detect(sample, "q3")) %>% 
  select(sample, batch, everything(), -class, -order)

##  Get mz_rt for later
solar_c18_qc_og_names <- colnames(solar_c18_qc)

##  Remove sample var and change to data.table
solar_c18_qc_mean1 <- solar_c18_qc %>% 
  select(-sample) %>% 
  as.data.table()

##  Calculate mean intensity for each feature in QC samples by batch
solar_c18_qc_mean2 <- 
  solar_c18_qc_mean1[ , lapply(.SD,mean), by=batch] 

##  Add important Vars
solar_c18_qc_mean2 <- solar_c18_qc_mean2 %>% 
  mutate(samples = "QC")

# Now, do the same with samples:
##  With Samples: Remove QC Samples and change to data.table
solar_c18_pre1 <- solar_c18_pre %>%
  filter(str_detect(sample, "q3", negate = TRUE)) %>% 
  select(batch, everything(),-sample, -class, -order) %>% 
  as.data.table()

##  With Samples: Calculate batch mean intensity for each feature in all samples (except QC) 
solar_c18_features_mean <- solar_c18_pre1[, lapply(.SD,mean), by=batch] 

##  With Samples: Add important vars
solar_c18_features_mean <- solar_c18_features_mean %>%
  mutate(samples = "All Samples") %>% 
  select(batch, samples, everything())

# Merge qc batch means and sample batch means  
solar_c18_batch_means <- bind_rows(solar_c18_features_mean, solar_c18_qc_mean2)

##  Calculate difference between qc batch means and batch means for all samples
solar_c18_batch_deltas <- solar_c18_batch_means %>%
  group_by(batch) %>%
  arrange(samples) %>%
  summarise(across(where(is.numeric), function(x){x[1] - x[2]}))

## Select qc sample names and batch
solar_c18_qc_ids <- solar_c18_qc %>% 
  select(sample, batch) 

##  Join batch deltas with ids (expands table)
solar_c18_batch_delta_for_id <- left_join(solar_c18_qc_ids, solar_c18_batch_deltas)

## Bind QC sample values with Batch Deltas to be able to take sum later
solar_c18_temp_delta_data <- bind_rows(solar_c18_qc, 
                                       solar_c18_batch_delta_for_id, 
                                       by = NULL) 

##  Change 0's to na before calculating batch means
solar_c18_temp_delta_data <- map_dfc(solar_c18_temp_delta_data, 
                                     ~na_if(., 0))

##  Correct all QC samples by average intensity of batch
solar_c18_qc_corrected <- solar_c18_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

## Get minimum values (excluding 0) for each feature from the corrected QC data 
min_values_qc <- solar_c18_qc_corrected %>%
  select(-c(sample, batch)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then get min value
  mutate(value = "qc_min") 

## Get minimum values (excluding 0) for each feature from the original data 
min_values <- solar_c18_pre %>% 
  select(-c(sample, batch, class, order)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then gets min value
  mutate(value = "og_min") 

##  Calculate difference between min qc and min original
min_values1 <- bind_rows(min_values_qc, min_values) %>% 
  summarise(across(where(is.numeric), ~(.[2]-.[1])))

##  Cross join to get full data for second sum
min_values2 <- left_join(solar_c18_qc_ids, min_values1, by = character())

##  Add difference between minimum QC value and original minimum value (prevents negative qc values)
solar_c18_qc_corrected <- solar_c18_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

##  Bind minimum values and qc values, and perform correction
solar_c18_qc_corrected_by_minimum <- bind_rows(solar_c18_qc_corrected, 
                                               min_values2, 
                                               by = NULL)  %>% 
  group_by(sample, batch) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

# Merge everything back together 
##  Get QC sample names
sample <- solar_c18_qc_corrected_by_minimum$sample

## Transpose Corrected QC samples
solar_c18_qc_corrected_t <- solar_c18_qc_corrected_by_minimum %>%
  select(-sample, -batch) %>%
  t() %>%
  as.data.frame()

##  Rename Columns with correct sample names
colnames(solar_c18_qc_corrected_t) <- sample

##  Rownames to variable
solar_c18_qc_corrected_t$name_1 <- rownames(solar_c18_qc_corrected_t)
rownames(solar_c18_qc_corrected_t) <- NULL
solar_c18_qc_corrected_t1 <- solar_c18_qc_corrected_t %>% select(name_1, everything())

# Merge qc_corrected with Original Data
##  Also, remove QC saples (corrected samples are merged in next step)
solar_c18_pre_samples <- solar_c18 %>%
  select(!contains("q3"))

##  Merge Samples and Batch Corrected QC's (left join also drops features which were only selected based on missingness )
solar_c18_pre_batch_correction <- left_join(solar_c18_qc_corrected_t1, 
                                            solar_c18_pre_samples)

# Reorder
solar_c18_pre_batch_correction <- solar_c18_pre_batch_correction %>%
  select(name_1, all_of(solar_c18_metadata$sample)) %>% 
  map_dfc(~na_if(., -Inf) %>% replace_na(., 0))

# Run solar c18 QC-RFSC -------------------------------------------------------
##  Set working directory
setwd(fs::path(dir_solar_data_secure_temp,
           "metabolomics_pre_batch_corrected", 
           "c18 batch correction"))

##  Save QC batch corrected LC-MS Data 
write_csv(solar_c18_pre_batch_correction, "solar_c18_raw_ft_for_qcsc.csv")

rm(list = ls())
##  Run Batch Correction
shiftCor(samPeno = "solar_c18_sample_metadata_for_qcrf.csv",
         samFile = "solar_c18_raw_ft_for_qcsc.csv",
         Frule = 0.8, # Removes features with >20% non-detects
         coCV = 30, # Removes features with QC sample CV >30
         MLmethod = 'QCRFSC',
         QCspan = 0,
         imputeM = "minHalf",
         plot = FALSE)


# Reset working directory and move corrected data to correct location ---------
setwd(here::here())
rm(list = ls())
source(here::here("!directories.R"))


##  Load batch corrected c18 data 
solar_c18_post_correction <- read_csv(fs::path(dir_solar_data_secure_temp, 
                                           "metabolomics_pre_batch_corrected", 
                                           "c18 batch correction",
                                           "statTarget", 
                                           "shiftCor", 
                                           "After_shiftCor", 
                                           "shift_sample_cor.csv"))

##  Log transform and scale features based on sample type (Class var)
solar_c18_post_correction_1 <- solar_c18_post_correction %>% 
  filter(str_detect(sample, "nist", negate = TRUE), 
         str_detect(sample, "_2", negate = TRUE)) %>% 
  group_by(class) %>% 
  mutate(across(where(is.numeric),
                ~log(.) %>% scale())) %>% 
  ungroup()

##  Add solar ID and remove sample ID
solar_c18_post_correction_2 <- solar_c18_post_correction_1 %>% 
  rename(id = sample) %>% 
  select(-class) %>% 
  mutate(id = str_sub(id, end = -3))

##  Save C18 data
write_csv(solar_c18_post_correction_2, 
          fs::path(dir_solar_data_secure_temp, 
               "metabolomics_data_cleaned",
               "solar_c18_final_ft.csv"))

# hilic Data
##  Read in hilic Data
solar_hilic_post_correction <- read_csv(fs::path(dir_solar_data_secure_temp, 
                                             "metabolomics_pre_batch_corrected", 
                                             "hilic batch correction",
                                             "statTarget", 
                                             "shiftCor", 
                                             "After_shiftCor", 
                                             "shift_sample_cor.csv"))

##  Log transform and scale features; remove sample which did not pass c18 QC's 
solar_hilic_post_correction_1 <- solar_hilic_post_correction %>% 
  filter(str_detect(sample, "nist", negate = TRUE), 
         str_detect(sample, "_2", negate = TRUE), 
         str_detect(sample, "00092", negate = TRUE)) %>% 
  rename(id = sample) %>% 
  group_by(class) %>% 
  mutate(across(where(is.numeric),
                ~log(.) %>% scale()), 
         id = str_sub(id, end = -3))

# Save HILIC data 
write_csv(solar_hilic_post_correction_1, 
          fs::path(dir_solar_data_secure_temp, 
               "metabolomics_data_cleaned",
               "solar_hilic_final_ft.csv"))
