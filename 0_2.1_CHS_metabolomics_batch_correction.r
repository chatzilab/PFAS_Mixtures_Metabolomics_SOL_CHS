# 0.2 CHS Metabolomics Data Cleaning
source(here::here("!directories.R"))
source(here::here("0.0.2_data_cleaning_funs.R"))

## Finalized Batch Correction Analysis
library(statTarget)
library(tidyverse)
library(data.table)
library(readxl)
library(gridExtra)
library(janitor)

# Read in Raw C18 and hilic LCMS Data
chs_c18 <- read_excel(path(dir_data_chs, 
                           "metabolomics_pre_batch_corrected", 
                           "metaair_c18_raw_ft.xlsx"))

chs_hilic <- read_excel(path(dir_data_chs, 
                             "metabolomics_pre_batch_corrected", 
                             "metaair_hilic_raw_ft.xlsx"))

# Read in Sample Metadata
chs_c18_metadata <- read_csv(path(dir_data_chs,  
                                  "metabolomics_pre_batch_corrected", 
                                  "c18 batch correction",
                                  "metaair_c18_sample_metadata_for_qcrf.csv"))

chs_hilic_metadata <- read_csv(path(dir_data_chs,  
                                    "metabolomics_pre_batch_corrected", 
                                    "hilic batch correction",
                                    "metaair_hilic_sample_metadata_for_qcrf.csv"))

# prep hilic data for qcrf-sc ---------------------------------------------

## Transpose hilic data and mutate id's to match metadata
chs_hilic_pre = transpose_ft(chs_hilic) %>% 
  mutate(sample = tolower(sample) %>% 
           str_replace(., "-", "_") %>%  
           str_sub(., end = -3),
         sample = if_else(str_detect(sample, "q3"), 
                          str_c(sample, "_QC", sep = ""), 
                          sample))

##  Merge with metadata
chs_hilic_pre <- left_join(chs_hilic_metadata, chs_hilic_pre, by = "sample")

# Scale QC Samples to the mean batch intensity 
##  Select QC samples from original data
chs_hilic_qc <- chs_hilic_pre %>%
  filter(str_detect(sample, "q3")) %>% 
  select(sample, batch, everything(), -class, -order)

##  Get mz_rt for later
chs_hilic_qc_og_names <- colnames(chs_hilic_qc)

##  Remove sample var and change to data.table
chs_hilic_qc_samples_overall_mean1 <- chs_hilic_qc %>% 
  select(-sample) %>% 
  as.data.table()

##  Calculate mean intensity for each feature in QC samples by batch
chs_hilic_qc_samples_overall_mean2 <- chs_hilic_qc_samples_overall_mean1[, lapply(.SD,mean), by=batch] 

##  Add important Vars
chs_hilic_qc_samples_overall_mean2 <- chs_hilic_qc_samples_overall_mean2 %>% 
  mutate(samples = "QC")

# Now, do the same with samples:
##  With Samples: Remove QC Samples and change to data.table
chs_hilic_pre1 <- chs_hilic_pre %>%
  filter(str_detect(sample, "q3", negate = TRUE)) %>% 
  select(batch, everything(),-sample, -class, -order) %>% 
  as.data.table()

##  With Samples: Calculate batch mean intensity for each feature in all samples (except QC) 
chs_hilic_features_overall_mean <- chs_hilic_pre1[, lapply(.SD,mean), by=batch] 

##  With Samples: Add important vars
chs_hilic_features_overall_mean <- chs_hilic_features_overall_mean %>%
  mutate(samples = "All Samples") %>% 
  select(batch, samples, everything())

# Merge qc batch means and sample batch means  
chs_hilic_batch_means <- bind_rows(chs_hilic_features_overall_mean, chs_hilic_qc_samples_overall_mean2)

##  Calculate difference between qc batch means and batch means for all samples
chs_hilic_batch_deltas <- chs_hilic_batch_means %>%
  group_by(batch) %>%
  arrange(samples) %>%
  summarise(across(where(is.numeric), function(x){x[1] - x[2]}))

## Select qc sample names and batch
chs_hilic_qc_ids <- chs_hilic_qc %>% select(sample, batch) 

##  Join batch deltas with ids (expands table)
chs_hilic_batch_delta_for_id <- left_join(chs_hilic_qc_ids, chs_hilic_batch_deltas)

## Bind QC sample values with Batch Deltas to be able to take sum later
chs_hilic_temp_delta_data <- bind_rows(chs_hilic_qc, chs_hilic_batch_delta_for_id, by = NULL)

##  Correct all QC samples by average intensity of batch
chs_hilic_qc_corrected <- chs_hilic_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

##############################################################################

# Center QC samples so that min value is the same as min from og data 
## Get minimum values (excluding 0) for each feature from the corrected QC data 
min_values_qc <- chs_hilic_qc_corrected %>%
  select(-c(sample, batch)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then get min value
  mutate(value = "qc_min") 

## Get minimum values (excluding 0) for each feature from the original data 
min_values <- chs_hilic_pre %>% 
  select(-c(sample, batch, class, order)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then get min value
  mutate(value = "og_min") 

##  Calculate difference between min qc and min original
min_values1 <- bind_rows(min_values_qc, min_values) %>% 
  # select(value, everything()) %>% 
  summarise(across(where(is.numeric), ~(.[2]-.[1])))

##  Cross join to get full data for second sum
min_values2 <- left_join(chs_hilic_qc_ids, min_values1, by = character())

##  Add difference between minimum QC value and original minimum value (prevents negative qc values)
chs_hilic_qc_corrected <- chs_hilic_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

chs_hilic_qc_corrected_by_minimum <- bind_rows(chs_hilic_qc_corrected, 
                                                 min_values2, 
                                                 by = NULL)  %>% 
  group_by(sample, batch) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 
##############################################################################


# Merge everything back together 
##  Get QC sample names
sample <- chs_hilic_qc_corrected_by_minimum$sample

## Transpose Corrected QC samples
chs_hilic_qc_corrected_t <- chs_hilic_qc_corrected_by_minimum %>%
  select(-sample, -batch) %>%
  t() %>%
  as.data.frame()

##  Rename Columns with correct sample names
colnames(chs_hilic_qc_corrected_t) <- sample

##  Rownames to variable
chs_hilic_qc_corrected_t$name <- rownames(chs_hilic_qc_corrected_t)
rownames(chs_hilic_qc_corrected_t) <- NULL
chs_hilic_qc_corrected_t1 <- chs_hilic_qc_corrected_t %>% select(name, everything())

# Merge qc_corrected with Original Data
##  Mutate chs_hilic to match metadata names (removes all of the "_1" at the end of sample names)
##  Also, remove QC saples (corrected samples are merged in next step)
chs_hilic_pre_samples <- chs_hilic %>%
  mutate(name = str_c(mz, time, sep = "_")) %>%
  select(name, everything(), -mz, -time) %>% 
  clean_names() %>% 
  rename_all(~ str_sub(.x, end = -3)) %>% 
  rename(name = na) %>% 
  select(!contains("q3"))

##  Merge Samples and Batch Corrected QC's
chs_hilic_pre_batch_correction <- left_join(chs_hilic_qc_corrected_t1, 
                                            chs_hilic_pre_samples)

# Reorder
chs_hilic_pre_batch_correction <- chs_hilic_pre_batch_correction %>%
  select(name, all_of(chs_hilic_metadata$sample))

# Run CHS hilic QC-RFSC -------------------------------------------------------
##  Set working directory
setwd(path(dir_data_chs,
           "metabolomics_pre_batch_corrected", 
           "hilic batch correction"))

##  Save QC batch corrected LC-MS Data 
write_csv(chs_hilic_pre_batch_correction, "metaair_hilic_raw_ft_for_qcsc.csv")


##  Run Batch Correction
shiftCor(samPeno = "metaair_hilic_sample_metadata_for_qcrf.csv",
         samFile = "metaair_hilic_raw_ft_for_qcsc.csv",
         Frule = 0.8, # Removes features with >20% non-detects
         coCV = 30, # Removes features with QC sample CV >30
         MLmethod = 'QCRFSC',
         QCspan = 0,
         imputeM = "minHalf",
         plot = FALSE)


# prep c18 data for qcrf-sc ---------------------------------------------

## Transpose c18 data and mutate id's to match metadata
chs_c18_pre = transpose_ft(chs_c18) %>% 
  mutate(sample = tolower(sample) %>% 
           str_replace(., "-", "_"),
         sample = if_else(str_detect(sample, "q3"), 
                          str_c(sample, "_QC", sep = ""), 
                          sample))

##  Merge with metadata
chs_c18_pre <- left_join(chs_c18_metadata, chs_c18_pre, by = "sample")

# Scale QC Samples to the mean batch intensity 
##  Select QC samples from original data
chs_c18_qc <- chs_c18_pre %>%
  filter(str_detect(sample, "q3")) %>% 
  select(sample, batch, everything(), -class, -order)

##  Get mz_rt for later
chs_c18_qc_og_names <- colnames(chs_c18_qc)

##  Remove sample var and change to data.table
chs_c18_qc_samples_overall_mean1 <- chs_c18_qc %>% 
  select(-sample) %>% 
  as.data.table()

##  Calculate mean intensity for each feature in QC samples by batch
chs_c18_qc_samples_overall_mean2 <- chs_c18_qc_samples_overall_mean1[, lapply(.SD,mean), by=batch] 

##  Add important Vars
chs_c18_qc_samples_overall_mean2 <- chs_c18_qc_samples_overall_mean2 %>% 
  mutate(samples = "QC")

# Now, do the same with samples:
##  With Samples: Remove QC Samples and change to data.table
chs_c18_pre1 <- chs_c18_pre %>%
  filter(str_detect(sample, "q3", negate = TRUE)) %>% 
  select(batch, everything(),-sample, -class, -order) %>% 
  as.data.table()

##  With Samples: Calculate batch mean intensity for each feature in all samples (except QC) 
chs_c18_features_overall_mean <- chs_c18_pre1[, lapply(.SD,mean), by=batch] 

##  With Samples: Add important vars
chs_c18_features_overall_mean <- chs_c18_features_overall_mean %>%
  mutate(samples = "All Samples") %>% 
  select(batch, samples, everything())

# Merge qc batch means and sample batch means  
chs_c18_batch_means <- bind_rows(chs_c18_features_overall_mean, chs_c18_qc_samples_overall_mean2)

##  Calculate difference between qc batch means and batch means for all samples
chs_c18_batch_deltas <- chs_c18_batch_means %>%
  group_by(batch) %>%
  arrange(samples) %>%
  summarise(across(where(is.numeric), function(x){x[1] - x[2]}))

## Select qc sample names and batch
chs_c18_qc_ids <- chs_c18_qc %>% select(sample, batch) 

##  Join batch deltas with ids (expands table)
chs_c18_batch_delta_for_id <- left_join(chs_c18_qc_ids, chs_c18_batch_deltas)

## Bind QC sample values with Batch Deltas to be able to take sum later
chs_c18_temp_delta_data <- bind_rows(chs_c18_qc, chs_c18_batch_delta_for_id, by = NULL)

##  Correct all QC samples by average intensity of batch
chs_c18_qc_corrected <- chs_c18_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

##############################################################################

# Center QC samples so that min value is the same as min from og data 
## Get minimum values (excluding 0) for each feature from the corrected QC data 
min_values_qc <- chs_c18_qc_corrected %>%
  select(-c(sample, batch)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then get min value
  mutate(value = "qc_min") 

## Get minimum values (excluding 0) for each feature from the original data 
min_values <- chs_c18_pre %>% 
  select(-c(sample, batch, class, order)) %>% 
  map_dfc(~na_if(., 0) %>% min(., na.rm = TRUE)) %>%   # Changes all 0's to NA then get min value
  mutate(value = "og_min") 

##  Calculate difference between min qc and min original
min_values1 <- bind_rows(min_values_qc, min_values) %>% 
  # select(value, everything()) %>% 
  summarise(across(where(is.numeric), ~(.[2]-.[1])))

##  Cross join to get full data for second sum
min_values2 <- left_join(chs_c18_qc_ids, min_values1, by = character())

##  Add difference between minimum QC value and original minimum value (prevents negative qc values)
chs_c18_qc_corrected <- chs_c18_temp_delta_data %>%
  group_by(batch, sample) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 

chs_c18_qc_corrected_by_minimum <- bind_rows(chs_c18_qc_corrected, 
                                               min_values2, 
                                               by = NULL)  %>% 
  group_by(sample, batch) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() 
##############################################################################


# Merge everything back together 
##  Get QC sample names
sample <- chs_c18_qc_corrected$sample

## Transpose Corrected QC samples
chs_c18_qc_corrected_t <- chs_c18_qc_corrected %>%
  select(-sample, -batch) %>%
  t() %>%
  as.data.frame()

##  Rename Columns with correct sample names
colnames(chs_c18_qc_corrected_t) <- sample

##  Rownames to variable
chs_c18_qc_corrected_t$name <- rownames(chs_c18_qc_corrected_t)
rownames(chs_c18_qc_corrected_t) <- NULL
chs_c18_qc_corrected_t1 <- chs_c18_qc_corrected_t %>% select(name, everything())

# Merge qc_corrected with Original Data
##  Mutate chs_c18 to match metadata names 
##  Also, remove QC saples (corrected samples are merged in next step)
chs_c18_pre_samples <- chs_c18 %>%
  mutate(name = str_c(mz, time, sep = "_")) %>%
  select(name, everything(), -mz, -time) %>% 
  clean_names() %>% 
  select(!contains("q3"))

##  Merge Samples and Batch Corrected QC's
chs_c18_pre_batch_correction <- left_join(chs_c18_qc_corrected_t1, 
                                          chs_c18_pre_samples)

# Reorder
chs_c18_pre_batch_correction <- chs_c18_pre_batch_correction %>%
  select(name, all_of(chs_c18_metadata$sample))

# Run CHS c18 QC-RFSC -------------------------------------------------------
##  Set working directory
setwd(path(dir_data_chs,
           "metabolomics_pre_batch_corrected", 
           "c18 batch correction"))

##  Save QC batch corrected LC-MS Data 
write_csv(chs_c18_pre_batch_correction, "metaair_c18_raw_ft_for_qcsc.csv")


##  Run Batch Correction
shiftCor(samPeno = "metaair_c18_sample_metadata_for_qcrf.csv",
         samFile = "metaair_c18_raw_ft_for_qcsc.csv",
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

##  Load linking dataset
linker <- read_csv(path(dir_data_chs,
                        "metabolomics_pre_batch_corrected",
                        "metaair_sample_participant_id_link.csv"))

##  Load batch corrected c18 data 
chs_c18_post_correction <- read_csv(path(dir_data_chs, 
                                         "metabolomics_pre_batch_corrected", 
                                         "c18 batch correction",
                                         "statTarget", 
                                         "shiftCor", 
                                         "After_shiftCor", 
                                         "shift_sample_cor.csv"))
##  Log transform and scale features
chs_c18_post_correction_1 <- chs_c18_post_correction %>% 
  select(-class) %>% 
  mutate(across(where(is.numeric), ~log(.) %>% scale()))

##  Add CHS ID and remove sample ID
chs_c18_post_correction_2 <- chs_c18_post_correction_1 %>% 
  inner_join(linker, by = 'sample') %>%
  arrange(id) %>% 
  select(id, everything(), -sample)

##  Save C18 data
write_csv(chs_c18_post_correction_2, 
          path(dir_data_chs, 
               "metabolomics_data_cleaned",
               "metaair_c18_final_ft.csv"))

# hilic Data
##  Read in hilic Data
chs_hilic_post_correction <- read_csv(path(dir_data_chs, 
                                             "metabolomics_pre_batch_corrected", 
                                             "hilic batch correction",
                                             "statTarget", 
                                             "shiftCor", 
                                             "After_shiftCor", 
                                             "shift_sample_cor.csv"))

##  Log transform and scale features
chs_hilic_post_correction_1 <- chs_hilic_post_correction %>% 
  select(-class) %>% 
  mutate(across(where(is.numeric), ~log(.) %>% scale()))

##  Add CHS ID and remove sample ID
chs_hilic_post_correction_2 <- chs_hilic_post_correction_1 %>% 
  inner_join(linker, by = 'sample') %>%
  arrange(id) %>% 
  select(id, everything(), -sample)

# Save HILIC data 
write_csv(chs_hilic_post_correction_2, 
          path(dir_data_chs, 
               "metabolomics_data_cleaned",
               "metaair_hilic_final_ft.csv"))