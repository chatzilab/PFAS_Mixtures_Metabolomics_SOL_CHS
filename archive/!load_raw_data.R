library(janitor)
library(jag2)

# load formats for variables
source(here::here("0_0.1_format_vars_funs.R"))
source(here::here("!directories.R"))

# Load Feature Tables --------------------------------------

readRDS()dir_metabolomics


# Get mz RT metadata ----------------------
c18neg_metadata   <- c18neg   %>% select(mz:max_intensity) %>% mutate(mode = "c18neg")
c18pos_metadata   <- c18pos   %>% select(mz:max_intensity) %>% mutate(mode = "c18pos")
hilicneg_metadata <- hilicneg %>% select(mz:max_intensity) %>% mutate(mode = "hilicneg")
hilicpos_metadata <- hilicpos %>% select(mz:max_intensity) %>% mutate(mode = "hilicpos")

# combine
feature_metadata <- bind_rows(c18neg_metadata,  
                              c18pos_metadata,  
                              hilicneg_metadata,
                              hilicpos_metadata)

# clean working environment
rm(c18neg_metadata, c18pos_metadata, hilicneg_metadata, hilicpos_metadata)


# Get clean feature tables (mz, rt, intensity) ----------------------
c18neg_ft   <- c18neg   %>% select(-(mz_min:max_intensity)) %>% transpose_ft()
c18pos_ft   <- c18pos   %>% select(-(mz_min:max_intensity)) %>% transpose_ft()
hilicneg_ft <- hilicneg %>% select(-(mz_min:max_intensity)) %>% transpose_ft()
hilicpos_ft <- hilicpos %>% select(-(mz_min:max_intensity)) %>% transpose_ft()


# Load mapfiles --------------------------------------
## C18 neg ---------------
c18neg_map <- read.delim(
  fs::path(dir_data, 
           "1_Mapfile",
           "HRE0026_C18neg_sample_mapfile_V2.txt")) %>% 
  janitor::clean_names() %>% 
  mutate(mode = "c18neg", order = row_number()) %>% 
  select(mode, everything())

## C18 pos ---------------
c18pos_map <- read.delim(
  fs::path(dir_data, 
           "1_Mapfile",
           "HRE0026_C18pos_sample_mapfile_V2.txt")) %>% 
  janitor::clean_names() %>% 
  mutate(mode = "c18pos", order = row_number()) %>% 
  select(mode, everything())

## HILIC neg  ---------------
hilicneg_map <- read.delim(
  fs::path(dir_data, 
           "1_Mapfile",
           "HRE0026_HILICneg_sample_mapfile_V2.txt")) %>% 
  janitor::clean_names() %>% 
  mutate(mode = "hilicneg", order = row_number()) %>% 
  select(mode, everything())

## HILIH pos ---------------
hilicpos_map <- read.delim(
  fs::path(dir_data, 
           "1_Mapfile",
           "HRE0026_HILICpos_sample_mapfile_V2.txt")) %>% 
  janitor::clean_names() %>% 
  mutate(mode = "hilicpos", order = row_number()) %>% 
  select(mode, everything())



# Join all map files ------------
map_files <- bind_rows(c18pos_map, c18neg_map, 
                       hilicneg_map, hilicpos_map) %>% 
  as_tibble() %>% 
  rename(lab_id = sample_id)   %>%
  mutate(lab_id = str_sub(lab_id, end = -3) %>% 
           str_remove_all(., "_02")) 

# Clean Environment
rm(c18neg_map, c18pos_map, hilicneg_map, hilicpos_map)


## Load MetaChem Sample Manifest --------------
chs_samplemanifest <- readxl::read_xlsx(
  fs::path(dir_data,
           "Pool_sample_manifest_MetaChem_02022021_SX.xlsx")) %>% 
  clean_names() %>% 
  select(id_on_label, labno, id) %>% 
  rename(lab_id = id_on_label)


# Get final MAP data Merge map_files with chs_samplemanifest and create study_id variable
sample_metadata <- left_join(map_files, chs_samplemanifest) %>% 
  mutate(id = case_when(sample_type == "SOLAR" ~ lab_id, 
                        sample_type == "CHS" ~ id, 
                        TRUE ~ sample_class), 
         file_name = tolower(file_name))


# Merge mapfiles and feature tables ---------------------------
c18neg_ft_2 <- sample_metadata %>% 
  select(mode, file_name, id, batch, order) %>% 
  filter(mode == "c18neg") %>%
  left_join(c18neg_ft ) %>% 
  select(-mode, -file_name)

c18pos_ft_2 <- sample_metadata %>% 
  select(mode, file_name, id, batch, order) %>% 
  filter(mode == "c18pos") %>%
  left_join(c18pos_ft ) %>% 
  select(-mode, -file_name)

hilicneg_ft_2 <- sample_metadata %>% 
  select(mode, file_name, id, batch, order) %>% 
  filter(mode == "hilicneg") %>%
  left_join(hilicneg_ft ) %>% 
  select(-mode, -file_name)

hilicpos_ft_2 <- sample_metadata %>% 
  select(mode, file_name, id, batch, order) %>% 
  filter(mode == "hilicpos") %>%
  left_join(hilicpos_ft ) %>% 
  select(-mode, -file_name)



# clean environment
# rm("c18neg", "c18neg_ft",
#    "c18pos", "c18pos_ft",
#    "chs_samplemanifest",
#    "hilicneg", "hilicneg_ft",
#    "hilicpos", "hilicpos_ft",
#    "map_files",
#    "transpose_ft")
# 
c18neg_ft   <- c18neg_ft_2  #; rm(c18neg_ft_2)
c18pos_ft   <- c18pos_ft_2  #; rm(c18pos_ft_2)
hilicneg_ft <- hilicneg_ft_2#; rm(hilicneg_ft_2)
hilicpos_ft <- hilicpos_ft_2#; rm(hilicpos_ft_2)
