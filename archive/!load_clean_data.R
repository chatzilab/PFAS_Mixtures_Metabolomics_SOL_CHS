
# load formats for variables
source(here::here("0_0.1_format_vars_funs.R"))
# source(here::here("!directories.R"))

# load separate data sets

##  SOLAR
solar_c18 <- read_csv(fs::path(dir_solar_metabolomics, 
                               "solar_c18_final_ft.csv"))

solar_hilic <- read_csv(fs::path(dir_solar_metabolomics, 
                                 "solar_hilic_final_ft.csv")) 

# solar_metab <- list(
#   metab = left_join(solar_c18, solar_hilic),
#   hilic_metab = colnames(solar_hilic %>% select(-id)), 
#   c18_metab = colnames(solar_c18 %>% select(-id))
# )

##  CHS
chs_c18 <- read_csv(fs::path(dir_chs_metabolomics, 
                             "metaair_c18_final_ft.csv"))

chs_hilic <- read_csv(fs::path(dir_chs_metabolomics,
                               "metaair_hilic_final_ft.csv"))

# Read in Exposure/outcome data
source(here::here("0_1.2_Exposure_outcome_data_cleaning.r"))



# Read in CV (RSD) data 
## Solar C18
solar_c18_rsd <- read_csv(fs::path(dir_solar_metabolomics, 
                                   "solar_c18_rsd.csv")) %>% 
  janitor::clean_names() %>% 
  rename(metabolite = x1) %>%
  mutate(metabolite = str_remove_all(metabolite, "X")) %>% 
  filter(metabolite %in% colnames(solar_c18))

## Solar HILIC
solar_hilic_rsd <- read_csv(fs::path(dir_solar_metabolomics, 
                                     "solar_hilic_rsd.csv")) %>% 
  janitor::clean_names() %>% 
  rename(metabolite = x1) %>%
  mutate(metabolite = str_remove_all(metabolite, "X")) %>% 
  filter(metabolite %in% colnames(solar_hilic))

## CHS C18
chs_c18_rsd <- read_csv(fs::path(dir_chs_metabolomics, 
                                 "metaair_c18_rsd.csv")) %>% 
  janitor::clean_names() %>% 
  rename(metabolite = x1) %>%
  mutate(metabolite = str_remove_all(metabolite, "X")) %>% 
  filter(metabolite %in% colnames(chs_c18))

## CHS HILIC
chs_hilic_rsd <- read_csv(fs::path(dir_chs_metabolomics, 
                                   "metaair_hilic_rsd.csv")) %>% 
  janitor::clean_names() %>% 
  rename(metabolite = x1) %>%
  mutate(metabolite = str_remove_all(metabolite, "X")) %>% 
  filter(metabolite %in% colnames(chs_hilic))