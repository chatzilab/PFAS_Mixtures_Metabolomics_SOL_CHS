# Set up data_for_mixtures_analysis on HPC. 
library("tidyverse")
library("purrr")

rm(list = ls())
source(here::here("!directories.R"))
source(here::here("!set_exposure_outcome_vars.R"))


# ALL PFAS ---------------------------
## SOLAR ------------------------------
rm(list = ls())
source(here::here("!set_exposure_outcome_vars.R"))
mixtures_name = "all_pfas"
mixtures_components = exposures_continuous
source(here::here("1_mixtures_analysis",
                  "0_1_SOLAR Set up datasets for mixtures analysis.R"))

## CHS ------------------------------
rm(list = ls())
source(here::here("!set_exposure_outcome_vars.R"))
mixtures_name = "all_pfas"
mixtures_components = exposures_continuous
source(here::here("1_mixtures_analysis",
                  "0_2_CHS Set up datasets for mixtures analysis.R"))


# PFSAs --------------------------------
## SOLAR ---------------------------
rm(list = ls())
source(here::here("!set_exposure_outcome_vars.R"))
mixtures_name = "pfsas"
mixtures_components = c("pfos", "pfhxs")
source(here::here("1_mixtures_analysis",
                  "0_1_SOLAR Set up datasets for mixtures analysis.R"))

## CHS ------------------------------
rm(list = ls())
source(here::here("!set_exposure_outcome_vars.R"))
mixtures_name = "pfsas"
mixtures_components = c("pfos", "pfhxs")
source(here::here("1_mixtures_analysis",
                  "0_2_CHS Set up datasets for mixtures analysis.R"))

# PFCAs --------------------------------
## SOLAR ---------------------------
rm(list = ls())
source(here::here("!set_exposure_outcome_vars.R"))
mixtures_name = "pfcas"
mixtures_components = c("pfda",  "pfhps", "pfna",  "pfoa")
source(here::here("1_mixtures_analysis",
                  "0_1_SOLAR Set up datasets for mixtures analysis.R"))

## CHS ------------------------------
rm(list = ls())
source(here::here("!set_exposure_outcome_vars.R"))
mixtures_name = "pfcas"
mixtures_components = c("pfda",  "pfhps", "pfna",  "pfoa")
source(here::here("1_mixtures_analysis",
                  "0_2_CHS Set up datasets for mixtures analysis.R"))
