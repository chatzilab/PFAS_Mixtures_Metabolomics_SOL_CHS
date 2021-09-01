# directory fs::paths for file architecture
# ignore file on github
library(fs)
# home directory for project
dir_home <- dirname(here::here())

# Project Directory
dir_projects <- dirname(dir_home)

# Metabolomics Data Cleaning Data Directory
dir_metabolomics <- fs::path(
  "C:/Users/jagoodri/Dropbox (USC Lab)/Project Directories", 
  "HRE Metabolomics Data Cleaning SOL CHS", 
  "3_clean_data")

# reports folder
dir_data <- fs::path(dir_home, "0_Data")

# SOLAR data folder on secure server 
dir_solar_data_secure <- fs::path("R:/SOLAR EDCs")

# CHS data folder on secure server
dir_chs_data_secure <- fs::path("R:/MetaChem")

# SOLAR Metabolomics Data Folder
# dir_data_solar <- 

# CHS Metabolomics Data Folder
# dir_data_chs <- 

# discussion folder
# dir_discuss <- paste0(dir_home, "/4 Discussions")

# devtools::install_github("JAGoodrich/jag2")

