# directory fs::paths for file architecture
# ignore file on github
library(fs)

# home directory for project
dir_home <- dirname(here::here())

# Project Directory
dir_projects <- dirname(dir_home)

# data folder
dir_data <- fs::path("G:", 
                             "My Drive", 
                             "SOL CHS PFAS Metabolomics", 
                             "0_Data_mirror_do_not_edit")

# Temp results folder
dir_temp <- fs::path(dir_home, "1_Code", "Temporary results")

# reports folder
dir_reports <- fs::path(dir_home, "2_Reports")

# SOLAR data folder on secure server 
dir_solar_data_secure <- fs::path("G:", 
                                  "My Drive", 
                                  "SOLAR CHS PFAS Metabolomics Data", 
                                  "0_Data_mirror_do_not_edit")

# CHS data folder on secure server
dir_chs_data_secure <- fs::path("R:/MetaChem")

# SOLAR Metabolomics Data Folder
# dir_data_solar <- 

# CHS Metabolomics Data Folder
# dir_data_chs <- 

# discussion folder
# dir_discuss <- paste0(dir_home, "/4 Discussions")

# devtools::install_github("JAGoodrich/jag2")

