# directory fs::paths for file architecture
# ignore file on github
library(fs)
# home directory for project
dir_home <- dirname(here::here())

dir_projects <- dirname(dir_home)

# SOLAR data folder on secure server 
dir_solar_data_secure <- fs::path("R:/SOLAR EDCs")

# SOLAR metabolomics data
dir_solar_metabolomics <- fs::path(dir_projects, 
                                   "Metabolomics Data Cleaning SOLAR", 
                                   "2 Cleaned Data")

# CHS metabolomics data
dir_chs_metabolomics <- fs::path(dir_projects, 
                                   "Metabolomics Data Cleaning CHS", 
                                   "2 Cleaned Data")

# CHS data folder on secure server
dir_chs_data_secure <- fs::path("R:/MetaChem")

# SOLAR Metabolomics Data Folder
dir_data_solar <- fs::path(paste0(dir_solar_data_secure,"/13.Untargeted Metabolomics 2021"))

# CHS Metabolomics Data Folder
dir_data_chs <- fs::path(paste0(dir_chs_data_secure,"/CHS Metabolomics Data"))

# reports folder
dir_report <- fs::path(dir_home, "3 Reports")

# discussion folder
# dir_discuss <- paste0(dir_home, "/4 Discussions")
 
