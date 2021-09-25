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
dir_temp <- fs::path(dir_home, "1_Code", "Temp results")

# reports folder
dir_reports <- fs::path(dir_home, "2_Reports")
