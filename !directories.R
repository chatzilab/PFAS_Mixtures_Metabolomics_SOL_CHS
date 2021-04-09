# directory fs::paths for file architecture
# ignore file on github
library(fs)
# home directory for project
dir_home <- here::here()

# SOLAR data folder on secure server
dir_solar_data_secure <- fs::path("R:/SOLAR EDCs")
dir_solar_data_secure_temp <- fs::path("C:/Users/jagoodri/Desktop", 
                                   "temp_lcms_cleaning/13.Untargeted Metabolomics 2021")
                              
# CHS data folder on secure server
# dir_chs_data_secure <- fs::path("R:/MetaChem")
dir_chs_data_secure <- fs::path("C:/Users/jagoodri/Desktop/temp_lcms_cleaning/")

# SOLAR Metabolomics Data Folder
dir_data_solar <- fs::path(paste0(dir_solar_data_secure,"/13.Untargeted Metabolomics 2021"))

# CHS Metabolomics Data Folder
dir_data_chs <- fs::path(paste0(dir_chs_data_secure,"/CHS Metabolomics Data"))

# # reports folder
dir_report <- fs::path(dir_home, "3 Reports")



# 
# # reports folder
# dir_program <- paste0("/Users/domusc/Documents/GitHub/multiomicsnafld")
# 

# 
# 
# # discussion folder
# dir_discuss <- paste0(dir_home, "/4 Discussions")
# 
