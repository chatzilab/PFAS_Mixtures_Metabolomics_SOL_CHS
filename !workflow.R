# Workflow for analysis overall
# 0) Exploratory Data Analysis
# 1)  Metabolome Wide Association Study (MWAS)
# 2) Targeted Selection
# 3) Sensitivity Analysis
# 4) Pathway Analysis
# 5) Integrated Analysis

# packages for creating html reports
library(knitr)
library(markdown)
library(rmarkdown)
source("!directories.R") 


# 0) Data Cleaning and QC's (Jesse's responsibilbity)
rmarkdown::render(input = paste0(dir_program,"/0_exploratory_data_analysis.Rmd"),
                  output_format = "html_document",
                  output_dir = dir_report)


# 1) Metabolome Wide Association Study (MWAS):
# 1.1: SOLAR  
# 1.1.1.1: SOLAR, PFAS, C18
# 1.1.1.2: SOLAR, PFAS, HELIC 
# 1.1.2.1: SOLAR, Lipophilic Chem, C18 
# 1.1.2.2: SOLAR, Lipophilic Chem, Helic
# 1.2: CHS
# 1.2.1.1: CHS, PFAS, C18
# 1.2.1.2: CHS, PFAS, HELIC
# 1.2.2.1: CHS, Lipophilic Chem, C18
# 1.2.2.2: CHS, Lipophilic Chem, Helic

# 2) Mummichog:
# 2.1: SOLAR
# 2.1.1.1: SOLAR, PFAS, C18
# 2.1.1.2: SOLAR, PFAS, HELIC
# 2.1.2.1: SOLAR, Lipophilic Chem, C18
# 2.1.2.2: SOLAR, Lipophilic Chem, Helic
# 2.2: CHS
# 2.2.1.1: CHS, PFAS, C18
# 2.2.1.2: CHS, PFAS, HELIC
# 2.2.2.1: CHS, Lipophilic Chem, C18
# 2.2.2.2: CHS, Lipophilic Chem, Helic



# 2) Targeted Selection
rmarkdown::render(input = paste0(dir_program,"/2_targeted_selection.Rmd"),
                  output_format = "html_document",
                  output_dir = dir_report)

# 5) LUCIDus Analysis
rmarkdown::render(input = paste0(dir_program,"/5_lucid_analysis.Rmd"),
                  output_format = "html_document",
                  output_dir = dir_report)
