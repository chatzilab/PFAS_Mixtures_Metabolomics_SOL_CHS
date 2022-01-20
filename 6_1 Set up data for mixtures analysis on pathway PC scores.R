# Set up data_for_mixtures_analysis on HPC
library("tidyverse")
library("purrr")

rm(list = ls())
source(here::here("!directories.R"))
source(here::here("!set_exposure_outcome_vars.R"))
source(here::here("0_0_1_format_vars_funs.R"))
source(here::here("!load_data.R"))

# Get BHRMA Function
source(here::here("1_mixtures_analysis",
                  "0_BHRMA.g_function.R"))

# Solar ------------------------------------------
# Read in PCA Data 
pca_dat <- read_rds(fs::path(dir_results_mum_mixtures,
                             "SOL CHS PC Scores of sig pathways.RDS"))

sol_pca_dat <- pca_dat$solar

# Create numberic Vars
solar_eo <- exposure_outcome$solar %>% 
  mutate(id = as.character(id), 
         sex.num = ifelse(sex == "Female",1,0),
         ses.num = recode(ses, "[3,11]" = 1, "(11,15.5]" = 2, 
                          "(15.5,22]" = 3, "(22,63.5]" = 4, "missing" = 0),
         wave.num = ifelse(wave == "first wave", 1 ,2))


# Join exposures and metabolites 
solar <- left_join(solar_eo, 
                   sol_pca_dat %>% mutate(id = as.character(id)), 
                   by = "id")

# PFAS
X.obs = solar[exposures_continuous] %>%
  mutate(across(everything(), ~scale(log2(.))))

# exclude outcome, leave only predictors:
Y = solar %>%
  dplyr::select(colnames(sol_pca_dat)[2]:ncol(solar)) %>% 
  scale(center = F, scale = F) 

dim(Y)

U = cbind.data.frame(age = as.numeric(solar$age), 
                     sex.num = as.numeric(solar$sex.num), 
                     bmi = as.numeric(solar$bmi), 
                     tanner = as.numeric(solar$tanner),
                     ses.num = as.numeric(solar$ses.num), 
                     wave.num = as.numeric(solar$wave.num))


P = ncol(X.obs)
LOD = c(0.01,0.05, 0.01, 0.01, 0.01, 0.43)
profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)
# exposure.Names = colnames(X.obs)

rm(list = setdiff(ls(), c("X.obs",
                          "Y",
                          "U",
                          "LOD",
                          "profiles", 
                          # "exposure.Names",
                          "ridge.BDL.model", 
                          "BHRMA.g")))

# Save data 
save.image(file = fs::path(dirname(here::here()),
                           "0_Data", 
                           "data_for_mixtures_pw_pca_analysis", 
                           "SOLAR_pw_pcs_mixtures_datasets_hyper_g_prior.Rdata"))

rm(list = ls())


# CHS -----------------------------------------------------------------
source(here::here("!directories.R"))
source(here::here("!set_exposure_outcome_vars.R"))
source(here::here("0_0_1_format_vars_funs.R"))
source(here::here("!load_data.R"))
source(here::here("1_mixtures_analysis",
                  "0_BHRMA.g_function.R"))

# Read in PCA Data 
pca_dat <- read_rds(fs::path(dir_results_mum_mixtures,
                             "SOL CHS PC Scores of sig pathways.RDS"))
# Get chs data only
chs_pca_dat <- pca_dat$chs

# Create numberic Vars
chs_eo <- exposure_outcome$chs %>% 
  mutate(id = as.character(id), 
         sex.num = ifelse(sex == "Female",1,0),
         ses.num = recode(ses, "[1,2]" = 1, "(2,4]" = 2, "(4,9]" = 3))


# Join exposures and metabolites 
chs <- left_join(chs_eo, 
                 chs_pca_dat %>% mutate(id = as.character(id)), 
                 by = "id")

# PFAS
X.obs = chs[exposures_continuous] %>%
  mutate(across(everything(), ~scale(log2(.))))

# exclude outcome, leave only predictors:
Y = chs %>%
  dplyr::select(colnames(chs_pca_dat)[2]:ncol(chs)) %>% 
  scale(center = F, scale = F) 

dim(Y)

U = cbind.data.frame(age = as.numeric(chs$age), 
                     sex.num = as.numeric(chs$sex.num), 
                     bmi = as.numeric(chs$bmi), 
                     ses.num = as.numeric(chs$ses.num))


P = ncol(X.obs)
LOD = c(0.01,0.05, 0.01, 0.01, 0.01, 0.43)
profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)


rm(list = setdiff(ls(), c("X.obs",
                          "Y",
                          "U",
                          "LOD",
                          "profiles", 
                          # "exposure.Names",
                          "ridge.BDL.model", 
                          "BHRMA.g")))

# Save data 
save.image(file = fs::path(dirname(here::here()),
                           "0_Data", 
                           "data_for_mixtures_pw_pca_analysis", 
                           "chs_pw_pcs_mixtures_datasets_hyper_g_prior.Rdata"))

rm(list = ls())
