### BHRMA
source(file=paste0("BHRMA.g.R")) 
exp_cov_dat = exposure_outcome$solar
metab_dat = ftdata$solar$c18neg
name_of_exposure = exposures_continuous

data <- left_join(exp_cov_dat %>% mutate(id = as.character(id)), 
                  metab_dat %>% mutate(id = as.character(id)), 
                  by = "id") %>%
        mutate(sex.num = ifelse(sex == "Female",1,0),
               ses.num = recode(ses, "[3,11]" = 1, "(11,15.5]" = 2, "(15.5,22]" = 3, "(22,63.5]" = 4, "missing" = 0),
               wave.num = ifelse(wave == "first wave", 1 ,2))

X.obs = data[name_of_exposure]
Y = data %>%
  dplyr::select(colnames(metab_dat)[2]:ncol(data)) %>% # exclude outcome, leave only predictors
  scale(center = F, scale = F) 
  

U = cbind.data.frame(as.numeric(data$age), 
                     as.numeric(data$sex.num), 
                     as.numeric(data$bmi), 
                     as.numeric(data$tanner),
                     as.numeric(data$ses.num), 
                     as.numeric(data$wave.num)
                     )
P = ncol(X.obs)
LOD = c(0.01, 0.01, 0.05, 0.01, 0.01, 0.01, 0.43, 0.01)
profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)
exposure.Names = colnames(X.obs)

#X=X.obs; Y=Y[,1]; U=U; LOD=LOD; profiles=profiles

ptm <- proc.time()
fit = BHRMA.g(X=X.obs, Y=Y[,1], U=U, LOD=LOD, profiles=profiles)
fit
(proc.time() - ptm)/60

### Bayesian ridge with non-informative prior
# source(file=paste0("BHRMA_non-informative_prior.R"))
# ptm <- proc.time()
# set.seed(2021)
# fit = BHRMA.noninfo(X=X.obs, Y=Y[,1], U=U, LOD=LOD, profiles=profiles)
# fit
# (proc.time() - ptm)/60
