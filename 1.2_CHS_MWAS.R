# script for MWAS analysis (regression)
# C8 and hilic metabolites separated
# Jesse Goodrich 040821

exposures = c("hexachlorobenzene_impute", 
              "dde_impute",
              "ocs",
              "pbde_154_impute",
              "pbde_47_impute", 
              "pbde_100_ngml_detect",
              "pbde_153_ngml_detect",    
              "pbde_85_ngml_detect",          
              "pcb_118_ngml_detect",   
              "pcb_138_ngml_detect",          
              "pcb_153_ngml_detect",     
              "pcb_180_ngml_detect", 
              "pcb_num_detect", 
              "pbde_num_detect")

# C18 --------------------------------------
C18_MWAS_reg = function(exposure, met_name){
  covar = 'chs_exposure_outcome$age + 
  chs_exposure_outcome$sex + 
  chs_exposure_outcome$bmi + 
  chs_exposure_outcome$ses + 
  chs_exposure_outcome$wave'
  formula = paste0('chs_c18$',
                   '`',met_name,'`',
                   '~',
                   'chs_exposure_outcome$',
                   exposure,'+', covar)
  lmfit = lm(formula = formula)
  
  #p-value for exposure in regression
  p.value =  coef(summary(lmfit))[2,4]
  #t.score for exposure in regression
  t.score = coef(summary(lmfit))[2,3]
  
  m.z = as.numeric(str_split(met_name,'_')[[1]][1])
  r.t = as.numeric(str_split(met_name,'_')[[1]][2])
  mode = 'negative'
  return(c(m.z, r.t, p.value, t.score, mode)) 
}

input_data = cbind.data.frame(rep(exposures, 
                                  each = ncol(chs_c18)-1 ), 
                              rep(colnames(chs_c18)[2:ncol(chs_c18)], 
                                  length(exposures)))
names(input_data) = c('exposure','met_name')

#need to use plyr function
#this will take some time
c18_mwas_results = plyr::mdply(input_data,.fun = C18_MWAS_reg,.progress = 'text')
names(c18_mwas_results) = c('exposure','metabolite','m.z','r.t', 'p.value', 't.score', 'mode' )

#version 1: new LC-MS Data
saveRDS(c18_mwas_results, file = here::here('Temporary results','1.1.0_chs_c18_mwas.rds'))



##HILIC mwas ---------------------------------------------
HILIC_MWAS_reg = function(exposure, met_name){
  covar = 'chs_exposure_outcome$age + 
  chs_exposure_outcome$sex + 
  chs_exposure_outcome$bmi + 
  chs_exposure_outcome$ses + 
  chs_exposure_outcome$wave'
  formula = paste0('chs_hilic$',
                   '`',met_name,'`',
                   '~',
                   'chs_exposure_outcome$',
                   exposure,'+', covar)
  lmfit = lm(formula = formula)
  
  #p-value for exposure in regression
  p.value =  coef(summary(lmfit))[2,4]
  #t.score for exposure in regression
  t.score = coef(summary(lmfit))[2,3]
  
  m.z = as.numeric(str_split(met_name,'_')[[1]][1])
  r.t = as.numeric(str_split(met_name,'_')[[1]][2])
  mode = 'positive'
  return(c(m.z,r.t, p.value, t.score,mode)) 
}


#this is slightly different from c18
input_data = cbind.data.frame(rep(exposures,
                                  each = ncol(chs_hilic)-1), 
                              rep(colnames(chs_hilic)[2:ncol(chs_hilic)], 
                                  length(exposures))  )
names(input_data) = c('exposure','met_name')

#need to use plyr function
#this will take some time
HILIC_mwas_results = plyr::mdply(input_data,.fun = HILIC_MWAS_reg,.progress = 'text')
names(HILIC_mwas_results) = c('exposure','metabolite','m.z','r.t', 'p.value', 't.score','mode' )

# Save HILIC MWAS 
saveRDS(HILIC_mwas_results, 
        file = here::here('Temporary results','1.1.0_chs_hilic_mwas.rds'))