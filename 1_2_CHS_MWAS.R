# script for MWAS analysis (regression)
# C8 and hilic metabolites separated
# Jesse Goodrich 040821


# C18 --------------------------------------
c18_mwas_reg_chs = function(exposure, met_name){
  covar = 'chs_exposure_outcome$age + 
  chs_exposure_outcome$sex + 
  chs_exposure_outcome$bmi + 
  chs_exposure_outcome$ses'
  formula = paste0('chs_c18$',
                   '`',met_name,'`',
                   '~',
                   'chs_exposure_outcome$',
                   exposure,'+', covar)
  lmfit = lm(formula = formula) %>% 
    tidy(., conf.int = TRUE)
  
  #p-value for exposure in regression
  p.value =  lmfit$p.value[2]
  #t.score for exposure in regression
  t.score = lmfit$statistic[2]
  #beta & confint for exposure in regression
  beta = lmfit$estimate[2]
  conf.low = lmfit$conf.low[2]
  conf.high = lmfit$conf.high[2]
  
  m.z = as.numeric(str_split(met_name,'_')[[1]][1])
  r.t = as.numeric(str_split(met_name,'_')[[1]][2])
  mode = 'negative'
  return(c(m.z, r.t, p.value, t.score, beta, conf.low, conf.high, mode)) 
}

input_data = cbind.data.frame(rep(exposures, 
                                  each = ncol(chs_c18)-1 ), 
                              rep(colnames(chs_c18)[2:ncol(chs_c18)], 
                                  length(exposures)))

names(input_data) = c('exposure','met_name')

#need to use plyr function
#this will take some time
c18_mwas_results = plyr::mdply(input_data,.fun = c18_mwas_reg_chs,.progress = 'text')
names(c18_mwas_results) = c('exposure','metabolite',"m.z", "r.t", 
                            "p.value", "t.score", "beta", 
                            "conf.low", "conf.high", "mode")

#version 1: new LC-MS Data
saveRDS(c18_mwas_results, file = here::here('Temporary results', exposure_type,
                                            '1_1_chs_c18_mwas.rds'))

##HILIC mwas ---------------------------------------------
hilic_mwas_reg = function(exposure, met_name){
  covar = 'chs_exposure_outcome$age + 
  chs_exposure_outcome$sex + 
  chs_exposure_outcome$bmi + 
  chs_exposure_outcome$ses'
  formula = paste0('chs_hilic$',
                   '`',met_name,'`',
                   '~',
                   'chs_exposure_outcome$',
                   exposure,'+', covar)
  lmfit = lm(formula = formula) %>% 
    tidy(., conf.int = TRUE)
  
  #p-value for exposure in regression
  p.value =  lmfit$p.value[2]
  #t.score for exposure in regression
  t.score = lmfit$statistic[2]
  #beta & confint for exposure in regression
  beta = lmfit$estimate[2]
  conf.low = lmfit$conf.low[2]
  conf.high = lmfit$conf.high[2]
  
  m.z = as.numeric(str_split(met_name,'_')[[1]][1])
  r.t = as.numeric(str_split(met_name,'_')[[1]][2])
  mode = 'positive'
  return(c(m.z, r.t, p.value, t.score, beta, conf.low, conf.high, mode)) 
}

#this is slightly different from c18
input_data = cbind.data.frame(rep(exposures,
                                  each = ncol(chs_hilic)-1), 
                              rep(colnames(chs_hilic)[2:ncol(chs_hilic)], 
                                  length(exposures))  )
names(input_data) = c('exposure','met_name')

#need to use plyr function
#this will take some time
hilic_mwas_results = plyr::mdply(input_data,
                                 .fun = hilic_mwas_reg,
                                 .progress = 'text')
names(hilic_mwas_results) = c('exposure','metabolite',"m.z", "r.t", 
                              "p.value", "t.score", "beta", 
                              "conf.low", "conf.high", "mode")

# Save HILIC MWAS 
saveRDS(hilic_mwas_results, 
        file = here::here('Temporary results', exposure_type, 
                          '1_1_chs_hilic_mwas.rds'))



# Save into separate folders---------------------------------------------
exposure_merge_n_save = function(exposure_name){
  c18_data = 
    c18_mwas_results %>% 
    dplyr::filter(exposure==exposure_name)  %>% 
    ungroup() %>%
    dplyr::select(m.z,r.t ,beta, t.score, p.value,beta, conf.low, conf.high, mode)
  
  hilic_data = 
    hilic_mwas_results %>% 
    dplyr::filter(exposure==exposure_name)  %>% 
    ungroup() %>%
    dplyr::select( m.z,r.t,beta, t.score, p.value,beta, conf.low, conf.high, mode)
  
  merged_data = dplyr::union(c18_data,hilic_data) 
  
  #save csv
  write_csv(merged_data,
            file = fs::path(here::here('Temporary results'), exposure_type, 
                            paste0("CHS_", exposure_type,'_MWAS_c18_hilic.csv')))
}

#run the merge and save for all exposure
plyr::l_ply(exposures, .fun = exposure_merge_n_save, .progress = 'text')
