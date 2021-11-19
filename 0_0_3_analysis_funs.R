
# MWAS FUNCTION -----------------------------
mwas <- function(metab_dat, exp_cov_dat, name_of_exposure, cohort){
  # merge data to ensure matching ids
  data <- left_join(exp_cov_dat %>% mutate(id = as.character(id)), 
                    metab_dat %>% mutate(id = as.character(id)), 
                    by = "id")
  # only run analysis if exposure var has more than 1 value
  if(length(unique(data[[name_of_exposure]])) >1){
    # Run MWAS
    if(cohort == "solar"){ 
      MWAS_output_lm <- data %>%
        dplyr::select(colnames(metab_dat)[2]:ncol(data)) %>%  # exclude outcome, leave only predictors
        as_tibble() %>%
        map(~lm(.x ~ data[[name_of_exposure]] +
                  data$age + data$sex + data$bmi + 
                  data$tanner + data$ses + data$wave, 
                data = data)) %>%
        map(~tidy(.x, conf.int = TRUE)) %>% 
        bind_rows(.id = "name")
    }
    
    if(cohort == "chs"){ 
      MWAS_output_lm <- data %>%
        dplyr::select(colnames(metab_dat)[2]:ncol(data)) %>%  # exclude outcome, leave only predictors
        as_tibble() %>%
        map(~lm(.x ~ data[[name_of_exposure]] +
                  data$age + data$sex + data$ses, 
                data = data)) %>%
        map(~tidy(.x, conf.int = TRUE)) %>% 
        bind_rows(.id = "name")
    }
    
    
    # Add new variables, filter only the exposure var
    MWAS_output <- MWAS_output_lm %>% 
      mutate(exposure = name_of_exposure, 
             term = if_else(str_detect(term, "name_of_exposure"), 
                            name_of_exposure, 
                            term)) %>% 
      filter(term == name_of_exposure) %>%
      janitor::clean_names()
    
    nfeat <- length(unique(MWAS_output$name))
    
    # Calculate FDR Adjustments
    MWAS_output <- MWAS_output %>% 
      group_by(term) %>%
      mutate(q_value = p.adjust(p_value, 
                                method = "fdr", 
                                n = nfeat), # Calculate FDR Adjustment  
             significance = ifelse(p_value < 0.05, "P<0.05", "Not Sig"), 
             significancefdr = ifelse(q_value < 0.05, "P<0.05", "Not Sig"), 
             names = if_else(p_value < 0.05, name, ""), 
             namesfdr = if_else(q_value < 0.05, name, "")) %>% 
      ungroup()
    
    return(MWAS_output)
  }
}





# Run mode specific MWAS for all exposures
# cohort = "solar"
# lcms_mode = "c18neg"
# metab_dat = temp_met
# exp_cov_dat = exposure_outcome
# analysis_exposures = exposures_for_analysis
# i = analysis_exposures[1]
# rm(cohort,lcms_mode, metab_dat, exp_cov_dat)

mwas_all_exposures <- function(cohort, lcms_mode, metab_dat, exp_cov_dat,
                          analysis_exposures){
  # Initialize list
  list_mwas_results <- vector("list", length(analysis_exposures))
  names(list_mwas_results) <- analysis_exposures
  
  # Run MWAS 
  for(i in analysis_exposures){
    print(paste(i))
    sol_mwas <- mwas(metab_dat = metab_dat[[cohort]][[lcms_mode]], 
                     exp_cov_dat = exp_cov_dat[[cohort]], 
                     name_of_exposure = i, 
                     cohort = cohort) %>% 
      list()
    
    list_mwas_results[i] = sol_mwas
    # print(paste("End", i, Sys.time()))
  }
  list_mwas_results <- bind_rows(list_mwas_results, .id = "exposure") %>% 
    select(exposure, everything(), -term)
  
  return(list_mwas_results)
}

