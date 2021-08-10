## 2_7_SOLAR CHS Create Annotated Feature Table

solar_single_matches_w_mwas <- readRDS(here::here("Temporary results", 
                                                  "solar_mz_cpd_pathway_key_w_mwas.rds"))


# Select top 50 features
tentative_annotations <- solar_single_matches_w_mwas %>% 
  group_by(metabolite) %>% 
  slice_head() %>% 
  ungroup()

tentative_annotations <- tentative_annotations %>% 
  select(metabolite:main_sp_pthw_categorized, mode, super_pathway, n_cpds_per_mz, n_mz_per_cpd)




# Get positive/neg mode 
solar_c18_metab <- tentative_annotations %>% filter(mode == "negative")
solar_hilic_metab <- solar_single_matches_w_mwas %>%  filter(mode == "positive")

solar_c18_temp <- solar_c18 %>% select(id, all_of(solar_c18_metab$metabolite)) 
solar_hilic_temp <- solar_hilic %>% select(id, all_of(solar_hilic_metab$metabolite))

solar_ft_anno <- left_join(solar_c18_temp, solar_hilic_temp, by = "id") 


##
# C18
#create a adjustment variable
solar_adjustvar_temp_3 = solar_exposure_outcome %>% 
  filter(!is.na(dde)) %>%
  arrange(id) %>%
  dplyr::select(id,age,sex,bmi,tanner,ses,wave)

solar_adjustvar_temp_2 = fastDummies::dummy_cols(solar_adjustvar_temp_3,
                                          select_columns = c("sex", "tanner","ses", "wave")) %>% 
  select(-c(sex, tanner, ses, wave))


## C18 variable
solar_ft_temp <- solar_ft_anno %>% 
  dplyr::filter(id %in% solar_adjustvar_temp_2$id) %>% 
  arrange(id) %>%  
  # select(-id) %>% 
  hciR::as_matrix() 

solar_ft_temp2 <- apply(solar_ft_temp, 2, function(x){ifelse(is.na(x), min(x, na.rm = TRUE), x)})

table(is.na(solar_ft_temp))
# Covariate matrix
solar_adjustVar = solar_adjustvar_temp_2 %>% 
  select(-id) %>% as.data.frame() %>%
  as.matrix()

#get c18 residual after regression
solar_reg = lm(solar_ft_temp2 ~ solar_adjustVar)
solar_resid <- solar_reg$residuals %>% 
  as.data.frame()

rownames(solar_resid) = rownames(solar_ft_temp)
colnames(c18_residual) = colnames(c18_2)[-1]