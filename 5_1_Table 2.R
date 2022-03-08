# Summarize SOLAR and CHS PFAS Data -----------------
sol_pfas <- exposure_outcome$solar %>% 
  select(all_of(exposures_continuous)) %>% 
  mutate(cohort = "SOLAR")

chs_pfas <- exposure_outcome$chs %>% 
  select(all_of(exposures_continuous)) %>% 
  mutate(cohort = "CHS")

# Bind rows and pivot longer
pfas_l <- bind_rows(sol_pfas, chs_pfas) %>% 
  pivot_longer(cols = all_of(exposures_continuous), 
               values_to = "concentration", 
               names_to = "pfas")


# By cohort
pfas_summary <- pfas_l %>% 
  group_by(cohort, pfas) %>%
  summarise(
    geometric_mean = fungm(concentration), 
    percentile_50 = qntle_fxn(concentration, .50), 
    percentile_75 = qntle_fxn(concentration, .75), 
    percentile_90 = qntle_fxn(concentration, .9), 
    pct_below_lod = npct(concentration, 
                         level_of_interest = 0, 
                         n.digits = 2)) %>% 
  ungroup() %>% 
  mutate(cohort = fct_relevel(cohort, "SOLAR")) %>% 
  arrange(cohort,pfas)


# Calculate p-values
(p_values <- pfas_l %>% 
    mutate(concentration = log(concentration)) %>%
    group_by(pfas) %>% 
    nest() %>% 
    mutate(ttest = map(data, 
                       ~t.test(concentration ~ cohort,
                               data = .x, paired = FALSE) %>%
                         tidy())) %>%
    select(pfas, ttest) %>% 
    unnest(cols = c(ttest)) %>% 
    select(pfas, p.value) %>% 
    mutate(p.value = formatC(p.value, digits = 2, format = 'g')))


## Pivot Summary Data Wider on Cohort
pfas_summary_w <- pfas_summary %>% 
  select(cohort, pfas, geometric_mean) %>%
  pivot_wider(id_cols = pfas,
              names_from = cohort, 
              values_from = geometric_mean) %>% 
  select(pfas, contains("SOLAR"), contains("CHS")) 

# join p values 
pfas_summary_final <- pfas_summary_w %>% 
  tidylog::left_join(p_values) %>% 
  mutate(pfas = rename_pfas(pfas) %>% 
           fct_relevel("PFOS", "PFHxS", "PFHpS", "PFOA", "PFNA", "PFDA")) %>% 
  rename(PFAS = pfas) %>% 
  arrange(PFAS)


# Save Data
write_csv(pfas_summary_final, 
          here::here(dir_reports,
                     "Table 2.csv"))
