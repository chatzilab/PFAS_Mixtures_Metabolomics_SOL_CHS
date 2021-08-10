

xxx <- single_matches_w_mwas_final %>% 
  filter(exposure2 == "HCB")

write_csv(xxx, here::here("hcb_metabolites.csv"))

length(unique(xxx$cpd))

xxx %>% 
  group_by(mz_solar) %>% 
  summarise(n_cpds_per_mz = length(cpd)) %>% 
  ggplot(aes(n)) + geom_histogram()
