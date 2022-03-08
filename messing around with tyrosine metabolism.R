
# t1 comes from 3_5
tyrosine <- t1 %>% 
  tidylog::filter(str_detect(pathway,
                             "Tyrosine"), 
                  mixture == "all pfas")


# length(unique(tyrosine$query_mass))


# There are 46 empirical compounds in this data set. I will average the pips across 
# each of these 46 compounds, and then plot the results

# Remove extra columns prior to reshaping data
tyrosine_wo_extracols <- tyrosine %>% 
  group_by(empirical_compound) %>%
  # summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  select(-contains("beta_ci"), -contains("lcl_beta"), -contains("ucl_beta"), 
         -contains("ucl_pip"), -contains("lcl_pip"), -contains("wald"), 
         -contains("sd"), -contains("q_value"), -contains("neg_log_p")) %>%
  rename_all(~ str_replace(., "p_value", "estimate_p") %>% 
               str_replace(., "mixture_effect", "mixture"))


# Pivot longer on all estimates, then wider ob
tyrosine_reshaped <- tyrosine_wo_extracols %>%
  pivot_longer(cols = estimate_beta_chs_mixture:estimate_pip_sol_pfos, 
               names_to = "variable_name", values_to = "estimate") %>% 
  separate(variable_name, sep = "_",
           into = c("empty", "term", "cohort", "pfas")) %>% 
  pivot_wider(values_from = estimate, names_from = term) %>%
  select(-empty) %>% 
  ungroup()


# create seperate DFs for mixtures vs. pfas effects
mixture_only <- tyrosine_reshaped %>% 
  filter(pfas == "mixture") %>% 
  select(mixture, empirical_compound, matched_form, name, matched_compound, cohort, beta, p) %>% 
  rename(beta_mixture = beta, p_mixture = p)

# Merge with pfas only
tyrosine_final <- tyrosine_reshaped %>% filter(pfas != "mixture") %>% 
  tidylog::full_join(mixture_only)


# Summarize: 
tyrosine_final_2 <- tyrosine_final %>% 
  tidylog::filter(p_mixture < 0.05) %>%
  group_by(empirical_compound, pfas, cohort) %>% 
  summarise(met_name = met_name[1],
            met_name_lc = tolower(met_name),
            pct_pip80 = sum(pip>0.8)/length(pip), 
            meanpip = mean(pip),
            max_pip = max(pip), 
            p_mixture = p_mixture[1], 
            mean_beta = mean(beta), 
            sign_beta = sign(mean_beta), 
            pip_times_sign = meanpip*sign_beta) %>% 
  ungroup() %>% 
  mutate(pfas = rename_pfas(pfas, arrange_by_class = TRUE), 
         met_name_first_only = str_split_fixed(met_name, "; ", n = 2)[,1])


# group compounds by whether or not they have a single pip > 0.5
ec_by_maxpip <- tyrosine_final_2 %>% 
  group_by(met_name) %>% 
  summarise(max_pip_for_cpd = max(meanpip, na.rm = TRUE), 
            max_pip_cpd_over_80 = max_pip_for_cpd > .8)

tyrosine_final_summarized <- tidylog::left_join(tyrosine_final_2, 
                                       ec_by_maxpip)





# Heatmap of mean PIPs ----------------------
ggplot(tyrosine_final_summarized,
                          aes(x = pfas,
                              y = met_name_first_only,
                              # fill = pip_times_sign)) + 
                              fill = meanpip)) + 
    geom_tile() +
    theme(#axis.text.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x =  element_text(angle = 90,
                                      hjust = 1, 
                                      vjust = 0.5)) + 
    scale_fill_gradient2(midpoint = 0,
                         low = "blue",
                         high = "red")+
    facet_grid(max_pip_cpd_over_80~cohort, space = "free_y", scales = "free_y")



tyrosine_final_summarized %>% 
  group_by(pfas, cohort) %>% 
  summarise(pips_above_80 = npct(meanpip>0.8, TRUE)) %>% 
  ungroup() %>%
  pivot_wider(names_from = cohort, values_from = pips_above_80)

#
#























# Linegraph of mean pips --------------------
ggplot(tyrosine_final_summarized, aes(x = pfas, 
                                      y = meanpip, 
                                      # shape = cohort, 
                                      color = cohort, 
                                      group = empirical_compound)) + 
  
  # geom_boxplot(position = position_dodge(width = .5)) +
  geom_point() +
  geom_line() +
  facet_wrap(~cohort, nrow = 2) +
  geom_hline(yintercept = 1/6) +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.text.x =  element_text(angle = 90), 
        strip.text = element_blank()) 





