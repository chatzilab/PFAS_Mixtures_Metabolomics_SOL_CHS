# plot PIPs for tyrosine metabolism

# Read in Data
annotated_sig_ecs_ee <- read_rds(
  here::here(dir_results_mixtures,
             "Sig annotated metabolite effect estimates sol chs.RDS") ) %>% 
  select(-contains("blank")) %>% 
  tidylog::filter(str_detect(pathway, "fatty acid"))


# Subset SOLAR and CHS
sol <- annotated_sig_ecs_ee %>% 
  select(met_name:mass_diff, contains("sol")) %>% 
  tidylog::filter(ucl_beta_sol_mixture <= 0 | 
                    lcl_beta_sol_mixture >= 0) %>% 
  group_by(met_name) %>% 
  filter(estimate_beta_sol_mixture == max(estimate_beta_sol_mixture)) %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                       estimate_beta_sol_mixture), 
         cohort = "SOLAR") %>% 
  rename_all(~str_remove(., "_sol"))


chs <- annotated_sig_ecs_ee %>% 
  select(met_name:mass_diff, contains("chs")) %>% 
  tidylog::filter(ucl_beta_chs_mixture <= 0 | 
                    lcl_beta_chs_mixture >= 0) %>% 
  group_by(met_name) %>% 
  tidylog::filter(estimate_beta_chs_mixture == max(estimate_beta_chs_mixture)) %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                       estimate_beta_chs_mixture), 
         cohort = "CHS") %>% 
  rename_all(~str_remove(., "_chs"))


# Combine
sig_ecs_l <- bind_rows(sol, chs) %>% 
  select(cohort, everything())


sig_ecs_l_2 <- sig_ecs_l %>% 
  select(-contains("beta_ci_"), -contains("neg_log_p"),
         -contains("sd_"), -contains("wald_")) %>% 
  mutate(across(estimate_beta_mixture:q_value_pfos, as.numeric)) %>%
  pivot_longer(cols = estimate_beta_mixture:q_value_pfos, 
               names_to = "effect_name")

# seperate effect and term
sig_ecs_l_3 <- sig_ecs_l_2 %>% 
  mutate(effect_1 = str_split_fixed(effect_name, "_", 3)[,1], 
         effect_2 = str_split_fixed(effect_name, "_", 3)[,2],
         exposure = str_split_fixed(effect_name, "_", 3)[,3], 
         effect = str_c(effect_1, effect_2, sep = "_")) %>% 
  select(-effect_1, -effect_2, effect_name)

# Pivot wider
fa_mets <- sig_ecs_l_3 %>% 
  pivot_wider(id_cols = c(cohort:mass_diff, exposure), 
              names_from = effect, 
              values_from = value)  %>% 
  select(cohort, exposure, everything()) %>% 
  mutate(estimate_b = estimate_beta*estimate_pip, 
         sig_beta = if_else(lcl_beta >= 0 |
                              ucl_beta <= 0, estimate_beta, 0), 
         pfas = rename_pfas(exposure,arrange_by_class = TRUE))


sol_data <- fa_mets %>% 
  filter(cohort == "SOLAR", 
         exposure != "mixture")

chs_data <-  fa_mets %>% 
  filter(cohort == "CHS", 
         exposure != "mixture")


# Plot -----------------------------------------

(sol <- ggplot(sol_data, 
                    aes(x = pfas,
                        y = met_name,
                        fill = estimate_pip)) + 
   geom_tile() +
   theme(axis.title.y = element_blank(), 
         axis.title.x = element_blank(), 
         strip.text.x = element_blank(), 
         axis.text.x = element_blank(), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         strip.text.y = element_text(angle = 0, hjust = 0)) +
   binned_scale(aesthetics = "fill",
                scale_name = "stepsn",
                palette = function(x) c("white", "grey90", "grey70","grey60","grey50","grey40"),
                breaks = signif(c(  1/6, 2/6, 3/6, 4/6, 5/6), digits = 2),
                labels = c("0.17",  "",  "0.5",  "",  "0.83"),
                limits = c(0, 1),
                show.limits = FALSE,
                guide = "colorsteps", 
                name = "Posterior inclusion probability"))
  

(chs <- ggplot(chs_data, 
               aes(x = pfas,
                   y = met_name,
                   fill = estimate_pip)) + 
    geom_tile() + 
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          strip.text.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom",
          legend.justification = c("center","bottom"),
          strip.text.y = element_text(angle = 0, hjust = 0), 
          legend.key.size = unit(.6, "cm")) +
    binned_scale(aesthetics = "fill",
                 scale_name = "stepsn",
                 palette = function(x) c("white", "grey90", "grey70","grey60","grey50","grey40"),
                 breaks = signif(c(  1/6, 2/6, 3/6, 4/6, 5/6), digits = 2),
                 labels = c("0.17",  "",  "0.5",  "",  "0.83"),
                 limits = c(0, 1),
                 show.limits = FALSE,
                 guide = "colorsteps", 
                 name = "Posterior inclusion probability"))






# Combine Plots -----------------------
(fa_met_pips <- plot_grid(NULL, sol, 
                      NULL, chs, 
                      ncol = 1, align = "v", 
                      rel_heights = c(.05, 1,.05, 1.2),
                      labels = c("A) SOLAR","",
                                 "B) CHS", "")))


# ggsave(fa_met_pips, 
#        filename = fs::path(dir_reports, 
#                            "Figure S2 PFAS and fa metabolite pips.jpg"), 
#        width = 6, height = 5)



# Summary of PIPs ------------------------------

quantile(sol_data$estimate_pip, 0.9, na.rm = TRUE)



sol_data %>% 
  group_by(pfas) %>% 
  summarise(n_pip_above_50 = sum(estimate_pip>.88)/length(estimate_pip))



summary(chs_data$estimate_pip)
quantile(chs_data$estimate_pip, 0.9, na.rm = TRUE)


chs_data %>% 
  group_by(pfas) %>% 
  summarise(n_pip_above_50 = sum(estimate_pip>.63)/length(estimate_pip))
