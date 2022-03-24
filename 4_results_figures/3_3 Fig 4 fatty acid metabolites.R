# Make coefficient plots for PFAS and FFA metabolites

# Read data
annotated_sig_ecs_ee <- read_rds(
  here::here(dir_results_mixtures,
             "Sig annotated metabolite effect estimates sol chs.RDS") ) %>% 
  select(-contains("blank")) %>% 
  tidylog::filter(str_detect(pathway, "fatty acid"))

# write_csv(annotated_sig_ecs_ee, 
#           here::here("De novo FA biosynthesis cpds.csv"))

# Change metabolite names to those associated with Tyrosine metabolism --------------
annotated_sig_ecs_ee <- annotated_sig_ecs_ee %>% 
  tidylog::mutate(met_name_fa_pw = case_when(
    met_name == "Oleic acid; Elaidic acid" ~ "Elaidic acid", 
    TRUE ~ met_name), 
    borderline_sig_sol = if_else(lcl_beta_sol_mixture == 0 | 
                                   ucl_beta_sol_mixture == 0, 
                                 "borderline", 
                                 "sig"), 
    borderline_sig_chs = if_else(lcl_beta_chs_mixture == 0 | 
                                   ucl_beta_chs_mixture == 0, 
                                 "borderline", 
                                 "sig"))

length(unique(annotated_sig_ecs_ee$met_name))
length(unique(annotated_sig_ecs_ee$empirical_compound))


# Make Plots ----------------------------------------------------------
# Subset SOLAR and CHS
sol <- annotated_sig_ecs_ee %>% 
  select(met_name_fa_pw, pathway:mass_diff, contains("sol")) %>% 
  select(met_name_fa_pw, pathway:mass_diff, 
         contains("estimate"),
         contains("beta"),
         everything()) %>%
  tidylog::filter(ucl_beta_sol_mixture <= 0 | 
                    lcl_beta_sol_mixture >= 0) %>% 
  group_by(met_name_fa_pw) %>% 
  filter(estimate_beta_sol_mixture == max(estimate_beta_sol_mixture)) %>% 
  ungroup() %>% 
  mutate(met_name_fa_pw = fct_reorder(met_name_fa_pw, 
                                      estimate_beta_sol_mixture))



chs <- annotated_sig_ecs_ee %>% 
  select(met_name_fa_pw, pathway:mass_diff, contains("chs")) %>% 
  select(met_name_fa_pw, pathway:mass_diff, 
         contains("estimate"),
         contains("beta"),
         everything()) %>%
  tidylog::filter(ucl_beta_chs_mixture <= 0 | 
                    lcl_beta_chs_mixture >= 0) %>% 
  group_by(met_name_fa_pw) %>% 
  tidylog::filter(estimate_beta_chs_mixture == max(estimate_beta_chs_mixture)) %>% 
  ungroup() %>% 
  mutate(met_name_fa_pw = fct_reorder(met_name_fa_pw, 
                                      estimate_beta_chs_mixture))

# Get data on overlap
sol_fa <- sol %>% 
  mutate(sig_both = met_name_fa_pw %in% chs$met_name_fa_pw, 
         met_name_fa_pw_with_overlap = if_else(
           sig_both, 
           str_c(met_name_fa_pw, "*"), 
           as.character(met_name_fa_pw)) %>% 
           fct_reorder(estimate_beta_sol_mixture))


chs_fa <- chs %>% 
  mutate(sig_both = met_name_fa_pw %in% sol$met_name_fa_pw, 
         met_name_fa_pw_with_overlap = if_else(sig_both, 
                                               str_c(met_name_fa_pw, "*"), 
                                               as.character(met_name_fa_pw)) %>% 
           fct_reorder(estimate_beta_chs_mixture))

# Plot solar------------------------------
(sol_met_plot <- ggplot(sol_fa, aes(x = met_name_fa_pw_with_overlap, 
                                    y = estimate_beta_sol_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_sol_mixture , 
                     ymax = ucl_beta_sol_mixture, 
                     linetype = borderline_sig_sol), 
                 width = 0) +
   geom_point() + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   coord_flip() + 
   scale_y_continuous(limits = c(-1.75, 3)) +
   scale_linetype_manual(values = c(5, 1)) +
   # ylab("PFAS Exposure Effect Estimate ψ (95% BCI)") +
   theme(axis.title.y = element_blank(), 
         axis.title.x = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         strip.text.y = element_text(angle = 0, hjust = 0)))



# Plot ------------------------------
(chs_metplot <- ggplot(chs_fa, aes(x = met_name_fa_pw_with_overlap, 
                                   y = estimate_beta_chs_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_chs_mixture , 
                     ymax = ucl_beta_chs_mixture, 
                     linetype = borderline_sig_chs), 
                 width = 0) + 
   geom_point() +  #aes(color = tyr_subpath)
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   coord_flip() + 
   scale_y_continuous(limits = c(-1.75, 3)) +
   scale_linetype_manual(values = c(5, 1)) +
   ylab("PFAS Mixture Effect ψ (95% CI)") +
   theme(axis.title.y = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey97"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         strip.text.y = element_text(angle = 0, hjust = 0)))


# Combine Plots
(figure_4 <- plot_grid(NULL, sol_met_plot, 
                       NULL, chs_metplot, 
                       ncol = 1, align = "v", 
                       rel_heights = c(.07, 1,.07, .9),
                       labels = c("A) SOLAR","",
                                  "B) CHS", "")))

# Save 
ggsave(figure_4, 
       filename = fs::path(dir_reports, 
                           "Figure 4 associations of PFAS and FA metabolites.jpg"), 
       width = 5, height = 3.5)
