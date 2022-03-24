# Plot Mummichog Pathway Results
library(colorspace)
library(janitor)
source(here::here("0_project_setup", "!directories.R"))

# Read in data
mum_pw_long <- read_rds(
  fs::path(dir_results_mum_mixtures, 
           "Mixture effect hyper_g", 
           "SOL CHS PFAS Mummichog long sig PW.RDS")) %>% 
  tidylog::filter(mixture == "all_pfas")


# Final Dataset for Bubbleplot: format "name" to remove _sol/_chs/_meta, then pivot wider
mum_pw_long2 <- mum_pw_long %>% 
  mutate(cohort = case_when(str_detect(cohort, "solar") ~ "SOLAR", 
                            str_detect(cohort, "chs") ~ "CHS", 
                            str_detect(cohort, "meta") ~ "Meta-analysis") %>% 
           as.factor() %>% 
           fct_relevel("SOLAR", "CHS", "Meta-analysis"), 
         hits_sig = if_else(cohort == "Meta-analysis", 
                            mean(hits_sig, na.rm = TRUE),
                            hits_sig), 
         shape = if_else(cohort == "Meta-analysis", 
                         "Meta p-value",
                         "Cohort Specific p-value"))

# Get vector to order pathways: 
pw_order <- mum_pw_long2 %>% 
  filter(cohort == "Meta-analysis",
         mixture == "all_pfas") %>% 
  rename(meta_neg_logp = neg_logp) %>% 
  select(path, meta_neg_logp)

# Join pw order and pw_data
mum_pw_long3 <- mum_pw_long2 %>% 
  tidylog::left_join(pw_order)

# Get number of significant pathways 
mum_summary <- mum_pw_long3 %>% 
  group_by(cohort, mixture) %>% 
  summarise(n_sig = sum(combined_pvals < 0.05, na.rm = TRUE), 
            length = length(combined_pvals))

# Adjust super pathways based on final results 
mum_pw_long4 <- mum_pw_long3 %>% 
  mutate(
    super_pathway_new = case_when(str_detect(path, "Alkaloid") ~ "Other", 
                                  str_detect(path, "Nitrogen") ~ "Other", 
                                  str_detect(path, "cytochrome P450") ~ "Other",
                                  str_detect(super_pathway, "Aromatic Amino Acids") ~ 
                                    "Aromatic Amino Acid Metabolism",
                                  TRUE ~ super_pathway), 
    path_new = path_2 %>% str_replace("metabolites formation", 
                                         "metabolite formation") %>% 
      str_replace("Putative anti-Inflammatory", 
                  "Anti-inflammatory") %>% 
      str_replace("met.", 
                  "metabolism"))

table(mum_pw_long4$super_pathway_new)


# Filter only sig. meta p-vals, order by meta p
mum_pw_long_final <- mum_pw_long4 %>%
  tidylog::filter(str_detect(sig_cohort, "Sig. ") | 
                    (sig_meta == "Sig.")) %>%
  arrange(-meta_neg_logp) %>%
  group_by(mixture) %>% 
  nest() %>% 
  mutate(data = map(data, 
                    ~mutate(.x, 
                            path_new = fct_reorder(path_new, 
                                                 meta_neg_logp,.desc = TRUE))))


# Sanity check:should have the same number of pathways in each cohort/meta
table(mum_pw_long$mixture, 
      mum_pw_long$cohort)

table(is.na(mum_pw_long_final$data[[1]]$hits_sig),
      mum_pw_long_final$data[[1]]$cohort)


# Bubble plot ----------------------------------------------------
bubbleplot_fxn <- function(data, mixture_name){
  data2 <- data[which(data$mixture == mixture_name), 2][[1]][[1]]
  
  ggplot(data2,
         aes(x = neg_logp,
             y = fct_rev(path_new),
             # color = super_pathway,
             size = hits_sig, 
             shape = shape)) +
    geom_point(alpha = .75) +
    geom_vline(xintercept = -log10(0.05),
               color = "grey50", 
               linetype = 2) + 
    facet_grid(super_pathway_new~cohort, 
               scales = "free_y", 
               space = "free") + 
    scale_size_continuous(name = "Sig. Metabolites",limits = c(1, 25), 
                          breaks = c(5,10,15,20)) +
    theme(panel.background = element_rect(fill="grey95"), 
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom",
          legend.justification = c("center","bottom"),
          strip.text.y = element_text(angle = 0, hjust = 0)) +
    scale_shape_manual(values = c(19,15), guide = "none") + 
    scale_x_continuous(limits = c(0, 6)) +
    guides(colour="none") +
    xlab("-log P") +
    ylab(NULL)
}

# (Min is 1); max is 24
min(mum_pw_long_final$data[[1]]$hits_sig)
max(mum_pw_long_final$data[[1]]$hits_sig)
# max(mum_pw_long4$data[[2]]$hits_sig)
# max(mum_pw_long4$data[[3]]$hits_sig)


# Run Figures
(fig_1_mum_bubble <- bubbleplot_fxn(mum_pw_long_final, "all_pfas")) #+ theme(legend.position = "none"))
# (pfsas <-    bubbleplot_fxn(mum_pw_long4, "pfsas") + theme(legend.position = "none"))
# (pfcas <-    bubbleplot_fxn(mum_pw_long4, "pfcas"))


# Save Fig --------------------
ggsave(fig_1_mum_bubble,
       filename = fs::path(dir_reports,
                           "Figure 1 Mummichog bubble plot PFAS Mixtures hyper g all_pfas mum_p05.jpeg"),
       width = 14, height = 6)
