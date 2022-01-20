# 2) Plot Mummichog Pathway Results
library(colorspace)
library(janitor)
exposure_type = "PFAS_Mixtures"

# Read in data
mum_pw_long <- read_rds(
  fs::path(dir_results_mum_mixtures, 
           "Mixture effect hyper_g", 
           "SOL CHS PFAS Mummichog long sig PW.RDS"))


# Final Dataset for Bubbleplot: format "name" to remove _sol/_chs/_meta, then pivot wider
mum_pw_long <- mum_pw_long %>% 
  mutate(cohort = case_when(str_detect(cohort, "solar") ~ "SOLAR", 
                            str_detect(cohort, "chs") ~ "CHS", 
                            str_detect(cohort, "meta") ~ "Combined") %>% 
           as.factor() %>% 
           fct_relevel("SOLAR", "CHS", "Combined"))


# Bubble plot ----------------------------------------------------
mum_pw_long2 <- mum_pw_long %>%
  # filter(sig_overall_p == "Sig.") #%>%
  filter(sig == "Sig. Both Cohorts") #%>%
# arrange(super_pathway, q_meta) %>%
# mutate(path_2 = fct_inorder(path_2))
table(mum_pw_long$sig)
mum_pw_long2

bubbleplot_fxn <- function(data, mixture_name){
  ggplot(data %>% filter(mixture == mixture_name),
         aes(x = neg_logp,
             y = fct_rev(path_2),
             color = super_pathway,
             size = hits_sig)) +
    geom_point(alpha = .75) +
    geom_vline(xintercept = -log10(0.05),
               color = "grey50", 
               linetype = 2) + 
    facet_grid(super_pathway~cohort, 
               scales = "free_y", 
               space = "free") + 
    scale_size(name = "Sig. Metabolites") +
    theme(panel.background = element_rect(fill="grey95"), 
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0, hjust = 0)) +
    guides(colour="none") +
    xlab("-log P") +
    ylab(NULL)
}



all_pfas <- bubbleplot_fxn(mum_pw_long2, "all_pfas")
pfsas <- bubbleplot_fxn(mum_pw_long2, "pfsas")
pfcas <- bubbleplot_fxn(mum_pw_long2, "pfcas")



# Save Figs
ggsave(all_pfas,
       filename = fs::path(dir_reports,
                           "mummichog_pw_bubbleplots",
                           "Mummichog bubble plot PFAS Mixtures hyper g all_pfas.jpeg"),
       width = 14, height = 6)

ggsave(pfsas,
       filename = fs::path(dir_reports,
                           "mummichog_pw_bubbleplots",
                           "Mummichog bubble plot PFAS Mixtures hyper g pfsas.jpeg"),
       width = 14, height = 6)

ggsave(pfcas,
       filename = fs::path(dir_reports,
                           "mummichog_pw_bubbleplots",
                           "Mummichog bubble plot PFAS Mixtures hyper g pfcas.jpeg"),
       width = 14, height = 6)
