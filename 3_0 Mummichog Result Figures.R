# 2) Plot Mummichog Pathway Results
library(colorspace)
library(janitor)
exposure_type = "PFAS"

mum_pw_wide <- read_rds(
  fs::path(dir_temp, 
             exposure_type, 
             "mum_pathway_results", 
             "SOL CHS PFAS Mummichog wide sig PW.RDS")) %>% 
  clean_names() %>% 
  mutate(q_meta = p.adjust(fet_meta), 
         sig_fdr = if_else(q_meta < 0.2, "Sig. FDR < 0.2", "Not. Sig"))

# #Pivot longer on all values
longer_df1 <- mum_pw_wide %>% 
  select(pfas, path, path_2, super_pathway, 
         functional_groupings, sig,sig_fdr,q_meta,
         fet_solar, fet_chs, fet_meta,
         enrichment_solar, enrichment_chs, enrichment_meta,
         neg_logp_solar,neg_logp_chs, neg_logp_meta) %>% 
  pivot_longer(names_to = "name", values_to = "value", 
               cols = fet_solar:neg_logp_meta) 

# Final Dataset for Bubbleplot
mum_pw_long <- longer_df1 %>% 
  mutate(cohort = case_when(str_detect(name, "solar") ~ "SOLAR", 
                            str_detect(name, "chs") ~ "CHS", 
                            str_detect(name, "meta") ~ "Combined") %>% 
           as.factor() %>% 
           fct_relevel("SOLAR", "CHS", "Combined"), 
         name = str_remove(name, "_solar") %>% 
           str_remove("_chs") %>% 
           str_remove("_meta")) %>% 
  pivot_wider(names_from = name, values_from = value) 





# Create dataset for first two figs ---------------------------------------



# Datasets for first two figures
# Select only overlapping significant pathways
mum_pw_w_onlysig <- mum_pw_wide %>% 
  group_by(pfas) %>% 
  slice_min(order_by = fet_meta, n = 3) %>% 
  slice_max(order_by = enrichment_meta, n = 3) %>%  
  filter(sig_meta == "Sig.") %>%
  mutate(label = path_2) %>% 
  ungroup()

# Final Dataset for Plots
mum_pw_wide2 <- left_join(mum_pw_wide, mum_pw_w_onlysig) %>% 
  mutate(label = if_else(is.na(label), "", label))

## Figure without Pathways Labeled ---------------------
# (pfas_mum_pathwayfig_nolabel <- ggplot(mum_pw_wide2,
#                                        aes(x = neg_logp_solar,
#                                            y = neg_logp_chs,
#                                            color = neg_logp_meta,
#                                            size = enrichment_meta,
#                                            shape = sig_meta,
#                                            group = path)) +
#    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey20") +
#    geom_point() +
#    scale_color_continuous_sequential(palette = "Oslo", rev = FALSE, begin = .1,
#                                      name = "Meta Analysis -Log P" )+
#    scale_size_continuous(name = "Weighted Enrichment", range = c(1, 4)) +
#    facet_wrap(~pfas) +
#    scale_shape_manual(values = c(16, 17), name = "Meta Analysis p<0.01")+
#    # theme(legend.title = element_blank()) +
#    xlab("-log P (SOLAR)") +
#    ylab("-log P (CHS)"))
# 
# # Save Fig
# ggsave(pfas_mum_pathwayfig_nolabel,
#        filename = fs::path(dir_reports,
#                            "Mummichog Pathways All PFAS SOL CHS no label.jpeg"),
#        width = 15, height = 11)
# 

# Figure with Pathways Labeled ---------------------------
# (pfas_mum_pathwayfig <- ggplot(mum_pw_wide2,
#                                aes(x = neg_logp_solar,
#                                    y = neg_logp_chs,
#                                    color = neg_logp_meta,
#                                    size = enrichment_meta,
#                                    shape = sig_meta,
#                                    group = path)) +
#    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey20") +
#    geom_point() +
#    geom_label_repel(data = mum_pw_w_onlysig, aes(label = label),
#                     color = "black",
#                     size = 3,
#                     box.padding = .75,
#                     force = 1,
#                     min.segment.length = 0,
#                     max.overlaps = 5) +
#    scale_color_continuous_sequential(palette = "Oslo",
#                                      rev = FALSE,begin = .1,
#                                      name = "Meta Analysis -Log P" )+
#    scale_size_continuous(name = "Weighted Enrichment", range = c(1, 4)) +
#    facet_wrap(~pfas) +
#    scale_shape_manual(values = c(16, 17), name = "Meta Analysis p<0.01")+
#    xlab("-log P (SOLAR)") +
#    ylab("-log P (CHS)"))
# 
# 
# # Save Fig
# ggsave(pfas_mum_pathwayfig,
#        filename = fs::path(dir_reports,
#                            "Mummichog Pathways All PFAS SOL CHS.jpeg"),
#        width = 15, height = 11)



# Bubble plot ---------------------------------------
mum_pw_long2 <- mum_pw_long %>%
  filter(sig == "Sig. Both Cohorts") %>%
  arrange(super_pathway, q_meta) %>%
  mutate(path_2 = fct_inorder(path_2))


(pfas_pw_bubbleplot <- ggplot(mum_pw_long2,
                              aes(x = neg_logp,
                                  y = fct_rev(path_2),
                                  color = super_pathway,
                                  size = enrichment)) +
    geom_point(alpha = .75) +
    facet_grid(pfas~cohort, scales = "free_y", space = "free") + 
    scale_size(name = "Enrichment") +
    
    # geom_vline(xintercept = -log(2.109881e-04), linetype = 2) +
    geom_vline(xintercept = 8.5) +
    theme(panel.background = element_rect(fill="grey95"), 
          strip.background = element_rect(fill = "white"), 
          strip.text.y = element_text(angle = 0, hjust = 0)) +
    xlab("-log P") +
    ylab(NULL)+
    scale_color_discrete_qualitative(name = "Super pathway"))



# Save Fig
ggsave(pfas_pw_bubbleplot,
       filename = fs::path(dir_reports,
                           "Mummichog buubble plot.jpeg"),
       width = 14, height = 11)

