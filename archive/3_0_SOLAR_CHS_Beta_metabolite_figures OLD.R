# All Metabolites

# Read in SOLAR and CHS mz/cpd key and mwas results
chs_single_matches_w_mwas <- readRDS(fs::path(dir_temp, 
                                                "chs_mz_cpd_pathway_key_w_mwas.rds"))


solar_single_matches_w_mwas <- readRDS(fs::path(dir_temp, 
                                                  "solar_mz_cpd_pathway_key_w_mwas.rds"))

single_matches_w_mwas_final <- inner_join(solar_single_matches_w_mwas, 
                                          chs_single_matches_w_mwas, 
                                          by = c("cpd", "exposure"
                                                 # "aa_cho_lipid"
                                          ), 
                                          suffix = c("_solar", "_chs")) %>% 
  mutate(sig = case_when((p_value_solar < 0.05) & 
                           (p_value_chs  < 0.05) ~ "Jointly Significant", 
                         p_value_solar < 0.05  ~ "Sig. SOLAR only",
                         p_value_chs   < 0.05  ~ "Sig. CHS only"  ,
                         TRUE ~ "Not-sig"))


# Coefficient plots -------
# for(i in exposures){
#   solar_exposure_data <- solar_single_matches_w_mwas %>% 
#     filter(exposure == i) %>% 
#     group_by(cpd) %>% 
#     filter(p_value == min(p_value)) %>% 
#     ungroup() %>% 
#     mutate(cpd = reorder(cpd, beta), 
#            main_super_pathway_categorized = 
#              fct_lump(f = super_pathway)) 
#   
#   ## Figure
#   solar_fig <- ggplot(solar_exposure_data,
#                       aes(x = cpd,
#                           y = beta,
#                           color = main_super_pathway_categorized)) +
#     geom_pointrange(aes(ymin = conf_low, 
#                         ymax = conf_high), 
#                     alpha = 0.85,
#                     size = 0.2, fatten = .75 ) +
#     geom_hline(yintercept = 0) +
#     theme(axis.text.x = element_blank(), 
#           # legend.position = "none", 
#           strip.background = element_blank(),
#           strip.text.x = element_blank()) + 
#     facet_grid(~main_super_pathway_categorized, 
#                scales = "free_x", space = "free_x") + 
#     ggtitle(paste("SOLAR, beta [95% ci] between", 
#                   i, 
#                   "and metabolites", sep = " "))
#   
#   
#   ggsave(solar_fig, 
#          filename = here::here("figures", 
#                                paste0("SOLAR_",i,"_metabolites.jpg")),
#          width = 11, height = 5)
#   
#   
#   
#   ## CHS   -------------------
#   chs_exposure_data <- chs_single_matches_w_mwas %>% 
#     filter(exposure == exposure[1]) %>% 
#     group_by(cpd) %>% 
#     filter(p_value == min(p_value)) %>% 
#     ungroup() %>% 
#     mutate(cpd = reorder(cpd, beta), 
#            main_super_pathway_categorized = 
#              fct_lump(f = super_pathway)) 
#   
#   ## Figure
#   chs_fig <- ggplot(chs_exposure_data,
#                     aes(x = cpd,
#                         y = beta,
#                         color = main_super_pathway_categorized)) +
#     geom_pointrange(aes(ymin = conf_low, 
#                         ymax = conf_high), 
#                     alpha = 0.85,
#                     size = 0.2, fatten = .75 ) +
#     geom_hline(yintercept = 0) +
#     theme(axis.text.x = element_blank(), 
#           # legend.position = "none", 
#           strip.background = element_blank(),
#           strip.text.x = element_blank()) + 
#     facet_grid(~main_super_pathway_categorized, 
#                scales = "free_x", space = "free_x") + 
#     ggtitle(paste("CHS, beta [95% ci] between", i, "and metabolites", sep = " "))
#   
#   
#   ggsave(chs_fig, 
#          filename = here::here("figures", 
#                                paste0("CHS_",
#                                       i,
#                                       "_metabolites.jpg")),
#          width = 11, height = 5)
# }


write_csv(single_matches_w_mwas_final, "singlenmatches.csv")
# Scatterplot of betas between matching compounds in SOLAR and CHS ---------------

## SOLAR_CHS HCB -------------------
(figure1 <- ggplot(single_matches_w_mwas_final,
                   aes(x = beta_solar,
                       y = beta_chs,
                       color = aa_cho_lipid_solar, 
                       group = cpd,
                       shape = sig,
                       alpha = sig
                   )) +
   geom_hline(yintercept = 0) + 
   geom_vline(xintercept = 0) + 
   scale_shape_manual(values = c(18, 20, 20, 20)) +
   scale_alpha_manual(values = c(1,.1, .1, .1)) +
   # geom_hline(yintercept = -log(0.05),
   #            color = "grey80",
   #            size = .5, linetype =2) +
   geom_point(size = 1.2) +
   # geom_text(aes(label = cpd_sig)) +
   # geom_smooth(aes(x = beta_solar,
   #                 y = beta_chs), inherit.aes = FALSE) +
   facet_wrap(~exposure2) + 
   xlab("Beta CHS") + 
   ylab("Beta Solar"))

ggsave(figure1,  
       filename = here::here("figures",  "matching_cpd_betas.jpg"),
       width = 11, height = 8)

plotly::ggplotly(figure1)

#######


(figure1_hcb <- ggplot(single_matches_w_mwas_final %>% 
                         filter(exposure2 == "HCB"),
                       aes(x = beta_solar,
                           y = beta_chs,
                           color = aa_cho_lipid_solar+, 
                           group = cpd, 
                           shape = sig, 
                           size  = sig
                       )) +
   geom_hline(yintercept = 0) + 
   geom_vline(xintercept = 0) + 
   geom_point() + 
   scale_shape_manual(values = c(16,20,15,17)) +
   scale_size_manual(values = c(1.6, 1, 1, 1)) +
   # geom_text(aes(label = cpd_sig)) +
   facet_wrap(~exposure2) + 
   xlab("Beta CHS") + 
   ylab("Beta Solar")) 

ggsave(figure1_hcb,  
       filename = here::here("figures",  
                             "matching_cpd_betas hcb.jpg"),
       width = 11, height = 8)







#### Type 2 manhattan plots  ----------------

ggplot(single_matches_w_mwas_final %>% 
         filter( exposure2 == "HCB"),
       aes(x = rt_solar,
           y = -log(p_value_solar),
           color = aa_cho_lipid_solar, 
           shape = aa_cho_lipid_solar 
           )) +
  geom_point(size = 2)+
  
  geom_hline(yintercept = -log(0.05),
             color = "grey80",
             size = .5, linetype =2)
