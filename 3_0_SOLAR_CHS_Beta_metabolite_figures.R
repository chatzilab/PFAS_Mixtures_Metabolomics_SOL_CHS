# All Metabolites

# Read in SOLAR and CHS mz/cpd key and mwas results
chs_single_matches_w_mwas <- readRDS(here::here("Temporary results", 
                                                "chs_mz_cpd_pathway_key_w_mwas.rds"))

solar_single_matches_w_mwas <- readRDS(here::here("Temporary results", 
                                                  "solar_mz_cpd_pathway_key_w_mwas.rds"))

single_matches_w_mwas_final <- inner_join(solar_single_matches_w_mwas, 
                                          chs_single_matches_w_mwas, 
                                          by = c("cpd", "exposure", "exposure2"), 
                                          suffix = c("_solar", "_chs")) %>% 
  mutate(sig = if_else((p_value_solar < 0.05) & 
                         (p_value_chs < 0.05), 
                       "Sig", 
                       "Not-sig"))



# Coefficient plots -------
for(i in exposures){
  solar_exposure_data <- solar_single_matches_w_mwas %>% 
    filter(exposure == i) %>% 
    group_by(cpd) %>% 
    filter(p_value == min(p_value)) %>% 
    ungroup() %>% 
    mutate(cpd = reorder(cpd, beta), 
           main_super_pathway_categorized = 
             fct_lump(f = super_pathway)) 
  
  ## Figure
  solar_fig <- ggplot(solar_exposure_data,
                      aes(x = cpd,
                          y = beta,
                          color = main_super_pathway_categorized)) +
    geom_pointrange(aes(ymin = conf_low, 
                        ymax = conf_high), 
                    alpha = 0.85,
                    size = 0.2, fatten = .75 ) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_blank(), 
          # legend.position = "none", 
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    facet_grid(~main_super_pathway_categorized, 
               scales = "free_x", space = "free_x") + 
    ggtitle(paste("SOLAR, beta [95% ci] between", 
                  i, 
                  "and metabolites", sep = " "))
  
  
  ggsave(solar_fig, 
         filename = here::here("figures", 
                               paste0("SOLAR_",i,"_metabolites.jpg")),
         width = 11, height = 5)
  
  
  
  ## CHS HCB -------------------
  chs_exposure_data <- chs_single_matches_w_mwas %>% 
    filter(exposure == exposure[1]) %>% 
    group_by(cpd) %>% 
    filter(p_value == min(p_value)) %>% 
    ungroup() %>% 
    mutate(cpd = reorder(cpd, beta), 
           main_super_pathway_categorized = 
             fct_lump(f = super_pathway)) 
  
  ## Figure
  chs_fig <- ggplot(chs_exposure_data,
                    aes(x = cpd,
                        y = beta,
                        color = main_super_pathway_categorized)) +
    geom_pointrange(aes(ymin = conf_low, 
                        ymax = conf_high), 
                    alpha = 0.85,
                    size = 0.2, fatten = .75 ) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_blank(), 
          # legend.position = "none", 
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    facet_grid(~main_super_pathway_categorized, 
               scales = "free_x", space = "free_x") + 
    ggtitle(paste("CHS, beta [95% ci] between", i, "and metabolites", sep = " "))
  
  
  ggsave(chs_fig, 
         filename = here::here("figures", 
                               paste0("CHS_",
                                      i,
                                      "_metabolites.jpg")),
         width = 11, height = 5)
}



# Scatterplot of betas between matching compounds in SOLAR and CHS ---------------

## SOLAR_CHS HCB -------------------
(figure1 <- ggplot(single_matches_w_mwas_final,
                   aes(x = beta_solar,
                       y = beta_chs,
                       color = super_pathway_solar, 
                       group = cpd, 
                       shape = sig, 
                       alpha = sig)) +
   geom_hline(yintercept = 0) + 
   geom_vline(xintercept = 0) + 
   # geom_hline(yintercept = -log(0.05),
   #            color = "grey80",
   #            size = .5, linetype =2) +
   geom_point(size = 1.2) +
   # geom_text(aes(label = cpd_sig)) +
   facet_wrap(~exposure2) + 
   xlab("Beta CHS") + 
   ylab("Beta Solar"))

ggsave(figure1,  
       filename = here::here("figures",  "matching_cpd_betas.jpg"),
       width = 11, height = 8)

#######


(figure1_hcb <- ggplot(single_matches_w_mwas_final %>% 
                         filter(exposure2 == "HCB"),
                       aes(x = beta_solar,
                           y = beta_chs,
                           color = super_pathway_solar, 
                           group = cpd, 
                           shape = sig, 
                           alpha = sig)) +
   geom_hline(yintercept = 0) + 
   geom_vline(xintercept = 0) + 
   # geom_hline(yintercept = -log(0.05),
   #            color = "grey80",
   #            size = .5, linetype =2) +
   geom_point(size = 1.2) +
   # geom_text(aes(label = cpd_sig)) +
   facet_wrap(~exposure2) + 
   xlab("Beta CHS") + 
   ylab("Beta Solar")) 

ggsave(figure1_hcb,  
       filename = here::here("figures",  
                             "matching_cpd_betas hcb.jpg"),
       width = 11, height = 8)







#### Betas between SOLAR and CHS Metabolites ----------------

ggplot(chs_single_matches_w_mwas,
       aes(x = beta,
           y = -log(p_value),
           color = super_pathway)) +
  # geom_hline(yintercept = -log(0.05),
  #            color = "grey80",
  #            size = .5, linetype =2) +
  geom_point(size = .75, shape = 21) +
  # geom_text(aes(label = cpd_sig)) + 
  facet_wrap(~exposure2)
xlab("Compound Number (ordered by p-value within superpathway)") + 
  ylab("-log P")
