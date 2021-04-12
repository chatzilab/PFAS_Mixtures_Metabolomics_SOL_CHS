xxx <- inner_join(solar_single_matches_w_mwas, 
                  chs_single_matches_w_mwas, 
                  by = c("cpd", "exposure", "exposure2"), suffix = c("_solar", "_chs"))


single_matches_w_mwas_final <- xxx %>% 
  filter(p_value_solar < 0.1, 
         p_value_chs < 0.1)

single_matches_w_mwas_final$super_pathway_solar


figure1 <- ggplot(single_matches_w_mwas_final,
       aes(x = beta_solar,
           y = beta_chs,
           color = super_pathway_solar, 
           group = cpd)) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  # geom_hline(yintercept = -log(0.05),
  #            color = "grey80",
  #            size = .5, linetype =2) +
  geom_point(size = 1.2) +
  # geom_text(aes(label = cpd_sig)) +
  facet_wrap(~exposure2) + 
  xlab("Beta CHS") + 
  ylab("Beta Solar")

library(plotly)

ggplotly(figure1)
