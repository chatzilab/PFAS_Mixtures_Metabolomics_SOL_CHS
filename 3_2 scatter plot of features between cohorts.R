

# Read in MWAS Beta Coefficients --------------------
mwas_beta_coefs <- read_rds(
   fs::path(dir_results, 
            "PFAS", 
            "SOL CHS all MWAS results long.rds")) %>% 
   bind_rows(.id = "cohort") 

# Pivot wider
mwas_beta_wide_chrt <- mwas_beta_coefs %>% 
   select(cohort, exposure, name, estimate, 
          p_value, q_value) %>% 
   pivot_wider(names_from = cohort, values_from = c(estimate, p_value, q_value ))


# Feature metaanalysis
wgt_sol = sqrt(310)
wgt_chs = sqrt(136)

mwas_beta_wide_chrt2 <- mwas_beta_wide_chrt %>%
   filter(!is.na(p_value_solar), !is.na(p_value_chs)) %>%
   mutate(sig = case_when(q_value_solar < 0.05 & q_value_chs < 0.05 ~ "Sig. Both Cohorts", 
                          q_value_solar < 0.05 ~ "Sig. SOLAR Only", 
                          q_value_chs   < 0.05 ~ "Sig. CHS Only", 
                          TRUE ~ "Not Significant"), 
          pfas = rename_pfas(str_remove(exposure, "lg2_"), include_asterisk = TRUE)) %>% 
   rowwise() %>%
   mutate(p_value_meta = metap::sumz(p = c_across(p_value_solar:p_value_chs),
                                     weights = c(wgt_sol, wgt_chs))$p) %>%
   ungroup() %>% 
   mutate(enrichment_meta =
             ((estimate_solar*wgt_sol)+(estimate_chs*wgt_chs))/(wgt_sol+wgt_chs),
          neg_logp_meta = -log(p_value_meta),
          sig_meta = if_else(p_value_meta < 0.01, "Sig.", "Not Sig."))




mwas_beta_wide_chrt3 <- mwas_beta_wide_chrt2 %>% 
   group_by(sig, pfas) %>% 
   slice_sample(n = 500) %>% 
   ungroup()




## SOLAR_CHS HCB -------------------
figure1 <- ggplot(mwas_beta_wide_chrt2,
                  aes(x = estimate_solar,
                      y = estimate_chs,
                      color = sig,
                      size = sig,
                      # group = cpd,
                      # shape = sig,
                      # alpha = sig
                  )) +
   geom_hline(yintercept = 0, color = "grey40") + 
   geom_vline(xintercept = 0, color = "grey40") + 
   geom_point() +
   scale_size_manual(values = c(.2, 1,.2,.2)) +
   scale_alpha_manual(values = c(.1,.7, .2, .2)) +
   # geom_hline(yintercept = -log(0.05),
   #            color = "grey80",
   #            size = .5, linetype =2) +
   
   # geom_text(aes(label = cpd_sig)) +
   # geom_smooth(aes(x = beta_solar,
   #                 y = beta_chs), inherit.aes = FALSE) +
   facet_wrap(~pfas, scales = "free") + 
   xlab("Beta SOLAR") + 
   ylab("Beta CHS") + 
   labs(color='FDR Q < 0.05', 
        size = 'FDR Q < 0.05', 
        alpha = 'FDR Q < 0.05')


ggsave(figure1, 
       filename = fs::path(dir_reports, "All Features Effect est scatterplots.jpg"), 
       height = 8, width = 10)




temp <- filter(mwas_beta_wide_chrt2, sig == "Sig. Both Cohorts")
