# Volcano Plots

# SOLAR Volcano Plots ----------------------------------------------
sol_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "sol_pfas_mixtures_results_final_v3.csv")) %>% 
  as_tibble()

# Get reduced dataset
sol_mwas_reduced <- sol_mwas %>% 
  filter(p  < 0.99)

# Volcano Plot
solar_volcano_plot <- ggplot(sol_mwas_reduced,
                             aes(x = estimate, y = -log(p))) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  facet_wrap(~exposure, scales = "free") +
  ylab("-log10 p") +
  xlab("Effect Estimate") +
  ggtitle("SOLAR")
# 
# # Save results
ggsave(solar_volcano_plot,
       filename =  fs::path(dir_reports,
                            "Volcano plots",
                            "SOLAR Mixtures analysis volcano plots_p_original_v3.jpg"), 
       width = 6, height = 5)

# Volcano Plot
solar_volcano_plot_og_p <- ggplot(sol_mwas_reduced, 
                             aes(x = estimate, y = -log(p_value))) +
  geom_point(size = 1, alpha = 0.5) + 
  facet_wrap(~exposure, scales = "free") + 
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  ylab("-log10 p") + 
  xlab("Effect Estimate") + 
  ggtitle("SOLAR")

# Save results
ggsave(solar_volcano_plot_og_p,
       filename =  fs::path(dir_reports, 
                            "Volcano plots",
                            "SOLAR Mixtures analysis volcano plots_p_from_chisq_v3.jpg"),
       width = 6.5, height = 5)

rm(solar_volcano_plot)


# CHS Volcano Plots ----------------------------------------------
chs_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "chs_pfas_mixtures_results_final_v3.csv")) 


chs_mwas <- chs_mwas %>% 
  mutate(wald = (abs(estimate)/sd),
         p = 2*(1-pnorm(wald,0,1)),
         p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
         p_value = if_else(wald^2 > 2000, 4.6*(10^ (-256)), p_value),
         neg_log_p = -log10(p_value))


# Get reduced dataset
chs_mwas_reduced <- chs_mwas %>% 
  filter(p  < 0.99)

# Volcano Plot
chs_volcano_plot <- ggplot(chs_mwas_reduced,
                           aes(x = estimate, y = -log(p))) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  facet_wrap(~exposure, scales = "free") +
  ylab("-log10 p") +
  xlab("Effect Estimate") +
  ggtitle("CHS")
 
# # Save results
ggsave(chs_volcano_plot,
       filename =  fs::path(dir_reports,
                            "Volcano plots",
                            "chs Mixtures analysis volcano plots_p_from_original_v3.jpg"),
       width = 6.5, height = 5)

# Volcano Plot
chs_volcano_plot_og_p <- ggplot(chs_mwas_reduced, 
                                aes(x = estimate, y = -log(p_value))) +
  geom_point(size = 0.5, alpha = 0.5) + 
  facet_wrap(~exposure, scales = "free") + 
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  ylab("-log10 p") + 
  xlab("Effect Estimate") + 
  ggtitle("CHS")

# Save results
ggsave(chs_volcano_plot_og_p,
       filename =  fs::path(dir_reports, 
                            "Volcano plots",
                            "chs Mixtures analysis volcano plots_p_from_chisq_v3.jpg"),
       width = 6.5, height = 5)