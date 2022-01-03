# Volcano Plots

# SOLAR Volcano Plots ----------------------------------------------
sol_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "sol_pfas_mixtures_results_w_09.csv")) %>% 
  as_tibble()

# Get reduced dataset
sol_mwas_reduced <- sol_mwas %>% 
  filter(p_value  < 0.99)

# Volcano Plot
(solar_volcano_plot <- ggplot(sol_mwas_reduced,
                              aes(x = estimate_beta, y = -log(p_value))) +
    geom_point(size = 1, alpha = 0.5) + 
    geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
    facet_wrap(~exposure, scales = "free") +
    ylab("-log10 p") +
    xlab("Effect Estimate") +
    ggtitle("SOLAR"))

# Save and clean env.
ggsave(solar_volcano_plot,
       filename =  fs::path(dir_reports,
                            "Volcano plots",
                            "SOLAR Mixtures analysis volcano plots_p_original_w_09.jpg"), 
       width = 6, height = 5)
rm(solar_volcano_plot)


# CHS Volcano Plots ----------------------------------------------
chs_mwas <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "chs_pfas_mixtures_results_w_09.csv")) 

# Get reduced dataset
chs_mwas_reduced <- chs_mwas %>% 
  filter(p_value  < 0.99)

# Volcano Plot
chs_volcano_plot <- ggplot(chs_mwas_reduced,
                           aes(x = estimate_beta, y = -log(p_value))) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  facet_wrap(~exposure, scales = "free") +
  ylab("-log10 p") +
  xlab("Effect Estimate") +
  ggtitle("CHS")

# # Save results
ggsave(chs_volcano_plot,
       filename =  fs::path(
         dir_reports,
         "Volcano plots",
         "chs Mixtures analysis volcano plots_p_from_original_w_09.jpg"),
       width = 6.5, height = 5)