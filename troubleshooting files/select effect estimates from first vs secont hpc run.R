
# Version 4 results
sol_mwas4 <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "sol_pfas_mixtures_results_final_v4.csv")) %>% 
  as_tibble() 


chs_mwas4 <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "chs_pfas_mixtures_results_final_v4.csv")) %>% 
  as_tibble()



# Version 3 results
sol_mwas3 <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "sol_pfas_mixtures_results_v3.csv")) %>% 
  as_tibble()  %>% 
  mutate(exposure = str_replace(exposure, "Mixture effect", "mixture"))


chs_mwas3 <- read_csv(
  file = fs::path(dir_results, 
                  'PFAS_Mixtures', 
                  "chs_pfas_mixtures_results_final_v3.csv")) %>% 
  as_tibble() %>% 
  mutate(exposure = str_replace(exposure, "Mixture effect", "mixture"))


# Join dfs 
solar_new_old <- tidylog::left_join(sol_mwas4, 
                                    sol_mwas3,
                                    by = c("mode", "feature", "exposure"),
                                    suffix = c("_v4", "_v3"))


# Mutate data 
solar_new_old <- solar_new_old %>%
  mutate(sig_both = case_when(p_value_v4 < 0.05 & p_value_v3 < 0.05 ~ "both sig",
                              p_value_v4 < 0.05 | p_value_v3 < 0.05 ~ "one sig",
                              TRUE ~ "not sig"), 
         delta_ee = abs(estimate_beta-estimate), 
         mean_ee = (estimate+estimate_beta)/2, 
         mean_abs_ee = (abs(estimate)+abs(estimate_beta))/2) %>% 
  select(mode, feature, exposure, mean_ee, mean_abs_ee, delta_ee, estimate_beta, estimate,
         everything())


table(solar_new_old$p_value_v3 < 0.05, solar_new_old$p_value_v4<0.05)

# Plot all effecct est.
# (first_v_second <- ggplot(solar_new_old %>% filter(exposure == "mixture"), 
#                           aes(x = estimate_beta, y = estimate, color = delta_ee)) + 
#     geom_point(alpha = 0.5, size = .75) + 
#     facet_wrap(~exposure) +
#     scale_color_distiller(type = "seq", palette = 1, direction = 1) +
#     ylab("effect estimate second run") + 
#     xlab("effect estimate first run"))

# ggsave(first_v_second, filename = here::here("first_versus_second_run.jpeg"), 
# height = 8, width = 10)


# Select metabolites with lowest mean effect, highest mean effect, and largest delta
lowest_mean <- solar_new_old %>% 
  filter(exposure == "mixture") %>%
  slice_min(n = 42, order_by = abs(mean_abs_ee)) %>% 
  mutate(cat = "lowest mean effect")
  

highest_mean <- solar_new_old %>% 
  filter(exposure == "mixture") %>%
  slice_max(n = 42, order_by = abs(mean_ee)) %>% 
  mutate(cat = "highest mean effect")

# Most discordant
most_discordant <- solar_new_old %>% 
  filter(exposure == "mixture") %>%
  slice_max(n = 42, order_by = delta_ee) %>% 
  mutate(cat = "most discordant")

# Combine 
met_for_followup <- bind_rows(lowest_mean, highest_mean, most_discordant)


# Plot
ggplot(met_for_followup, 
                          aes(x = estimate_beta,
                              y = estimate, 
                              color = cat)) + 
    geom_point(alpha = 0.75) + 
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +
    facet_wrap(~exposure) +
    ylab("effect estimate second run") + 
    xlab("effect estimate first run")


# Save 
write_rds(met_for_followup,fs::path(dir_))