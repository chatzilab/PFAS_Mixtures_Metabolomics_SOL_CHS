library(UpSetR)


# read data from 2_5 
upplot_data <- read_rds(
  here::here(
    dir_results_mixtures,
    "Ecd effect estimates and pathways sol chs for upset.RDS"))

# SOLAR --------------------------------------------

# NEW SCRIPT FEB 22 2022: -----------------------------
ecds_for_upset <- upplot_data %>% 
  group_by(met_name) %>% 
  mutate(empirical_compound_all = str_c(empirical_compound, collapse = ";"), 
         empirical_compound = empirical_compound[1]) %>% 
  ungroup() %>% 
  tidylog::filter(
    # From Doug: This is a weird adduct that probably is not correct
    matched_form != "M-HCOOH+H[1+]", 
    # This is not homovanilian, since it does not match the other 3 adduct retention times:
    !(met_name == "Homovanillin" & name == "185.080994057582_72.2749441021845")) 


# -------------------------------------------------------------------
# Calculate max estimate and sig
ecds_for_upset_sol <- ecds_for_upset %>% 
  filter(p_value_sol_mixture < 0.05) %>%
  group_by(pathway, empirical_compound) %>% 
  summarise(max_beta = max(estimate_beta_sol_mixture), 
            min_beta = min(estimate_beta_sol_mixture), 
            max_est = if_else(abs(max_beta) > abs(min_beta), 
                              max_beta, 
                              min_beta), 
            ec = str_c(empirical_compound, collapse = ";"), 
            met_name = str_c(met_name, collapse = ";")) %>% 
  ungroup()

# Pivot wider to get table of 0/1s
upplot_data_sol <- ecds_for_upset_sol %>% 
  select(-contains("beta")) %>%
  mutate(inpath = 1, 
         pathway = str_replace_all(pathway, 
                                   "Metabolism", 
                                   "metabolism") ) %>% 
  pivot_wider(values_from = inpath, 
              names_from = pathway,
              id_cols = c(empirical_compound, met_name, ec, "max_est")) %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  as.data.frame()

## UpSet graph showing the number of empirical cpds from each pathway --------------
upset_sol <- upset(upplot_data_sol,
                   sets =c("Tyrosine metabolism", 
                           "Urea cycle/amino group metabolism",
                           "Arginine and Proline metabolism",
                           "Lysine metabolism", 
                           "Glutathione metabolism", 
                           "Porphyrin metabolism",
                           "De novo fatty acid biosynthesis"),
                   keep.order = TRUE, 
                   text.scale = 1.4, 
                   mb.ratio = c(2/3, 1/3), # Ratio of size between bar graphs and matrix graph
                   show.numbers = FALSE,
                   mainbar.y.label = "Metabolite Intersections", 
                   sets.x.label = "Sig. Metabolites", 
                   nsets = length(unique(ecds_for_upset_sol$pathway)),
                   order.by = "freq", # either degree or freq
                   queries = list(
                     list(
                       query = elements, 
                       params = list(
                         "empirical_compound", 
                         unique(upplot_data_sol$empirical_compound)), 
                       color = "purple", active = TRUE),
                     list(
                       query = elements, 
                       params = list(
                         "empirical_compound", 
                         unique(upplot_data_sol$empirical_compound[
                           upplot_data_sol$max_est>0
                         ])), 
                       color = "orange", active = TRUE))); upset_sol

#Save plot
jpeg(file = here::here(dir_reports, "SOLAR Metabolite pathway upset.jpeg"), 
     width = 6, height = 6, units = "in", res = 350) 
upset_sol
dev.off()


# CHS --------------------------------------------
ecds_for_upset_chs <- ecds_for_upset %>% 
  tidylog::filter(p_value_chs_mixture < 0.05 | 
                    lcl_beta_chs_mixture == 0 | 
                    ucl_beta_chs_mixture == 0) %>%
  group_by(pathway, empirical_compound) %>% 
  summarise(max_beta = max(estimate_beta_chs_mixture), 
            min_beta = min(estimate_beta_chs_mixture), 
            max_est = if_else(abs(max_beta) > abs(min_beta), 
                              max_beta, 
                              min_beta)) %>% 
  ungroup()

# Pivot wider to get table of 0/1s
upplot_data_chs <- ecds_for_upset_chs %>% 
  select(-contains("beta")) %>%
  mutate(inpath = 1, 
         pathway = str_replace_all(pathway, "Metabolism", "metabolism")) %>% 
  pivot_wider(values_from = inpath, 
              names_from = pathway,
              id_cols = c(empirical_compound, "max_est")) %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  as.data.frame()


## UpSet graph showing the number of empirical cpds from each pathway --------------
upset_chs <- upset(upplot_data_chs,
                   sets =c("Tyrosine metabolism", 
                           "Urea cycle/amino group metabolism",
                           "Arginine and Proline metabolism",
                           "Lysine metabolism", 
                           "Glutathione metabolism", 
                           "Porphyrin metabolism",
                           "De novo fatty acid biosynthesis"),
                   keep.order = TRUE, 
                   text.scale = 1.4, 
                   mb.ratio = c(2/3, 1/3), # Ratio of size between bar and matrix
                   show.numbers = FALSE,
                   mainbar.y.label = "Metabolite Intersections", 
                   sets.x.label = "Sig. Metabolites", 
                   nsets = length(unique(ecds_for_upset_chs$pathway)),
                   order.by = "freq", # either degree or freq
                   queries = list(
                     list(
                       query = elements, 
                       params = list(
                         "empirical_compound", 
                         unique(upplot_data_chs$empirical_compound)), 
                       color = "purple", active = TRUE),
                     list(
                       query = elements, 
                       params = list(
                         "empirical_compound", 
                         unique(upplot_data_chs$empirical_compound[
                           upplot_data_chs$max_est>0
                         ])), 
                       color = "orange", active = TRUE)
                   )
); upset_chs

# Save figure
jpeg(file = here::here(dir_reports, "CHS Metabolite pathway upset.jpeg"), 
     width = 5.5, height = 6, units = "in", res = 350) 
upset_chs
dev.off()