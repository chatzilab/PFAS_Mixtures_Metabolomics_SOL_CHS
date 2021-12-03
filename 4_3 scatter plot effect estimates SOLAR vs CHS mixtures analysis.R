# Get most significant pathway from each PFAS
library(ggExtra)
library(ggrepel)
library(tidylog)
# 2) Plot Mummichog Pathway Results
library(colorspace)
library(janitor)

# Set vars
pfas_name <- "Mixture effect"
cohort_name <- "solar"
pw_name = "Tyrosine metabolism"


# 1) read in data  ------------------------
## Mummichog Data ------------------------------------
mum_pw_wide <- read_rds(
  fs::path(dir_results_mum_mixtures,
           "mum_pathway_results", 
           "SOL CHS PFAS Mummichog wide sig PW.RDS")) %>% 
  clean_names() %>% 
  mutate(q_meta = p.adjust(fet_meta), 
         sig_fdr = if_else(q_meta < 0.2, "Sig. FDR < 0.2", "Not. Sig"), 
         mean_num_sig = (hits_sig_solar + hits_sig_chs)/2)


## MWAS Data ------------------------------------
mwas_beta_coefs <- read_rds(
  fs::path(dir_results_mixtures, 
           "SOL CHS all Mixtures MWAS results long.rds")) %>% 
  modify(~.x %>% mutate(exposure2 = exposure %>% 
                          str_remove(., "lg2_") %>% 
                          rename_pfas(include_asterisk = TRUE)))


## Read in mzrt key--------------------------------------
mzrt_key <- read_rds(fs::path(dir_data_local,
                              "Supporting Files", 
                              "mummichog_pw_ec_feature_key_with_cpd_names.rds"))

# 2) Data Cleaning: filter and merge metabolites from specific pathways ----------------------
# Get mzrt key for specific pathway 
key_pw_met <- mzrt_key %>% 
  # filter(pathway == pw_name) %>% 
  select(mode, pathway, met_name, everything())

# Filter only primary ions; clean names of cpds with multiple names 
key_pw_met2 <- key_pw_met %>% 
  filter(matched_form == "M+H[1+]" | matched_form == "M-H[-]") %>%
  filter(!(matched_compound == met_name)) %>% # This drops metabolites without a name
  group_by(name) %>% 
  summarise(matched_compound = toString(unique(matched_compound)), 
            ec_mode = toString(unique(ec_mode)), 
            mean_mass = mean(query_mass), 
            sd_mass = sd(query_mass) %>% replace_na(., 0), 
            met_name_first_only  = toString(unique(met_name)[1]),
            met_name_all = toString(unique(met_name)), 
            met_name  = str_c(unique(met_name)[1], 
                              " (", 
                              length(unique(met_name)), 
                              ")") %>% 
              str_remove(" \\(1\\)") ) %>% 
  mutate(met_name_first_only = 
           case_when(
             str_detect(met_name_all,  
                        "3-methyl pyruvic acid, Acetoacetic acid") ~ "Acetoacetic acid", 
             str_detect(met_name_all,  
                        "Ascorbate, Glucuronolactone") ~ "Glucuronolactone", 
             str_detect(met_name_all, 
                        "Phosphorylcholine, Vanylglycol") ~ "Vanylglycol",
             str_detect(met_name_all, 
                        "Morphine, 2-Hydroxycarbamazepine") ~ "2-Hydroxycarbamazepine",
             TRUE ~ met_name_first_only))



# Get pathway specific effect estimates for main exposure
pw_eff_est_annotated_allfts <- mwas_beta_coefs %>%
  modify(~.x %>% 
           left_join(key_pw_met2, ., by = c("name" = "feature")))  %>% 
  bind_rows(.id = "cohort")


# select only a single feature from each compound
most_altered_fts <- pw_eff_est_annotated_allfts %>% 
  group_by(name, matched_compound, exposure, exposure2) %>% 
  summarise(estimate_mean = mean(estimate)) %>%
  # mutate(rownum = row_number()) %>%
  ungroup() %>% 
  # select(name, matched_compound, exposure, exposure2, rownum)
  group_by(matched_compound, exposure, exposure2) %>%
  filter(abs(estimate_mean) == max(abs(estimate_mean))) %>%
  ungroup()


# Join back with pw_eff_est_annotated_allfts
pw_eff_est_annotated <- left_join(most_altered_fts, 
                                  pw_eff_est_annotated_allfts, 
                                  by = c("name", 
                                         "matched_compound",
                                         "exposure", 
                                         "exposure2"))


# Pivot annotated dataset wider 
pw_eff_est_wide <- pw_eff_est_annotated %>%
  pivot_wider(., 
              names_from = cohort,
              values_from = c("estimate", "p_value", 
                              "q_value", "lcl", "ucl"),
              id_cols = c("name", "exposure", "exposure2",
                          "met_name",
                          "met_name_first_only",
                          "met_name_all", 
                          "matched_compound")) 


## Perform meta-analysis and calculate weighted mean effect estimate ------------
# set weights 
wgt_sol = sqrt(310); wgt_chs = sqrt(136)

# Run analysis
pw_eff_est_wide <- pw_eff_est_wide %>% 
  mutate(p_value_solar = if_else(p_value_solar > 0.999, 0.999, p_value_solar), 
         p_value_chs = if_else(p_value_chs > 0.999, 0.999, p_value_chs)) %>%
  rowwise() %>% 
  mutate(p_meta = metap::sumz(p = c_across(p_value_solar:p_value_chs), 
                              weights = c(wgt_sol, wgt_chs))$p, 
         neg_logp_meta = -log(p_meta), 
         sig_meta = if_else(p_meta < 0.05, "p<0.05", "p>0.05")) %>% 
  ungroup() 

# Create Significance variables 
pw_eff_est_wide <- pw_eff_est_wide %>% 
  mutate(sig = case_when(q_value_solar < 0.2 & q_value_chs < 0.2 ~ "q<0.2 both cohorts", 
                         p_value_solar < 0.05 & p_value_chs < 0.05 ~ 
                           "p<0.05 Both Cohorts",
                         p_value_solar < 0.05 ~ "p<0.05 (SOLAR only)",
                         p_value_chs < 0.05 ~ "p<0.05 (CHS only)",
                         TRUE ~ "Non. Sig"), 
         sig_solar = case_when(q_value_solar < 0.2 ~  "q<0.02", 
                               p_value_solar < 0.05 ~ "p<0.05",
                               # p_meta < 0.05 ~ "Meta-analysis: p<0.05",
                               TRUE ~ "Not Sig."),
         sig_chs = case_when(q_value_chs < 0.2 ~ "q<0.02", 
                             p_value_chs < 0.05 ~ "p<0.05",
                             # p_meta < 0.05 ~ "Meta-analysis: p<0.05",
                             TRUE ~ "Not Sig."),
         mean_effect = ((estimate_solar*wgt_sol)+
                          (estimate_chs*wgt_chs))/(wgt_sol+wgt_chs)) %>% 
  mutate(across(c(sig_solar, sig_chs), 
                ~fct_relevel(., "q<0.02", "p<0.05",
                             # "Meta-analysis: p<0.05",
                             "Not Sig.")))


# Create dataframe for single exposure of interest (mixture effect)
single_pfas_eff_est <- pw_eff_est_wide %>% 
  filter(exposure2 == pfas_name)



# Create figure of all metabolites -------------------------------------------
pw_eff_est_wide <- pw_eff_est_wide %>% 
  mutate(name_sig_only = if_else(sig == "q<0.2 both cohorts" | 
                                   sig == "p<0.05 Both Cohorts", 
                                 met_name_first_only, 
                                 NA_character_))


(figure_all_metabolites_coef_comp <- ggplot(pw_eff_est_wide, 
                                            aes(x = estimate_solar, 
                                                y = estimate_chs, 
                                                color = sig, 
                                                alpha = sig, 
                                                label = name_sig_only)) + 
    geom_point() + 
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey30") +
    facet_wrap(~exposure2,scales = "free") + 
    scale_alpha_manual(values = c(.6, .9,.9, 1,1)) +
    scale_color_manual(values = c("grey40", "green", "blue", "darkorange", "red")) +
    ylab("Beta (CHS)") + 
    xlab("Beta (SOLAR)"))


# Save
ggsave(figure_all_metabolites_coef_comp, 
       filename = fs::path(dir_reports, 
                           "SOL CHS Scatter plot comparison of all mixtures effect estimates.jpg"), 
       width = 10, height = 7, units = "in")

# Add labels
(fig_labs <- ggplot(pw_eff_est_wide, aes(x = estimate_solar, y = estimate_chs, 
                                         label = name_sig_only)) + 
    geom_point(aes(color = sig)) + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = .6) +
    geom_label_repel(size = 2.5) +
    facet_wrap(~exposure2, scales = "free") + 
    scale_alpha_manual(values = c(.6, .9,.9, 1,1), guide = "none") +
    scale_color_manual(values = c("grey40", "green", "blue", "darkorange", "red")) +
    ylab("Beta (CHS)") + 
    xlab("Beta (SOLAR)"))



ggsave(fig_labs, 
       filename = fs::path(dir_reports, 
                           "SOL CHS Scatter plot comparison of all mixtures effect estimates with labs.jpg"), 
       width = 10, height = 7, units = "in")

