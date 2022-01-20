# Get most significant pathway from each PFAS
library(ggExtra)
library(tidylog)
# 2) Plot Mummichog Pathway Results
library(colorspace)
library(janitor)

# Set vars
# cohort_name <- "solar"
pw_name = "Bile acid biosynthesis"


source(here::here("0_0_1_format_vars_funs.R"))

# 1) read in data  ------------------------
## Mummichog Data ------------------------------------
mum_pw_wide <- read_rds(
  fs::path(dir_results_mum_mixtures,
           "Mixture effect hyper_g", 
           "SOL CHS PFAS Mummichog wide sig PW.RDS")) %>% 
  clean_names() %>% 
  mutate(q_meta = p.adjust(pval_meta), 
         sig_fdr = if_else(q_meta < 0.2, 
                           "Sig. FDR < 0.2",
                           "Not. Sig"), 
         mean_num_sig = (hits_sig_solar + hits_sig_chs)/2)


## MWAS Data ------------------------------------
mwas_beta_coefs <- read_rds(
  fs::path(dir_results_mixtures, 
           "SOL CHS all Mixtures MWAS results long hyper_g.rds")) %>% 
  modify(~.x %>% mutate(exposure = exposure %>% 
                          rename_pfas(include_asterisk = TRUE)))

## Read in mzrt key--------------------------------------
mzrt_key <- read_rds(fs::path(dir_data_local,
                              "Supporting Files", 
                              "mummichog_pw_ec_feature_key_with_cpd_names.rds"))

# 2) Data Cleaning: filter and merge metabolites from specific pathways ------------
# Get mzrt key for specific pathway 
key_pw_met <- mzrt_key %>% 
  filter(pathway == pw_name) %>% 
  select(pathway, met_name, everything())

# Filter only primary ions; clean names of cpds with multiple names 
key_pw_met2 <- key_pw_met %>% 
  filter(matched_form == "M+H[1+]" | matched_form == "M-H[-]") %>%
  filter(!(matched_compound == met_name)) %>% # This drops metabolites without a name
  group_by(name) %>% 
  summarise(matched_compound = toString(unique(matched_compound)), 
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
             str_detect(
               met_name_all,  
               "3-methyl pyruvic acid, Acetoacetic acid") ~ "Acetoacetic acid", 
             str_detect(
               met_name_all,  
               "Ascorbate, Glucuronolactone") ~ "Glucuronolactone", 
             str_detect(
               met_name_all, 
               "Phosphorylcholine, Vanylglycol") ~ "Vanylglycol",
             str_detect(
               met_name_all, 
               "Morphine, 2-Hydroxycarbamazepine") ~ 
               "2-Hydroxycarbamazepine",
             TRUE ~ met_name_first_only))


# Get pathway specific effect estimates for main exposure
pw_eff_est_annotated_allfts <- mwas_beta_coefs %>%
  modify(~.x %>% 
           left_join(key_pw_met2, ., by = c("name" = "feature")))  %>% 
  bind_rows(.id = "cohort")


# select only a single feature from each compound
most_altered_fts <- pw_eff_est_annotated_allfts %>% 
  group_by(name, matched_compound, exposure) %>% 
  summarise(estimate_mean = mean(estimate_beta)) %>%
  # mutate(rownum = row_number()) %>%
  ungroup() %>% 
  # select(name, matched_compound, exposure, rownum)
  group_by(matched_compound, exposure) %>%
  filter(abs(estimate_mean) == max(abs(estimate_mean))) %>%
  ungroup()


# Join back with pw_eff_est_annotated_allfts
pw_eff_est_annotated <- left_join(most_altered_fts, 
                                  pw_eff_est_annotated_allfts, 
                                  by = c("name", 
                                         "matched_compound",
                                         "exposure"))


# Pivot annotated dataset wider 
pw_eff_est_wide <- pw_eff_est_annotated %>%
  pivot_wider(names_from = "cohort",
              values_from = c("estimate_beta", "p_value", 
                              "q_value", 
                              "lcl_beta", "ucl_beta", 
                              "estimate_pip", ),
              id_cols = c("name", "exposure", 
                          "met_name",
                          "met_name_first_only",
                          "met_name_all", 
                          "matched_compound")) 


## Perform meta-analysis and calculate weighted mean effect estimate ------------
# set weights 
wgt_sol = sqrt(310); wgt_chs = sqrt(136)

# calculate meta-analysis p-values
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
pw_eff_est_wide_bile_acids <- pw_eff_est_wide %>% 
  mutate(sig = case_when(q_value_solar < 0.2 & q_value_chs < 0.2 ~ "q<0.2 both cohorts", 
                         p_value_solar < 0.05 & p_value_chs < 0.05 ~ 
                           "p<0.05 Both Cohorts",
                         p_value_solar < 0.05 ~ "p<0.05 (SOLAR only)",
                         p_value_chs < 0.05 ~ "p<0.05 (CHS only)",
                         TRUE ~ "Non. Sig"), 
         sig_solar = case_when(p_value_solar < 0.05 ~ "p<0.05",
                               p_meta < 0.05 ~ "Meta-analysis: p<0.05",
                               TRUE ~ "Not Sig."),
         sig_chs = case_when(p_value_chs < 0.05 ~ "p<0.05",
                             p_meta < 0.05 ~ "Meta-analysis: p<0.05",
                             TRUE ~ "Not Sig."),
         mean_effect = ((estimate_beta_solar*wgt_sol)+
                          (estimate_beta_chs*wgt_chs))/(wgt_sol+wgt_chs)) %>% 
  mutate(across(c(sig_solar, sig_chs), 
                ~fct_relevel(., 
                             "p<0.05",
                             "Meta-analysis: p<0.05",
                             "Not Sig.")))


# Create dataframe for single exposure of interest (mixture effect)
single_pfas_eff_est <- pw_eff_est_wide_bile_acids %>% 
  filter(exposure == "Mixture effect")


# 3) SOLAR Figure ------------------
## Reorder metabolites by mean effect  --------------
single_pfas_eff_est <- single_pfas_eff_est %>% 
  arrange(mean_effect) %>%
  mutate(matched_compound = fct_reorder(matched_compound,mean_effect), 
         met_name = fct_reorder(met_name,mean_effect), 
         met_name_first_only = fct_reorder(met_name_first_only,
                                           # mean_effect)) 
                                           estimate_beta_solar))



## Panel a: Coefficient plots -------------------------------
(solar_efest_fig <- ggplot(single_pfas_eff_est,
                           aes(x = met_name_first_only,
                               y = estimate_beta_solar,
                               color = sig_solar)) +
   geom_point(size = 1) +
   geom_errorbar(aes(ymin = lcl_beta_solar,
                     ymax = ucl_beta_solar),
                 width = 0) +
   geom_hline(yintercept = 0, linetype = 2) +
   scale_y_continuous(name = "Beta (95% CI)",
                      # limits = c(-1*max_min_effect, max_min_effect)
   ) +
   xlab("Metabolite") +
   ggtitle("A) PFAS Mixture") + 
   # scale_color_brewer(palette = "Set1") +
   # scale_color_manual(values = c("red", "black", "grey50")) +
   theme(axis.title.y =  element_blank(),
         legend.position = "bottom",
         # legend.justification='left',
         legend.direction='horizontal',
         legend.title = element_blank()) +
   # guides(color=guide_legend(nrow=2,byrow=TRUE, )) +
   coord_flip())


## Panel B & C: Heatmap of annotated features -------------------------

# Reorder metabolites by mixture effect estimates
pw_eff_est_w_2 <- pw_eff_est_wide %>% 
  filter(name %in% single_pfas_eff_est$name) %>%
  mutate(met_name_first_only = factor(
    met_name_first_only,
    levels = levels(single_pfas_eff_est$met_name_first_only), 
    ordered = TRUE), 
    estimate_sol_sig_only = if_else(p_value_solar < 0.5, 
                                    estimate_beta_solar, 
                                    NA_real_))


# Create heatmap effect estimates
(sol_heatmap_ee <- ggplot(pw_eff_est_w_2 %>% 
                            filter(exposure != "Mixture effect"), 
                          aes(x = exposure,
                              y = met_name_first_only,
                              fill = estimate_beta_solar)) + 
    geom_tile() +
    labs(title = "B) Effect estimates of individual PFAS",
         fill="Effect\nEstimate") +
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x =  element_text(angle = 90,
                                      hjust = 1, 
                                      vjust = 0.5)) + 
    scale_fill_gradient2(midpoint = 0,
                         low = "blue",
                         high = "red") )

# Create heatmap PIPS
(sol_heatmap_pips <- ggplot(pw_eff_est_w_2 %>% 
                              filter(exposure != "Mixture effect"), 
                            aes(x = exposure,
                                y = met_name_first_only,
                                fill = estimate_pip_solar)) + 
    geom_tile() +
    labs(title = "C) PIP's of individual PFAS",
         fill="Posterior\ninclusion\nprobability") +
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x =  element_text(angle = 90,
                                      hjust = 1, 
                                      vjust = 0.5)) + 
    scale_fill_gradient2(midpoint = 1/6,
                         mid = "white",
                         low = "blue",
                         high = "black") )

# Combine coef plot and heatmap
(solar_figure <- cowplot::plot_grid(solar_efest_fig,
                                    sol_heatmap_ee,
                                    sol_heatmap_pips, 
                                    nrow = 1,
                                    rel_widths = c(1.1, .9, .9),
                                    axis = "tb",
                                    align = "h"))


# Save Figure
ggsave(solar_figure,
       filename =  fs::path(
         dir_reports,
         "SOLAR PFAS Mixtures hyper_g bile_acid metabolites.jpg"),
       width = 14, height = 5, units = "in")


# 4) CHS Figure ------------------------
## Reorder metabolites by mean effect  --------------
single_pfas_eff_est <- single_pfas_eff_est %>%
  arrange(mean_effect) %>%
  mutate(matched_compound = fct_reorder(matched_compound,mean_effect),
         met_name = fct_reorder(met_name,mean_effect),
         met_name_first_only = fct_reorder(met_name_first_only,
                                           estimate_beta_chs))

## Panel a: Coefficient plots -------------------------------
(chs_efest_fig <- ggplot(single_pfas_eff_est,
                         aes(x = met_name_first_only,
                             y = estimate_beta_chs,
                             color = sig_chs)) +
   geom_point(size = 1) +
   geom_errorbar(aes(ymin = lcl_beta_chs,
                     ymax = ucl_beta_chs),
                 width = 0) +
   geom_hline(yintercept = 0, linetype = 2) +
   scale_y_continuous(name = "Beta (95% CI)",
                      # limits = c(-1*max_min_effect, max_min_effect)
   ) +
   xlab("Metabolite") +
   ggtitle("A) PFAS Mixture") + 
   # scale_color_brewer(palette = "Set1") +
   # scale_color_manual(values = c("red", "black", "grey50")) +
   theme(axis.title.y =  element_blank(),
         legend.position = "bottom",
         # legend.justification='left',
         legend.direction='horizontal',
         legend.title = element_blank()) +
   # guides(color=guide_legend(nrow=2,byrow=TRUE, )) +
   coord_flip())


## Panel B & C: Heatmap of annotated features -------------------------

# Reorder metabolites by mixture effect estimates
pw_eff_est_w_2 <- pw_eff_est_wide %>% 
  filter(name %in% single_pfas_eff_est$name) %>%
  mutate(
    met_name_first_only = factor(
      met_name_first_only,
      levels = levels(single_pfas_eff_est$met_name_first_only), 
      ordered = TRUE), 
    estimate_beta_chs_sig_only = if_else(p_value_chs < 0.5, 
                                         estimate_beta_chs, 
                                         NA_real_))


# Create heatmap of effect estimates
(chs_heatmap_ee <- ggplot(pw_eff_est_w_2 %>%
                            filter(exposure != "Mixture effect"), 
                          aes(x = exposure,
                              y = met_name_first_only,
                              fill = estimate_beta_chs)) + 
    geom_tile() +
    labs(title = "B) Contribution of individual PFAS",
         fill="Effect\nEstimate") +
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x =  element_text(angle = 90,
                                      hjust = 1, 
                                      vjust = 0.5)) + 
    scale_fill_gradient2(midpoint = 0,
                         low = "blue",
                         high = "red"))

# Create heatmap PIPS
(chs_heatmap_pips <- ggplot(pw_eff_est_w_2 %>% 
                              filter(exposure != "Mixture effect"), 
                            aes(x = exposure,
                                y = met_name_first_only,
                                fill = estimate_pip_chs)) + 
    geom_tile() +
    labs(title = "C) PIP's of individual PFAS",
         fill="Posterior\ninclusion\nprobability") +
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x =  element_text(angle = 90,
                                      hjust = 1, 
                                      vjust = 0.5)) + 
    scale_fill_gradient2(midpoint = 1/6,
                         mid = "white",
                         low = "blue",
                         high = "black") )


# Combine
(chs_figure <- cowplot::plot_grid(chs_efest_fig,
                                  chs_heatmap_ee,
                                  chs_heatmap_pips,
                                  nrow = 1,
                                  rel_widths = c(1.1, 0.9, 0.9),
                                  axis = "tb",
                                  align = "h"))

ggsave(chs_figure,
       filename =  fs::path(dir_reports,
                            "CHS PFAS Mixtures hyper_g bile_acid metabolites q020.jpg"),
       width = 14, height = 5, units = "in")


# 5) Combined Figure --------------------------------
# pw_eff_est_sig_only <- pw_eff_est_wide %>% 
#   filter(p_meta < 0.05, exposure == "Mixture effect")
# 
# 
# pw_eff_est_wide_bile_acid <- pw_eff_est_wide
# 
# (tyr_met_comp <- ggplot(pw_eff_est_wide_bile_acid, aes(x = estimate_solar, y = estimate_beta_chs, 
#                                                       color = sig)) + 
#     geom_point() + 
#     geom_hline(yintercept = 0) + 
#     geom_vline(xintercept = 0) + 
#     geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey30") +
#     facet_wrap(~exposure) + 
#     ylab("Beta (CHS)") + 
#     xlab("Beta (SOLAR)") + 
#     ggtitle("bile_acid Metabolites"))
# 
# 
# 
# ggsave(tyr_met_comp, 
#        filename = fs::path(dir_reports, 
#                            "SOL CHS Scatter plot of bile_acid metabolites.jpg"), 
#        width = 8, height = 7, units = "in")
# 
# 
# 
# 
# 
# 
# 
# 
# # Create figure of all metabolites 
# pw_eff_est_wide
# 
# 
# (temp_fig <- ggplot(pw_eff_est_wide, aes(x = estimate_solar, y = estimate_beta_chs, 
#                                          color = sig)) + 
#     geom_point() + 
#     geom_hline(yintercept = 0) + 
#     geom_vline(xintercept = 0) + 
#     geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey30") +
#     facet_wrap(~exposure) + 
#     ylab("Beta (CHS)") + 
#     xlab("Beta (SOLAR)") + 
#     ggtitle("bile_acid Metabolites"))
