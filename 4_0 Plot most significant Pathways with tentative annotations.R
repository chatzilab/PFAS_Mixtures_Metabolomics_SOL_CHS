# Get most significant pathway from each PFAS
library(ggExtra)
library(tidylog)
# 2) Plot Mummichog Pathway Results
library(colorspace)
library(janitor)
exposure_type = "PFAS"

mum_pw_wide <- read_rds(
  fs::path(dir_temp, 
           exposure_type, 
           "mum_pathway_results", 
           "SOL CHS PFAS Mummichog wide sig PW.RDS")) %>% 
  clean_names() %>% 
  mutate(q_meta = p.adjust(fet_meta), 
         sig_fdr = if_else(q_meta < 0.2, "Sig. FDR < 0.2", "Not. Sig"), 
         mean_num_sig = (hits_sig_solar +hits_sig_chs)/2)



# Get most significant for each PFAS
mum_most_sig_pws <- mum_pw_wide %>% 
  filter(fet_meta < 0.05, 
         sig == "Sig. Both Cohorts") %>% 
  mutate(mean_hits = (hits_sig_solar + hits_sig_chs)/2) %>% 
  arrange(-1*mean_hits)





# Read in MWAS Beta Coefficients --------------------
mwas_beta_coefs <- read_rds(
  fs::path(dir_temp, 
           "PFAS", 
           "SOL CHS all MWAS results long.rds")) %>% 
  modify(~.x %>% mutate(exposure2 = exposure %>% 
                          str_remove(., "lg2_") %>% 
                          rename_pfas(include_asterisk = TRUE)))


# Read in mzrt key
mzrt_key <- read_rds(fs::path(dir_data_local,
                              "Supporting Files",
                              "mummichog_pw_ec_feature_key_with_cpd_names.rds"))

View(mzrt_key)
# Set Pathway
n = 26 # Tyrosine metabolism

pfas_name <- "PFHxS"
cohort_name <- "solar"
(pathway_name <- mum_most_sig_pws$path[n])
pathway_name = "Tyrosine metabolism"


# FUNCTION for creating figure---------------------------------

# make_pathway_fig <- function(pathway_name, pfas_name){
  
  # Get mzrt key for specific pathway ----------------------
  key_pathway_met <- mzrt_key %>% 
    filter(pathway == pathway_name) %>% 
    select(mode, pathway, met_name, everything())
  
  # Group by mz/rt, create unique compound names for each feature
  ## Important: name cpds with multiple names -----------------------
  key_pathway_met2 <- key_pathway_met %>% 
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
               TRUE ~ met_name_first_only))
  
  
  
  # Get pathway specific effect estimates for main exposure
  pathway_eff_est_annotated_allfts <- mwas_beta_coefs %>%
    modify(~.x %>% 
             left_join(key_pathway_met2, ., by = "name"))  %>% 
    bind_rows(.id = "cohort") #%>% 
  # filter(exposure2 == pfas_name)
  
  
  # get single feature from each compound
  most_altered_fts <- pathway_eff_est_annotated_allfts %>% 
    group_by(name, matched_compound, exposure, exposure2) %>% 
    summarise(estimate_mean = mean(estimate)) %>%
    # mutate(rownum = row_number()) %>%
    ungroup() %>% 
    # select(name, matched_compound, exposure, exposure2, rownum)
    group_by(matched_compound, exposure, exposure2) %>%
    filter(abs(estimate_mean) == max(abs(estimate_mean))) %>%
    ungroup()
  # 
  
  # Join back with pathway_eff_est_annotated_allfts
  pathway_eff_est_annotated <- left_join(most_altered_fts, 
                                         pathway_eff_est_annotated_allfts, 
                                         by = c("name", 
                                                "matched_compound",
                                                "exposure", 
                                                "exposure2"))
  
  
  # Pivot annotated dataset wider 
  pathway_eff_est_wide <- pathway_eff_est_annotated %>%
    pivot_wider(., 
                names_from = cohort,
                values_from = c("estimate", "p_value", 
                                "q_value", "conf_low", "conf_high"),
                id_cols = c("name", "exposure", "exposure2",
                            "met_name",
                            "met_name_first_only",
                            "met_name_all", 
                            "matched_compound")) 
  
  
  # Subset specific PFAS of interest
  single_pfas_effect_est <- pathway_eff_est_wide %>% 
    filter(exposure2 == pfas_name)
  
  
  # Perform Meta analysis
  wgt_sol = sqrt(310); wgt_chs = sqrt(136)
  
  # Perform Meta analysis ---------------------------
  single_pfas_effect_est <- single_pfas_effect_est %>% 
    mutate(sig = case_when(p_value_solar < 0.05 & p_value_chs < 0.05 ~ 
                             "Sig. Both Cohorts",
                           p_value_solar < 0.05 ~ "Sig. SOLAR",
                           p_value_chs < 0.05 ~ "Sig. CHS",
                           TRUE ~ "Non. Sig"), 
           mean_effect = ((estimate_solar*wgt_sol)+
                            (estimate_chs*wgt_chs))/(wgt_sol+wgt_chs)) %>% 
    rowwise() %>% 
    mutate(p_meta = metap::sumz(p = c_across(p_value_solar:p_value_chs), 
                                weights = c(wgt_sol, wgt_chs))$p, 
           neg_logp_meta = -log(p_meta), 
           sig_meta = if_else(p_meta < 0.05, "p<0.05", "p>0.05")) %>% 
    ungroup() %>% 
    arrange(mean_effect) %>%
    mutate(matched_compound = fct_reorder(matched_compound,mean_effect), 
           met_name = fct_reorder(met_name,mean_effect), 
           met_name_first_only = fct_reorder(met_name_first_only,
                                             mean_effect)) 
  
  # write_csv(data.frame(metabolite = unique(pathway_eff_est_wide$matched_compound)), 
  #           here::here("Aspartate and asparagine metabolism metabolites.csv"))
  
  # Create Figures -------------------------------
  ## SOLAR ------------------
  (solar_efest_fig <- ggplot(single_pfas_effect_est,
                             aes(x = met_name_first_only,
                                 y = estimate_solar,
                                 color = sig_meta)) +
     geom_point(size = .9) +
     geom_errorbar(aes(ymin = conf_low_solar,
                       ymax = conf_high_solar),
                   width = 0) +
     geom_hline(yintercept = 0, linetype = 2) +
     scale_y_continuous(name = "Beta (95% CI)",
                        # limits = c(-1*max_min_effect, max_min_effect)
     ) +
     xlab("Metabolite") +
     theme(axis.title.y =  element_blank(),
           legend.position = "none") +
     coord_flip()
  )
  
  ## CHS ----------------------------
  (chs_efest_fig <- ggplot(single_pfas_effect_est,
                           aes(x = met_name_first_only,
                               y = estimate_chs,
                               color = sig_meta)) +
     geom_point(size = .9) +
     geom_errorbar(aes(ymin =  conf_low_chs,
                       ymax = conf_high_chs),
                   width = 0) +
     geom_hline(yintercept = 0, linetype = 2) +
     scale_y_continuous(name = "Beta (95% CI)",
                        # limits = c(-1*max_min_effect, max_min_effect)
     ) +
     guides(color = guide_legend(title = "Meta p-value")) +
     theme(axis.title.y =  element_blank(),
           axis.text.y = element_blank(),
           legend.position = "right") +
     coord_flip()) #+
  
  
  (row2 <- plot_grid(NULL, NULL,
                     solar_efest_fig, chs_efest_fig,
                     labels = c('A. SOLAR',
                                'B. CHS'),
                     nrow = 2,
                     label_size = 12,
                     rel_widths = c(1, .75),
                     rel_heights = c(.02, 1),
                     label_x = 0,
                     hjust = 0))
  
  
  ## Plot both cohorts ---------------------------------------------
  
  ## Title
  title <- ggdraw() +
    draw_label(
      str_c("Associations of ",
            pfas_name,
            " exposure with ",
            pathway_name),
      fontface = 'bold',
      x = 0,size = 14,
      hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  
  
  # Figure
  (outfig <- plot_grid(title,
                       row2,
                       labels = c("", ""),
                       nrow = 2, 
                       label_size = 12,
                       rel_heights = c(.1, 1),
                       label_x = 0,
                       hjust = 0))
  
  return(outfig)
  
# }

outfig <- make_pathway_fig("Tyrosine metabolism", "PFHxS")



ggsave(outfig,
       filename =
         fs::path(dir_reports,
                  "PFAS Pathway Metabolite Figures",
                  str_c(str_replace_all(pfas_name, "[[:punct:]]", " "),
                        " with ",
                        str_replace_all(pathway_name, "[[:punct:]]", " "),
                        ".jpg")),
       height = 13, width = 9)



# Function to create all plots (not complete) -------------------------------------
# # PFOS
# for(p in unique(mum_most_sig_pws$pfas)[4:11]){
#   plot_sol_chs(p)
#   
# }
# plot_sol_chs("PFDS*")
# 
# 
