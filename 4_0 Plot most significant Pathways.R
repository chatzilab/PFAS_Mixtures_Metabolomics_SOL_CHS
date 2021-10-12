# Get most significant pathway from each PFAS
library(ggExtra)

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
         mean_num_sig = (hits_sig_solar + hits_sig_chs)/2)



# Get most significant for each PFAS
mum_most_sig_pws <- mum_pw_wide %>% 
  filter(fet_meta < 0.05, 
         sig == "Sig. Both Cohorts")

#   group_by(pfas) %>% 
#   filter(neg_logp_meta == max(neg_logp_meta)) %>%
#   ungroup()



# Read in MWAS Beta Coefficients --------------------
mwas_beta_coefs <- read_rds(
  fs::path(dir_temp, 
           "PFAS", 
           "SOL CHS all MWAS results long.rds")) %>% 
  modify(~.x %>% mutate(exposure2 = exposure %>% 
                          str_remove(., "lg2_") %>% 
                          rename_pfas(include_asterisk = TRUE)))


# Read in mzrt key
mzrt_key <- read_rds(fs::path(dir_temp, 
                              exposure_type, 
                              "mummichog_pw_ec_feature_key.rds")) 



# Create Functions --------------------------------------------------
# Set vars to test fxn
n = 17
data <- mzrt_key
pfas_name <- "PFDS*"
cohort_name <- "solar"
pathway_name <- mum_most_sig_pws$path[n]
# rm(n, data, pfas_name, cohort_name, pathway_name)

# Function for creating 
plot_pathway <- function(data,pfas_name, cohort_name, pathway_name){ 
  # Filter metabolites in the specified pathway
  key_pathway_met <- data %>% filter(pathway == pathway_name)
  # Get MWAS Results
  pathway_mwas_results <- mwas_beta_coefs %>%
    modify(~.x %>% 
             filter(name %in% unique(key_pathway_met$name)) )  %>% 
    bind_rows(.id = "cohort") 
  
  # Select PFAS
  pathway_mwas_pfas_cohort_1 <- pathway_mwas_results %>% 
    filter(exposure2 == pfas_name)
  
  # Get max and min effect
  max_ef <- max(pathway_mwas_pfas_cohort_1$conf_high)
  min_ef <- max(abs(pathway_mwas_pfas_cohort_1$conf_low))
  max_min_effect = max(c(max_ef, min_ef))        
  
  # Filter cohort 
  pathway_mwas_pfas_cohort <- pathway_mwas_pfas_cohort_1 %>% 
    filter(cohort == cohort_name) %>% 
    arrange(estimate) %>%
    mutate(name = fct_reorder(name,estimate)) 
  #plot
  
  # Create effect_est_fig Figure
  plotout <- ggplot(pathway_mwas_pfas_cohort, 
                    aes(x = name, 
                        y = estimate, 
                        color = significancefdr)) +
    geom_point(alpha = .7, size = .9) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(name = "Beta (95% CI)", 
                       limits = c(-1*max_min_effect, max_min_effect)) +
    theme(axis.text.x = element_blank(), 
          legend.position = "none") +
    xlab("Feature") +
    geom_errorbar(aes(ymin = conf_low, 
                      ymax = conf_high),
                  alpha = .5,
                  width = 0) #+ 
  # ggtitle(str_c(cohort_name, ", ", pfas_name, ", ", key_pathway_met$pathway[1])) 
  
  
  # Pivot wider to create 
  pathway_mwas_pfas_cohort_w <- pathway_mwas_pfas_cohort_1 %>% 
    pivot_wider(., names_from = cohort, 
                values_from = c("estimate", "p_value", "q_value"), 
                id_cols = c("name")) %>% 
    mutate(sig = case_when(p_value_solar < 0.05 & p_value_chs < 0.05 ~ "Sig. Both Cohorts", 
                           p_value_solar < 0.05 ~ "Sig. SOLAR", 
                           p_value_chs < 0.05 ~ "Sig. CHS", 
                           TRUE ~ "Non. Sig"))
  
  
  # Plot Scatter plots 
  p3 <- ggplot(pathway_mwas_pfas_cohort_w,
               aes(x = estimate_solar,
                   y = estimate_chs,
                   color = sig, 
                   shape = sig)) +
    geom_hline(yintercept = 0, linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    stat_smooth(method = "lm", color = "grey50", aes(group = "1")) +
    scale_shape_manual(values = c(20, 19, 20, 20)) +
    geom_point(size = 3) +
    theme(legend.title=element_blank()) +
    xlab("Beta SOLAR") + 
    ylab("Beta CHS")
  
  
  return(list(ggMarginal(plotout, type = "boxplot", margins = "y"), 
              p3))
}


# Plot both cohorts -------------------------------------------------------------
pfas_name <- "PFHxS"
plot_sol_chs <- function(pfas_name, return = FALSE){ 
  
  pfas_pathways <- which(mum_most_sig_pws$pfas == pfas_name)
  
  
  for(n in pfas_pathways){
    p1 <- plot_pathway(mzrt_key, pfas_name = pfas_name, "solar", "Aspartate and asparagine metabolism")[[1]] #mum_most_sig_pws$path[n])[[1]]
    p2 <- plot_pathway(mzrt_key, pfas_name = pfas_name, "chs", "Aspartate and asparagine metabolism")[[1]] #mum_most_sig_pws$path[n])[[1]]
    
    p3 <- plot_pathway(mzrt_key, pfas_name = pfas_name, "solar", "Aspartate and asparagine metabolism")[[2]] #mum_most_sig_pws$path[n])[[2]]
    # Row 2
    row2 <- plot_grid(NULL, NULL,
                      p1, p2, 
                      labels = c('A. SOLAR', 
                                 'B. CHS', 
                                 "",
                                 ""),
                      nrow = 2,label_size = 12,
                      rel_heights = c(.1, 1), 
                      label_x = 0,
                      hjust = 0) 
    
    ## Title
    title <- ggdraw() + 
      draw_label(
        str_c("Associations of ", 
              pfas_name,
              " exposure with ",
              mum_most_sig_pws$path_2[n]),
        fontface = 'bold',
        x = 0,size = 14,
        hjust = 0) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7))
    
    
    # Figure 
    outfig <- plot_grid(title, 
                        row2,
                        NULL,
                        p3,
                        labels = c("", "", "C. Effect Estimates SOLAR vs. CHS", ""),
                        nrow = 4, label_size = 12,
                        rel_heights = c(.1, 1,.1, 1), 
                        label_x = 0,
                        hjust = 0)
    
    
    ggsave(outfig,filename =  
             fs::path(dir_reports, 
                      "PFAS Pathway Metabolite Figures", 
                      str_c(str_replace_all(pfas_name, "[[:punct:]]", " "), 
                            " with ", 
                            str_replace_all(mum_most_sig_pws$path[n], "[[:punct:]]", " "),
                            ".jpg")), 
           height = 6.5, width = 6.5)
    
    
    if(return == TRUE){return(outfig)}
  }
}


# Create Plots: -------------------------------------
# PFOS
for(p in unique(mum_most_sig_pws$pfas)[4:11]){
  plot_sol_chs(p)
}

plot_sol_chs("PFDS*")


