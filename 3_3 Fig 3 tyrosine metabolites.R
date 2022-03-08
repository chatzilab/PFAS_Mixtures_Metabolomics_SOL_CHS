# Make coefficient plots based on dougs comments

# Read data
annotated_sig_ecs_ee <- read_rds(
  here::here(dir_results_mixtures,
             "Sig annotated metabolite effect estimates sol chs.RDS") ) %>% 
  select(-contains("blank")) %>% 
  tidylog::filter(str_detect(pathway, "Tyrosine"))


# Change metabolite names to those associated with Tyrosine metabolism --------------
annotated_sig_ecs_ee <- annotated_sig_ecs_ee %>% 
  tidylog::mutate(met_name_tyr_pw = case_when(
    met_name == "Ascorbate; Glucuronolactone; D-glucurono-3,6-lactone" ~ 
      "Ascorbate", 
    met_name == "Pyridoxine; Norepinephrine" ~ "Norepinephrine",
    met_name == "Dopaquinone; Leucodopachrome" ~ "Dopaquinone",
    met_name == "Formylanthranilic acid; Noradrenochrome" ~ 
      "Noradrenochrome",
    met_name == "Vanylglycol; Phosphorylcholine" ~ "Vanylglycol",
    met_name == "Adrenochrome; Hippuric acid" ~ "Hippuric acid", 
    met_name == "Homovanillic acid; 3-Methoxy-4-hydroxyphenylglycolaldehyde" ~
      "Homovanillic acid", 
    met_name == "3-Methoxytyramine; Epinine" ~ "3-Methoxytyramine", 
    met_name == "5-Hydroxytryptophol; 1,2-dehydrosalsolinol" ~ 
      "1,2-dehydrosalsolinol", 
    met_name == "Phenylacetylglutamine; Acetyl-N-formyl-5-methoxykynurenamine" ~ 
      "Phenylacetylglutamine", 
    met_name == "L-Glutamic acid; D-Glutamic acid; L-4-Hydroxyglutamate semialdehyde" ~
      "L-Glutamic acid",
    met_name == "L-Glutamic acid; D-Glutamic acid; DL-Glutamate; L-4-Hydroxyglutamate semialdehyde" ~ 
      "L-Glutamic acid",
    met_name == "sulfuric acid [4-(2-aminoethyl)phenyl] ester" ~ 
      "Tyramine-O-sulfate", 
    met_name == "4-Hydroxyphenylacetaldehyde; Phenylacetic acid" ~ 
      "4-Hydroxyphenylacetaldehyde", 
    met_name == "3-methyl pyruvic acid; Acetoacetic acid; Succinic acid semialdehyde; 2-Methyl-3-oxopropanoic acid; (S)-Methylmalonic acid semialdehyde" ~ "Acetoacetic acid",
    met_name == "Pyruvic acid; Malonic semialdehyde" ~ "Pyruvic acid",
    TRUE ~ met_name))

# Add tyrosine sub pathways
annotated_sig_ecs_ee <- annotated_sig_ecs_ee %>% 
  tidylog::mutate(
    tyr_subpath = case_when(
      met_name_tyr_pw %in% c("3,4-Dihydroxyphenylglycol","Metanephrine",
                             "Noradrenochrome", "Vanylglycol","Ascorbate",
                             "Adrenochrome","Norepinephrine",
                             "Norepinephrine sulfate", "L-Dopa", 
                             "1,2-dehydrosalsolinol","3-Methoxytyramine",
                             "3-O-methyldopa","Homovanillic acid",
                             "Homovanillin") ~
        "Catecholamine biosynthesis\nand degredation", 
      
      met_name_tyr_pw %in%  c("Phenylacetaldehyde",
                              "Phenylacetylglutamine", 
                              "Hippuric acid") ~ 
        "Phenylalanine metabolism",
      
      met_name_tyr_pw %in%  c("4-Hydroxyphenylacetaldehyde", 
                              "L-Glutamic acid", 
                              "Tyramine-O-sulfate", 
                              "Acetoacetic acid",
                              "Pyruvic acid") ~
        "Tyrosine metabolism and\ndegredation", 
      
      met_name_tyr_pw == "Dopaquinone" ~ "Melanin biosynthesis", 
      met_name_tyr_pw == "Thyroxine" ~ "Thyroid Hormone Biosynthesis"))


# Modify table based on Dougs feedback (comments are directly from Doug)
sig_ecs_reduced <- annotated_sig_ecs_ee %>% 
  tidylog::filter(
    #For query_mass == 167.0701851, this is homovanillin, not vanyglycl-
    # multiple adducts at the same retention time corresponding to this compound (~125 sec)
    !(met_name == "Vanylglycol; Phosphorylcholine" & query_mass == 167.0701851), # Removes 1 row
    
    # This is not homovanilian, since it does not match the other 3 adduct retention times:
    !(met_name == "Homovanillin" & retention_time  == "72.3"), # Removes 1 row
    
    #This is a very weird adduct, I would say this is not correct:
    matched_form != "M-HCOOH+H[1+]" # Removes 1 row
  ) 


# Get Mass difference for the tyrosine metabolites 
sig_ecs_reduced <- sig_ecs_reduced %>% 
  group_by(met_name, name) %>% 
  tidylog::mutate(
    mass_diff_ls = str_split(mass_diff, "; "), 
    mass_diff_tyr_met = case_when(
      met_name == "Formylanthranilic acid; Noradrenochrome" ~ mass_diff_ls[[1]][2],
      met_name == "Vanylglycol; Phosphorylcholine" ~ mass_diff_ls[[1]][2],
      met_name == "Ascorbate; Glucuronolactone; D-glucurono-3,6-lactone" ~ 
        mass_diff_ls[[1]][1], 
      TRUE ~ mass_diff))  %>% 
  ungroup() %>% 
  select(met_name_tyr_pw,  tyr_subpath, mass_diff_tyr_met,  everything(), 
         -mass_diff_ls,  -casnum, -HRE_standard)   


# Reorder tyrosine sub path by number of cpds within path
sig_ecs_reduced <- sig_ecs_reduced %>% 
  mutate(tyr_subpath = fct_infreq(tyr_subpath), 
         borderline_sig_sol = if_else(lcl_beta_sol_mixture == 0 | 
                                        ucl_beta_sol_mixture == 0, 
                                      "borderline", 
                                      "sig"), 
         borderline_sig_chs = if_else(lcl_beta_chs_mixture == 0 | 
                                        ucl_beta_chs_mixture == 0, 
                                      "borderline", 
                                      "sig"))

# Save Data -------------------------------
# write_rds(sig_ecs_reduced,
#           here::here(dir_results_mixtures,
#                     "Sig tyrosine metabolite effect estimates sol chs.RDS") )



# Make Plots ----------------------------------------------------------
# Subset SOL (filter sig features)
sol <- sig_ecs_reduced %>% 
  select(met_name_tyr_pw:mass_diff,empirical_compound, 
         contains("sol")) %>% 
  tidylog::filter(ucl_beta_sol_mixture <= 0 | 
                    lcl_beta_sol_mixture >= 0) %>% 
  group_by(met_name_tyr_pw) %>% 
  filter(estimate_beta_sol_mixture == max(estimate_beta_sol_mixture)) %>% 
  ungroup() %>% 
  mutate(met_name_tyr_pw = fct_reorder(met_name_tyr_pw, 
                                       estimate_beta_sol_mixture), 
         borderline_sig = if_else(lcl_beta_sol_mixture == 0 | ucl_beta_sol_mixture == 0, 
                                  "borderline", 
                                  "sig"))


# Subset CHS (filter sig features)
chs <- sig_ecs_reduced %>% 
  select(met_name_tyr_pw:mass_diff,empirical_compound, 
         contains("chs")) %>% 
  tidylog::filter(p_value_chs_mixture < 0.05 |  
                    ucl_beta_chs_mixture <= 0 | 
                    lcl_beta_chs_mixture >= 0) %>% 
  group_by(met_name_tyr_pw) %>% 
  tidylog::filter(estimate_beta_chs_mixture == max(estimate_beta_chs_mixture)) %>% 
  ungroup() %>% 
  mutate(met_name_tyr_pw = fct_reorder(met_name_tyr_pw, 
                                       estimate_beta_chs_mixture))

# Get data on overlap
sol <- sol %>% 
  mutate(sig_both = met_name_tyr_pw %in% chs$met_name_tyr_pw, 
         met_name_tyr_pw_with_overlap = if_else(sig_both, 
                                                str_c(met_name_tyr_pw, "*"), 
                                                as.character(met_name_tyr_pw)) %>% 
           fct_reorder(estimate_beta_sol_mixture))

chs <- chs %>% 
  mutate(sig_both = met_name_tyr_pw %in% sol$met_name_tyr_pw, 
         met_name_tyr_pw_with_overlap = if_else(sig_both, 
                                                str_c(met_name_tyr_pw, "*"), 
                                                as.character(met_name_tyr_pw)) %>% 
           fct_reorder(estimate_beta_chs_mixture))

# Plot solar------------------------------
(sol_met_plot <- ggplot(sol, aes(x = met_name_tyr_pw_with_overlap, 
                                 y = estimate_beta_sol_mixture, 
                                 color = borderline_sig_sol)) + 
   geom_errorbar(aes(ymin = lcl_beta_sol_mixture , 
                     ymax = ucl_beta_sol_mixture), 
                 width = 0) + 
   geom_point() +
   geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
   facet_grid(tyr_subpath ~ .,  scales = "free_y", space = "free_y") +
   scale_color_manual(values = c("grey60", "black")) +
   coord_flip(clip = "off") + 
   scale_y_continuous(limits = c(-1.75, 3)) +
   # ylab("PFAS Exposure Effect Estimate ψ (95% BCI)") +
   theme(axis.title.y = element_blank(), 
         axis.title.x = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         strip.text.y = element_text(angle = 0, hjust = 0))) 

# Plot ------------------------------
(chs_metplot <- ggplot(chs, aes(x = met_name_tyr_pw_with_overlap, 
                                y = estimate_beta_chs_mixture, 
                                color = borderline_sig_chs)) + 
   geom_errorbar(aes(ymin = lcl_beta_chs_mixture , 
                     ymax = ucl_beta_chs_mixture), 
                 width = 0) + 
   geom_point() +  #aes(color = tyr_subpath)
   geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
   facet_grid(tyr_subpath~"1",  scales = "free_y", space = "free_y") +
   scale_color_manual(values = c("grey60", "black"), 
                      name = "Tyrosine Sub pathway") +
   coord_flip() + 
   scale_y_continuous(limits = c(-1.75, 3)) +
   ylab("PFAS Mixture Effect ψ (95% CI)") +
   theme(axis.title.y = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         strip.text.y = element_text(angle = 0, hjust = 0)))



# Combine Plots
figure_3 <- plot_grid(NULL, sol_met_plot, 
                      NULL, chs_metplot, 
                      ncol = 1, align = "v", 
                      rel_heights = c(.05, 1,.05, .5),
                      labels = c("A) SOLAR","",
                                 "B) CHS", ""))


ggsave(figure_3, 
       filename = fs::path(dir_reports, 
                           "Figure 4 associations of PFAS and tyr metabolites.jpg"), 
       width = 8, height = 9)
