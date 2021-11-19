# Identify main pathways of interest: Tyrosine and PFHxS
# library(MetaboAnalystR)
library(tidylog)
library(RColorBrewer)
library(dendextend)
library(gplots)


# Read in MWAS Beta Coefficients --------------------
# mwas_results_long <- read_rds(
#   fs::path(dir_results, exposure_type,
#            "SOL CHS all MWAS results long.rds")) 
# 
# 
# 
# 
# #Read annotated data 
# ann_fts <- read_rds(fs::path(dir_data,"sol_chs_batch_cor_scaled_tentative_annotated_fts.rds"))
# 
# 
# 
# 
# 
# # Scale estimates 
# mwas_results_long <- mwas_results_long %>% 
#   modify(~.x %>% 
#            group_by(exposure) %>% 
#            mutate(estimate_scaled = scale(estimate, center = FALSE))) 
# 
# 
# # Pivot wider, change non-sig estimates to 0
# mwas_results_wide <- mwas_results_long %>% 
#   modify(~mutate(.x, estimate_scaled = if_else(p_value < 0.05, 
#                                                estimate_scaled, 
#                                                0)) %>%
#            # filter(str_detect(exposure, "lg2"), 
#            #        str_detect(exposure, "pfpes", negate = TRUE)) %>%
#            select(exposure, name, estimate_scaled) %>% 
#            pivot_wider(id_cols = name, 
#                        names_from = exposure, 
#                        values_from = estimate_scaled))
# 
# 
# 
# # Read in Dougs annotations ------------------------------------------------
# # annotations <- read_rds(
# #   fs::path(dir_data, 
# #            "4_Common_Metabolites_Annotation", 
# #            "Common_Metabolites_SOLAR_CHS_V1.RDS")) %>% 
# #   bind_rows(.id = "mode")
# # 
# # 
# # tem <- read_csv(here::here("Temporary results", 
# #                            "Pathway Results Annotated Metabolites", 
# #                            "All Annotated Metabolites Hand Curated.csv"))
# # 
# # tyrosine_pw <- readxl::read_xlsx(here::here("Temporary results", 
# #                                             "Pathway Results Annotated Metabolites", 
# #                                             "Pathway Metabolites",
# #                                             "Tyrosine Pathway.xlsx"), col_names = FALSE) %>% 
# #   clean_names()
# 
# 
# 
# 
# # Tyrosine pathway
# cpdnms <- c("2,5-Dihydroxybenzoate",
#             "3-(4-Hydroxyphenyl)pyruvate",
#             "3,4-Dihydroxy-L-phenylalanine",
#             "3-Methoxy-4-hydroxymandelate",
#             "3-Methoxytyramine",
#             "4-Hydroxyphenylacetate",
#             "Acetoacetate",
#             "Dopamine",
#             "Fumarate",
#             "Homogentisate",
#             "Homovanillate",
#             "L-Adrenaline",
#             "L-Normetanephrine",
#             "L-Tyrosine",
#             "Pyruvate",
#             "Thyroxine",
#             "Tyramine")
# 
# 
# ann_fts$solar$kegg
# 
# ## Look at 
# tyr_fts <- ann_fts$solar %>% 
#   # filter(refmet_name %in% cpdnms) %>% 
#   ungroup()
# 
# 
# 
# # Select only study samples
# tyr_fts2 <- tyr_fts %>% 
#   select(refmet_name, contains("sol")) 
# 
# # Transpose DF
# con_df <- tyr_fts2 %>% 
#   mutate(across(everything(), as.character)) %>%
#   pivot_longer(names_to = "id", values_to = "conc", cols = c(2:ncol(tyr_fts2)) ) %>% 
#   group_by(refmet_name, id) %>%
#   mutate(row = row_number(),
#          refmet_name = str_c(refmet_name, row, sep = "__")) %>%
#   select(-row) %>%
#   ungroup() %>%
#   mutate(conc = as.numeric(conc)) #%>% 
# # pivot_wider(names_from = refmet_name, values_from = conc) 
# 
# 
# 
# library(tidylog)
# 
# 
# # Merge with outcomes 
# full_data <- left_join(exposure_outcome$solar, con_df) %>% 
#   mutate(pfhxs_dicot = case_when(pfhxs > quantile(pfhxs, .75) ~ "High", 
#                                  pfhxs < quantile(pfhxs, .25) ~ "Low") %>% 
#            fct_relevel("Low", "High"), 
#          refmetnamelwr = tolower(refmet_name)) 
# 
# 
# ## PLOT
# ggplot(full_data %>% 
#          filter(!is.na(pfhxs_dicot), 
#                 str_detect(refmetnamelwr, "tyrosine") | 
#                   str_detect(refmetnamelwr, "dopamine")  | 
#                   str_detect(refmetnamelwr, "homovanillate") | 
#                   str_detect(refmetnamelwr, "metanephrin") | 
#                   str_detect(refmetnamelwr, "ate") | 
#                   str_detect(refmetnamelwr, "metanephrin") | 
#                   str_detect(refmetnamelwr, "metanephrin") | 
#                   str_detect(refmetnamelwr, "metanephrin") | 
#                   str_detect(refmetnamelwr, "metanephrin")  
#          ), 
#        aes(x = pfhxs_dicot, y = conc)) +
#   geom_boxplot() + 
#   facet_wrap(~refmet_name)
# 
# 
# cpdnms <- c("2,5-Dihydroxybenzoate",
#             "3-(4-Hydroxyphenyl)pyruvate",
#             "3,4-Dihydroxy-L-phenylalanine",
#             "3-Methoxy-4-hydroxymandelate",
#             "3-Methoxytyramine",
#             "4-Hydroxyphenylacetate",
#             "Acetoacetate",
#             "Dopamine",
#             "Fumarate",
#             "Homogentisate",
#             "Homovanillate",
#             "L-Adrenaline",
#             "L-Normetanephrine",
#             "L-Tyrosine",
#             "Pyruvate",
#             "Thyroxine",
#             "Tyramine")
# 




















# New, 9/25 ---------------------------------------------------------------

# Read in MWAS Beta Coefficients --------------------
mwas_beta_coefs <- read_rds(
  fs::path(dir_results, 
           "PFAS", 
           "SOL CHS all MWAS results long.rds"))


# Read in mzrt key
mzrt_key <- read_rds(fs::path(dir_results, exposure_type,  "mummichog_pw_ec_feature_key.rds")) 






# Select most impacted pathwaysnv_snapshot_fixup_renv(records)


plot_pathway <- function(data,pfas_name, cohort_name, pathway_name){ 
  
  key_pathway_met <- data %>% filter(str_detect(pathway, pathway_name))
  
  # Get MWAS Results
  pathway_mwas_results <- mwas_beta_coefs %>%
    modify(~.x %>% 
             filter(name %in% unique(key_pathway_met$name)) )%>% 
    bind_rows(.id = "cohort") 
  
  
  pathway_mwas_pfas_cohort <- pathway_mwas_results %>% 
    filter(str_detect(exposure, pfas_name), 
           cohort == cohort_name) %>% 
    arrange(estimate) %>%
    mutate(name = fct_reorder(name,estimate)) 
  #plot
  
  
  plotout <- ggplot(pathway_mwas_pfas_cohort, 
                      aes(x = name, y = estimate, color = significancefdr)) +
    geom_point(alpha = .7, size = .9) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme(axis.text.x = element_blank()) +
    geom_errorbar(aes(ymin = conf_low, 
                      ymax = conf_high), width = 0) + 
    ggtitle(str_c(cohort_name, ", ", pfas_name, ", ", key_pathway_met$pathway[1]))
  
  return(plotout)
  
}



## PFDA, SOLAR, Aspartate Metabolites
plot_pathway(mzrt_key,"pfda", "solar", "Aspartate")


# Tyrosine metabolism
plot_pathway(mzrt_key,"pfos", "solar", "Tyrosine")


# Arginine 

plot_pathway(mzrt_key,"pfhxs", "solar", "Arginine")
