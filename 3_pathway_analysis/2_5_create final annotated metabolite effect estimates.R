# Get most significant pathway from each PFAS
library(tidyverse)
library(ggExtra)

# 2) Plot Mummichog Pathway Results
# library(colorspace)
library(janitor)

# Set vars
# cohort_name <- "solar"
source(here::here("0_project_setup", "!functions.R"))
source(here::here("0_project_setup", "!directories.R"))


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
         mean_num_sig = (hits_sig_sol + hits_sig_chs)/2)

# Get pathways which are significant in both cohorts
sig_pw_both_chrt <- mum_pw_wide %>% 
  filter(#sig == "Sig. Both Cohorts",
    pval_meta < 0.05,
    mixture == "all_pfas")


## MWAS Data ------------------------------------
mwas_results <- read_rds(
  fs::path(dir_results_mixtures, 
           "SOL CHS all Mixtures MWAS results long hyper_g_v2.rds")) %>% 
  modify(~.x %>% mutate(exposure = exposure %>% 
                          rename_pfas(include_asterisk = TRUE)))


## Read in mzrt key--------------------------------------
mzrt_key <- read_rds(fs::path(dir_data_local,
                              "Supporting Files", 
                              "mummichog_pw_ec_feature_key_with_cpd_names.rds"))

# 2) modify MWAS data ------------
# modify mwas data
# Create dataframe and reduce number of rows in data 
mwas_beta_coefs_df  <- mwas_results %>% 
  bind_rows() %>% 
  mutate(beta_ci = jag2::effest_ci(estimate_beta,
                                   lcl_beta,
                                   ucl_beta, n.digits = 2)) %>% 
  filter(mixture == "all pfas")


# Determine the number pips>0.8 for each PFAS compound
# (mwas_beta_coefs_df %>% 
#     group_by(mixture, cohort, exposure) %>% 
#     summarise(pip80 = sum(estimate_pip>0.8), 
#               pct = pip80/length(estimate_pip)) %>% 
#     filter(exposure != "Mixture effect", 
#            mixture == "all pfas") %>% 
#     arrange(-pip80))


# Pivot data wider to get a single row for each mz/rt for each mixture analysis
mwas_beta_coef_w <- mwas_beta_coefs_df %>% 
  pivot_wider(id_cols = c(mode, mixture, feature), 
              names_from = c(cohort, exposure), 
              values_from = c(estimate_beta:q_value, beta_ci)) %>% 
  rename_all(~tolower(.) %>% str_replace_all(" ", "_")) %>% 
  janitor::remove_empty(which = "cols") %>% 
  janitor::clean_names()


# Filter only mzrt in annotated set
mwas_beta_coef_annotated_set <- mwas_beta_coef_w %>% 
  tidylog::filter(feature %in% mzrt_key$name) %>% 
  mutate(sig = case_when(p_value_sol_mixture   < 0.25 & 
                           p_value_chs_mixture < 0.25 ~ "Sig. Both cohorts", 
                         p_value_sol_mixture   < 0.25 ~ "Sig. Solar",  
                         p_value_chs_mixture   < 0.25 ~ "Sig. CHS", 
                         TRUE ~ "Not sig.")) %>% 
  select(mode, mixture, feature, sig, contains("mixture"), everything()) 



# 3) Merge mzrt/ EC key with mwas results -------------------------
# change format for mzrt key pathways for next step 
mzrt_key_all_pathways <- mzrt_key$pathway %>% 
  str_split(., "; ") %>% 
  enframe() %>% 
  rename(pathway_dfs = value) %>% select(-name)

# Filter mzrt associated with pathways sig. in both cohorts to get mzrt_ec_pw key
mzrt_key_sig_pws_only_1 <- bind_cols(mzrt_key,
                                     mzrt_key_all_pathways) %>% 
  unnest(pathway_dfs) %>% 
  tidylog::filter(pathway_dfs %in% sig_pw_both_chrt$path) %>% 
  ungroup() %>%
  select(-pathway, -path_2) %>% 
  rename(pathway = pathway_dfs)

# Create data for upset plot diagrams 
upplot_data <- tidylog::left_join(mzrt_key_sig_pws_only_1, 
                                mwas_beta_coef_annotated_set, 
                                by = c("name" = "feature"))

# Save upset plot data 
write_rds(upplot_data,
          here::here(dir_results_mixtures,
                     "Ecd effect estimates and pathways sol chs for upset.RDS"))

# Reduce across ECs 
mzrt_key_sig_pws_only_2 <- mzrt_key_sig_pws_only_1 %>% 
  tidylog::group_by(empirical_compound, 
                    matched_form,
                    query_mass,
                    retention_time,
                    name) %>% 
  tidylog::summarise(met_name = unique(met_name) %>% 
                       str_c(collapse = "; "), 
                     matched_compound = unique(matched_compound) %>% 
                       str_c(collapse = "; "), 
                     mass_diff =  unique(mass_diff) %>% 
                       str_c(collapse = "; "), 
                     pathway =  unique(pathway) %>% 
                       str_c(collapse = "; "), ) %>% 
  select(empirical_compound, met_name, everything()) %>% 
  ungroup()


## Merge mzrt_key_sig_pws_only with mwas results ----------------------------
t1 <- tidylog::left_join(mzrt_key_sig_pws_only_2,      
                         mwas_beta_coef_annotated_set,   
                         by = c("name" = "feature")) %>% 
  select(mode, mixture, met_name, pathway, sig, everything()) %>% 
  arrange(mixture, pathway, empirical_compound)

# Select only mixture containing all PFAS
all_pfas_t1 <- t1 %>% 
  filter(mixture == "all pfas")

# Get summary of each EC (including min p value for each EC)
summary_annotated_ee <- all_pfas_t1 %>% 
  group_by(mixture, empirical_compound) %>% 
  summarise(n_met_names = length(unique(met_name)), 
            n_pathway = length(unique(pathway)), 
            n_matched_form = length(unique(matched_form)), 
            n_features = length(unique(name)), 
            min_p_sol = min(p_value_sol_mixture),
            min_p_chs = min(p_value_chs_mixture)) %>% 
  mutate(sig = case_when(min_p_sol   < 0.25 & 
                           min_p_chs < 0.25 ~ "Sig. Both cohorts", 
                         min_p_sol   < 0.25 ~ "Sig. Solar",  
                         min_p_chs   < 0.25 ~ "Sig. CHS", 
                         TRUE ~ "Not sig.")) %>% 
  ungroup()

table(all_pfas_t1$pathway)


# Filter only empirical compounds which have at least one mz/rt significant
annotated_sig_ecs_ee <- all_pfas_t1 %>%  
  tidylog::filter(empirical_compound %in% 
                    summary_annotated_ee$empirical_compound) %>% 
  tidylog::filter(sig != "Not sig.") %>% 
  mutate(blank = NA, 
         blank1 = NA, 
         blank2 = NA, 
         casnum = NA, 
         HRE_standard = NA, 
         # query_mass = formatC(query_mass, format = "f", digits = 3), 
         retention_time = formatC(retention_time, format = "f", digits = 1), 
         across(contains("estimate_pip_sol"),
                ~formatC(., format = "f", digits = 2)), 
         across(contains("estimate_pip_chs"),
                ~formatC(., format = "f", digits = 2)), 
         mode = toupper(mode) %>% 
           str_replace("NEG", " Negative") %>% 
           str_replace("POS", " Positive"))

# Save data
write_rds(annotated_sig_ecs_ee,
          here::here(dir_results_mixtures, "Sig annotated metabolite effect estimates sol chs.RDS") )


# Get final table 3 ----------------------
tyrosine_metabolites_pfas_mixture <- annotated_sig_ecs_ee %>% 
  ungroup() %>%
  select(mode, empirical_compound, met_name, casnum, HRE_standard, 
         matched_form, query_mass,  retention_time, sig, 
         beta_ci_sol_mixture, 
         sd_beta_sol_mixture,
         p_value_sol_mixture, 
         blank, 
         beta_ci_chs_mixture, 
         sd_beta_chs_mixture,
         p_value_chs_mixture,
         blank1,
         contains("estimate_pip_sol"), 
         blank2,
         contains("estimate_pip_chs"), 
         everything())