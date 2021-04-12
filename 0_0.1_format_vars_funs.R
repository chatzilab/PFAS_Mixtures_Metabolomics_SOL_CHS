## functions that format variables so always the same over multiple data sets






















# # colors
# colors_fct_liv_inj  <- c("Normal" = "#054ea0", "Injury" = "#f77f0a")
# # long cohort names for presentation
# colors_fct_cohorts <- c(
#   "BIB United Kingdom"  = "#09863f", # forest green
#   "EDEN France" = "#06a7b2", # light sea green
#   "KANC Lithuania" = "#fdb913", # orange
#   "MOBA Norway" = "#ef2b2d", # crimson
#   "RHEA Greece" = "#054ea0", # teal
#   "INMA Spain" = "#c0f451"  # green yellow
#   )
# # short cohort names for coding
# colors_fct_cohorts_short <- colors_fct_cohorts
# names(colors_fct_cohorts_short) <- stringr::word(names(colors_fct_cohorts), 1)
# 
# colors_fct_omics_light <- c(
#   "CpGs" = "#63a3cf", 
#   "Transcripts" = "#ffc283", 
#   "Proteins" = "#c5afe8", 
#   "Serum Metabolites" = "#b0dfab", 
#   "Urine Metabolites" = "#fff89a"
#   )
# colors_fct_omics_dark <- c(
#   "CpGs" = "#5291bd", 
#   "Transcripts" = "#edb071", 
#   "Proteins" = "#966fd6", 
#   "Serum Metabolites" = "#9ecd99", 
#   "Urine Metabolites" = "#ebe687"
# )
# colors_fct_diffexpr <- c(
#   "Significantly Higher in Normal" = "blue", 
#   "Higher in Normal" = "lightblue", 
#   "Significantly No Difference" = "black", 
#   "No Difference" = "gray", 
#   "Higher in Liver Injury" = "lightcoral",
#   "Significantly Higher in Liver Injury" = "red"
# )
# 
# ## helix cohort variables
# fmt_h_cohort <- function(x){factor(x, levels = c("BIB", "EDEN", "KANC","MOBA", "RHEA", "INMA"))}
# 
# 
# ## outcome variables for liver enzymes
# fmt_liv_inj <- function(x){factor(x, levels = 0:1, labels = c("Normal", "Injury"))}
# fmt_liv_enz_ord <- function(x){factor(x, ordered = TRUE, c("alt_sex", "q90_any", "q90_alt", "q90_ast", "q90_ggt", "q90_ck18"))}
# 
# ## genetic SNP exposures
# fmt_snp_type <- function(x){factor(x, 1:2, c("NAFLD-specific", "General Liver Injury"))}
# 
# 
# ## omics groups derived from string prefix
# group_feature_prefix <- function(df, feature_group){
#   dplyr_feature_group = enquo(feature_group)
#   df %>% 
#     mutate(omics_group = case_when(
#       str_detect(!!dplyr_feature_group, "pro_") ~ 1,
#       str_detect(!!dplyr_feature_group, "metser_") ~ 2,
#       str_detect(!!dplyr_feature_group, "meturi_") ~ 3
#     ) %>% factor(1:3, c("Proteins", "Serum Metabolites", "Urine Metabolites"))
#     )
# }
# 
# # remove extraneous info from labels of each feature
# clean_feat_labs <- function(long_label){
#   long_label %>% 
#     str_remove("pro_") %>%     # proteomics prefix
#     str_remove("metser_") %>%  # serum metabolomics prefix
#     str_remove("meturi_") %>%  # urine metabolomics prefix
#     str_remove("log.") %>%     # log prefix
#     str_remove("PC.")       # unknown serum metab prefix
# }
# 
# # remove extraneous sample names "aa" and "ae"
# clean_feat_aa_ae <- function(long_label){
#   long_label %>% 
#   str_remove("aa.") %>%      # unknown serum metab prefix
#   str_remove("ae.")          # unknown serum metab prefix
# }
# 
# # clean excess punctuation
# clean_punct <- function(messy_label){
#   messy_label %>% 
#     str_replace_all("\\.+", "\\.") %>% 
#     str_replace_all("\\.", " ")
# }
# clean_punct("this..that..HeY")
# # scale omics biological features
# 
#   ## normalize biological features
#   norm_bio_features <- function(df, prefix){
#     df %>% mutate_at(vars(starts_with(prefix)), function(x) as.double(scale(x)))
#   }
#   # lst_proteome$proteome %>% norm_bio_features("pro_") %>% glimpse
#   
#   ## "Pareto scaling...uses the square root of the standard deviation as the scaling factor."
#   ## "...Pareto scaling...is recommended for metabolomics data." - Grace and Hudson (2016) 
#   paretoscale <- function(df, prefix) {
#     # function adapted from Grace and Hudson (2016) 
#     paretoscale_col <- function(i_col){
#       col_center <- i_col - mean(i_col, na.rm = TRUE) # mean center
#       col_scale <- col_center / sqrt(sd(i_col, na.rm = TRUE)) # divide by sqrt sd 
#       return(col_scale)
#     }
#     # apply Pareto scaling to all variables matching prefix
#     df %>% mutate_at(vars(starts_with(prefix)), paretoscale_col)
#   }
#   # lst_proteome$proteome %>% paretoscale("pro_") %>% glimpse
# 
