# Set exposure outcome vars

# Set important vars
exposure_type = "PFAS"

exposures = c("netfosaa","nmefosaab","pfbs","pfda","pfdoa","pfds","pfhpa",
              "pfhps","pfhxa","pfhxs","pfna","pfns","pfoa","pfos","pfpes",
              "pfuda", "x82fts")

# exposures_continuous_old <- c("pfda","pfhpa",
#                           "pfhps","pfhxs","pfna","pfoa","pfos","pfpes")

exposures_continuous <- c("pfda","pfhps","pfhxs","pfna","pfoa","pfos")


exposures_for_analysis <- c("lg2_pfda",
                            "lg2_pfhpa",
                            "lg2_pfhps",
                            "lg2_pfhxs",
                            "lg2_pfna",
                            "lg2_pfoa",
                            "lg2_pfos",
                            "lg2_pfpes",
                            "nmefosaab",
                            "netfosaa", 
                            "pfds", 
                            "pfbs", 
                            "x82fts") 

modes = c("c18pos","c18neg", "hilicpos", "hilicneg")

cohort = c("solar", "chs")

exp_cont_below_lod_na_not_trns <- str_c(exposures_continuous, "_w_na")


# exposures_oc <- c("hexachlorobenzene_impute", "dde_impute", "ocs",
#                   "pbde_154_impute", "pbde_47_impute",
#                   "pbde_100_ngml_detect", "pbde_153_ngml_detect",
#                   "pbde_85_ngml_detect", "pcb_118_ngml_detect",
#                   "pcb_138_ngml_detect", "pcb_153_ngml_detect",
#                   "pcb_180_ngml_detect",
#                   "pcb_num_detect", 
#                   "pbde_num_detect")

# outcomes = list(
#   glucose_outcomes   = c("og_glu_5","og_glu30","og_glu60","og_glu120","guac"),
#   insulin_outcomes   = c("og_ins_5","og_ins30","og_ins60","og_ins120","iuac","a1c"),
#   lipid_outcomes     = c("tag","tot_chol","og_ffa_5","og_ffa120"),
#   glucose_homeostasis= c("homa","matsuda","cubert_igi","si","sg","air","di"),
#   body_composition   = c("bmi","saat","iaat","tot_pf","tot_fat_kg","tot_lm_kg"),
#   proteins           = c("ttesto","ftesto","estrad","dhs",
#                          "shbg","il6","tnfa","igf","il1b","il8",
#                          "leptin","cortisol","adiponec"))