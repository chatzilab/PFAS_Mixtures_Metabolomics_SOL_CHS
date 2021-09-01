# Functions

# Create Functions for summarizing Data ---------------
transpose_ft <- function(ft) {
  dataout <- ft %>%
    mutate(name = str_c(mz, time, sep = "_")) %>%
    select(name, everything(), -mz, -time) %>%
    gather(file_name, val, 2:ncol(.)) %>%
    spread(name, val) %>% 
    mutate(file_name = str_remove(file_name, "_mz_xml"))
  
  return(dataout)
}

# # Rename Compounds -----------------------
rename_pfas <- function(pfas_names, include_asterisk = FALSE, arrange_by_class = FALSE){
  x <- tibble(pfas = pfas_names)
  pfas2 <-  x %>%
    mutate(pfas2 = case_when(
      pfas == "pfhxs" ~ "PFHxS",
      pfas == "pfhps" ~ "PFHpS",
      pfas == "pfpes" ~ "PFPeS",
      pfas == "pfhpa" ~ "PFHpA",
      pfas == "nmefosaab" ~ "N-MeFOSAA-b†", 
      pfas == "pfuda" ~ "PFUnDA†",
      pfas == "pfds" ~ "PFDS†",
      pfas == "netfosaa" ~ "N-EtFOSAA†",
      pfas == "pfns" ~ "PFNS†",
      pfas == "pfbs" ~ "PFBS†",
      pfas == "x82fts" ~ "8:2 FTS†", 
      pfas == "pfhxa" ~ "PFHxA†", 
      pfas == "pfdoa" ~ "PFDoDA†",
      TRUE ~ toupper(pfas)) %>% 
        as.factor() %>% 
        fct_relevel(., 
                    "PFOS", 
                    "PFOA", 
                    "PFHxS", 
                    "PFNA", 
                    "PFHpS",
                    "PFDA", 
                    "PFPeS", 
                    "PFHpA",
                    "N-MeFOSAA-b†",
                    "N-EtFOSAA†",
                    "PFDS†",
                    "PFBS†", 
                    "8:2 FTS†", 
                    "PFDoDA†", 
                    "PFUnDA†",
                    "PFNS†",
                    "PFHxA†")) 
  
  if(include_asterisk == TRUE){ 
    pfas2 <-  x %>%
      mutate(pfas2 = case_when(
        pfas == "pfhxs" ~ "PFHxS",
        pfas == "pfhps" ~ "PFHpS",
        pfas == "pfpes" ~ "PFPeS",
        pfas == "pfhpa" ~ "PFHpA",
        pfas == "nmefosaab" ~ "N-MeFOSAA-b*", 
        pfas == "pfuda" ~ "PFUnDA*",
        pfas == "pfds" ~ "PFDS*",
        pfas == "netfosaa" ~ "N-EtFOSAA*",
        pfas == "pfns" ~ "PFNS*",
        pfas == "pfbs" ~ "PFBS*",
        pfas == "x82fts" ~ "8:2 FTS*", 
        pfas == "pfhxa" ~ "PFHxA*", 
        pfas == "pfdoa" ~ "PFDoDA*",
        TRUE ~ toupper(pfas)) %>% 
          as.factor() %>% 
          fct_relevel(., 
                      "PFOS", 
                      "PFOA", 
                      "PFHxS", 
                      "PFNA", 
                      "PFHpS",
                      "PFDA", 
                      "PFPeS", 
                      "PFHpA",
                      "N-MeFOSAA-b*",
                      "N-EtFOSAA*",
                      "PFDS*",
                      "PFBS*", 
                      "8:2 FTS*", 
                      "PFDoDA*", 
                      "PFUnDA*",
                      "PFNS*",
                      "PFHxA*")) 
  }
  
  if(arrange_by_class == TRUE){ 
    pfas2 <-  pfas2 %>% 
      left_join(lod, by = "pfas") %>% 
      mutate(pfas2 = fct_reorder(pfas2, order_by_class))
  }
  
  return(pfas2$pfas2)
}





# # Rename outcomes -----------------------
# rename_outcomes <- function(outcome, include_all_ogtt_timepoints=FALSE ){
#   outcome2 = case_when(
#     outcome == "og_glu_5"	 ~ "Fasting glu",
#     outcome == "glu15ave"   ~ "15 min glu",
#     outcome == "og_glu30"   ~ "30 min glu",
#     outcome == "glu45ave"   ~ "45 min glu",
#     outcome == "og_glu60"   ~ "60 min glu",
#     outcome == "og_glu90"   ~ "90 min glu",
#     outcome == "og_glu120"  ~ "2-hour glu",
#     outcome == "max_glu"	   ~ "Max Glucose",
#     outcome == "guac"	     ~ "Glu AUC",
#     outcome == "glu_inc_auc"~ "Glu incremental AUC",
#     outcome == "og_ins_5"	 ~ "Fasting ins",
#     outcome == "ins15"      ~ "15 min ins",
#     outcome == "og_ins30"   ~ "30 min ins",
#     outcome == "ins45"      ~ "45 min ins",
#     outcome == "og_ins60"   ~ "60 min ins",
#     outcome == "og_ins120"	 ~ "2-hour ins",
#     outcome == "max_ins"	   ~ "Max Insulin",
#     outcome == "iuac"	     ~ "Ins AUC",
#     outcome == "ins_inc_auc"~ "Ins incremental AUC",
#     outcome == "matsuda"	   ~ "ISI Matsuda",
#     outcome == "homa"	     ~ "HOMA-IR",
#     outcome == "cubert_igi" ~ "IGI", 
#     outcome == "a1c"	        ~ "HbA1c") %>% as.ordered()
#   
#   if(include_all_ogtt_timepoints == TRUE){
#     outcome3 = fct_relevel(
#       outcome2, 
#       "Fasting glu", "15 min glu",
#       "30 min glu", "45 min glu",
#       "60 min glu", "90 min glu",
#       "2-hour glu","Max Glucose",
#       "Glu AUC", "Glu incremental AUC",
#       "Fasting ins","15 min ins","30 min ins","45 min ins",
#       "60 min ins","2-hour ins","Max Insulin",
#       "Ins AUC", 
#       "ISI Matsuda",
#       "HOMA-IR",
#       "IGI", 
#       "HbA1c")
#   } 
#   else{
#     outcome3 = fct_relevel(
#       outcome2, 
#       "Fasting glu",
#       "30 min glu",
#       "60 min glu",
#       "2-hour glu",
#       "Glu AUC", 
#       "Ins AUC",
#       "ISI Matsuda",
#       "HOMA-IR",
#       "IGI", 
#       "HbA1c")
#   }
#   
#   
#   
# }
