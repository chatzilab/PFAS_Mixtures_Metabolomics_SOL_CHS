# Read in mummichog pathway results, get tables of mummichog superpathway results

key = readxl::read_xlsx(here::here("Supporting files",  
                                   "superpathway_key.xlsx"))

# CHS ---------------

chs_mummichog_result_files <- list.files(here("Temporary results", 
                                              "mum_pathway_results", 
                                              "chs"))

chs_mum_results <- map2_dfr(fs::path(here("Temporary results", 
                                          "mum_pathway_results", 
                                          "chs"), 
                                     chs_mummichog_result_files), 
                            chs_mummichog_result_files,
                            ~read_csv(.x) %>% 
                               mutate(exposure = .y)) 


chs_mum_results <- left_join(chs_mum_results, key) %>% 
   mutate(cohort = "CHS")   %>%
   mutate(exposure = exposure %>% 
             str_remove("_mum_pathway_c18.csv") %>% 
             str_remove("_mum_pathway_hilic.csv") %>% 
             str_remove("_ngml") %>%  
             str_remove("_impute") %>% 
             str_replace("hexachlorobenzene", "HCB") %>% 
             str_replace("_detect", "*") %>% 
             toupper() %>% 
             str_replace_all("_", " ")) 


# SOLAR -----------------
solar_mummichog_result_files <-  list.files(here("Temporary results", 
                                                 "mum_pathway_results", 
                                                 "solar"))


solar_mum_results <- map2_dfr(fs::path(here("Temporary results", 
                                            "mum_pathway_results", 
                                            "solar"), 
                                       solar_mummichog_result_files), 
                              solar_mummichog_result_files,
                              ~read_csv(.x) %>% 
                                 mutate(exposure = .y))

rm(solar_mummichog_result_files, chs_mummichog_result_files)


solar_mum_results <- left_join(solar_mum_results, key) %>%  
   mutate(cohort = "solar")  %>%
   mutate(exposure = exposure %>% 
             str_remove("_mum_pathway_c18.csv") %>% 
             str_remove("_mum_pathway_hilic.csv") %>% 
             str_remove("_ngml") %>%  
             str_remove("_impute") %>% 
             str_replace("hexachlorobenzene", "HCB") %>% 
             str_replace("_detect", "*") %>% 
             toupper() %>% 
             str_replace_all("_", " ")) 


matching_pathways <- inner_join(solar_mum_results, chs_mum_results, by = c("exposure", "pathway", "mode")) %>% 
   mutate(exp_pth = paste(exposure, pathway, sep = "_")) %>% 
   select(exp_pth)


# Results Tables --------------------
pathway_results1 <- bind_rows(solar_mum_results, 
                              chs_mum_results) %>% 
   filter(exposure %in% c("DDE", "HCB", "PBDE NUM*", "PCB NUM*")) %>%
   mutate(exp_pth = paste(exposure, pathway, sep = "_")) %>% 
   filter(exp_pth %in% matching_pathways$exp_pth)


pathway_results1 <- pathway_results1 %>% 
   mutate(enrichment = hits_sig/pathway_total)



# filter out gamma > 0.05
pathway_results2 <- pathway_results1 %>% 
   filter(gamma < 0.05, 
          hits_sig > 2)  %>%
   mutate(negLogP = -1 * log(gamma)) %>% 
   janitor::clean_names()


# get names

# Merge with superpathway key
pathway_results3 <- left_join(pathway_results2, key)

# Mutate superpathway 
pathway_results <- pathway_results3 %>% 
   mutate(super_pathway = fct_infreq(super_pathway) %>% 
             fct_lump_lowfreq(), 
          pathway = fct_infreq(pathway) %>% 
             fct_lump_lowfreq())





(bubble_plot<-ggplot(pathway_results %>% 
                        filter(exposure == "HCB"), 
                     aes(x = neg_log_p, y = pathway)) + 
      geom_point(aes(size = hits_sig, color = enrichment)) +
      geom_vline(xintercept = -1 * log(.05)) +
      # scale_colour_brewer(   type = "seq",
      #                        palette = 1,
      #                        direction = 1,
      #                        aesthetics = "colour"
      # ) +
   theme_bw() +
      # scale_x_continuous(limits = c(1.2, 3.2)) +  
      xlab(expression(-log(p))) +
      ylab(NULL) +
      ggtitle("") +
      facet_grid(super_pathway ~ cohort, scales = "free", space = "free")
) 





###############




















# SOLAR C18
table1::table1(~ super_pathway | exposure, 
               pathway_results3 %>% 
                  filter(mode == "c18", 
                         cohort == "solar"), overall = FALSE)


# SOLAR HILIC
table1::table1(~ super_pathway | exposure, 
               pathway_results1 %>% 
                  filter(mode == "hilic", 
                         cohort == "solar"), 
               overall = FALSE)



# chs C18
table1::table1(~ super_pathway | exposure, 
               pathway_results1 %>% 
                  filter(mode == "c18", 
                         cohort == "CHS"), overall = FALSE)


# chs HILIC
table1::table1(~ super_pathway | exposure, 
               pathway_results1 %>% 
                  filter(mode == "hilic", 
                         cohort == "CHS"), 
               overall = FALSE)



# CHS AND SOLAR
table1::table1(~ pathway | exposure,
               pathway_results1 %>% 
                  filter(exposure %in% 
                            c("OCS", 
                              "PBDE NUM*", 
                              "PCB NUM*")))



# How much overlap between c18 and hilic? 
xxx <- pathway_results1 %>% 
   group_by(exposure) %>% 
   summarise(n_total_pathways = length(pathway), 
             n_unique_pathways =length(unique(pathway)),
             n_overlapping_pathways = n_total_pathways - n_unique_pathways)

