
View(single_pfas_effect_est)

tyrosine_met <- single_pfas_effect_est %>% 
  mutate(mz = str_split_fixed(name, "_", n = 2)[,1] %>% 
           as.numeric() %>%
           formatC(., format = "f", digits = 5), 
         rt = str_split_fixed(name, "_", n = 2)[,2] %>% 
           as.numeric() %>%
           formatC(., format = "f", digits = 1)) %>% 
  pivot_longer(cols = estimate_solar:conf_high_chs, 
               names_to = "variable_name", 
               values_to = "value") %>% 
  mutate(cohort = case_when(str_detect(variable_name, "solar") ~ "SOLAR", 
                            str_detect(variable_name, "chs") ~ "CHS") %>% 
           fct_relevel("SOLAR", "CHS"), 
         variable_name = str_remove(variable_name, "_solar") %>% 
           str_remove(., "_chs")) %>% 
  pivot_wider(names_from = variable_name, 
              values_from = value) 



# Hand Edit some metabolites
tyrosine_met2 <- tyrosine_met %>% 




table(tyrosine_met2$variable_name
)
colnames(tyrosine_met2)


# L-Tyrosine ------------------------
# met_specific_data <- tyrosine_met %>% 
#   filter(met_name_first_only == met_name_first_only)


met_name_first_only_for_fig <- i 

min(tyrosine_met2$conf_low)

# Create Funciton
create_pfhxs_tyrosine_figs <- function(data, met_name_first_only_for_fig){
  
  met_specific_data <- data %>% 
    filter(met_name_first_only == met_name_first_only_for_fig)
  
  (met_specific_fig <- ggplot(met_specific_data, 
                              aes(x = cohort, 
                                  y = estimate) ) + 
      geom_point(size = 2.0) +
      geom_errorbar(aes(ymin =  conf_low,
                        ymax = conf_high),
                    width = 0) + 
      geom_hline(yintercept = 0, linetype = 2) +
      scale_y_continuous(limits = c(-.41, 0.41), name = "Beta (95% CI)") + 
      theme(axis.title.x =  element_blank(),
            legend.position = "right") + 
      labs(title = met_name_first_only_for_fig, 
           subtitle = str_c("mz = ", met_specific_data$mz, "\n", 
                            "retention time = ", met_specific_data$rt, " (s)")))
  
  ggsave(met_specific_fig, 
         filename =
           fs::path(dir_reports,
                    "pfhxs tyrosine metabolites",
                    str_c(str_replace_all(met_name_first_only_for_fig, 
                                          "[[:punct:]]", " "),
                          " and pfhxs",
                          ".jpg")),
         height = 5, width = 3)
}

# Run for loop

for(i in unique(tyrosine_met2$met_name_first_only)){
  create_pfhxs_tyrosine_figs(tyrosine_met2, i)
}


# 1 off figures
create_pfhxs_tyrosine_figs(tyrosine_met2, "Acetoacetic acid")
create_pfhxs_tyrosine_figs(tyrosine_met2, "Vanylglycol")
