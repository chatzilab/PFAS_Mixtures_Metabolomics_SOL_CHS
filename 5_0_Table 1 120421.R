# Table 1
source(here::here("other_scripts", "Create Parental Education SES vars.R"))

library(jag2) # Contains personal functions

# Combine baseline data from SOLAR and CHS -----------------------------

## Get baseline data in SOLAR -------------------

# Rejoin with baseline data 
solar_eo <- tidylog::left_join(exposure_outcome$solar, solar_edu)

solar_eo$hispanic <- "Hispanic"
solar_eo$cohort <- "SOLAR"

# Select Key Vars
solar_eo <- solar_eo %>% 
  select(cohort, sex, age, tot_pf, wave,
         bmi_og, hispanic, edu_house, modified2factor_score, puberty, 
         og_glu_5, og_glu120, guac, iuac, a1c) %>%
  rename(bmi = bmi_og)



## Baseline data in CHS --------------------

# Subset CHS Data
chs_eo <- exposure_outcome$chs %>% 
  tidylog::left_join(., chs_edu) %>%
  select(cohort, sex, age, 
         bmi, hispanic, edu_house, 
         og_glu_5, og_glu120, guac, iuac, a1c) #pfoa, 


chs_eo$cohort <- "CHS" 
chs_eo$puberty <- ""
chs_eo$delta_years <- 0
chs_eo$puberty <- "" 
chs_eo$n_visits <- 1


# Bind rows, and change pfhxs < lod to NA
t1_data <- bind_rows(solar_eo, chs_eo) 


# Table 1 ------------------------------------------
## Means and geom_means --------------------------
table1_a <- t1_data %>% 
  group_by(cohort) %>% 
  summarise(
    "General Characteristics" = "", 
    "Sample Size" = length(sex) %>% as.character(), 
    "Sex, Female [n (%)]" = npct(sex, "Female", n.digits = 1), 
    "Age, years" = avg_sd_fxn(age, n.digits = 1), 
    "BMI, kg/m2" = avg_sd_fxn(bmi, n.digits = 1), 
    # "Body Fat (%)" = avg_sd_fxn(tot_pf, n.digits = 1, include.n = TRUE), 
    "Puberty status [n (%)]" = "",
    "Pre Puberty (Tanner stage 1)" = 
      paste0(sum(puberty== "Pre puberty"), " (",
             100*round(sum(puberty== "Pre puberty")/length(puberty), digits = 2), 
             "%)"), 
    "Puberty (Tanner stages 2-4)" =
      paste0(sum(puberty== "Puberty"), " (",
             100*round(sum(puberty== "Puberty")/length(puberty), digits = 2),
             "%)"),
    "Post Puberty (Tanner stage 5)" =  
      paste0(sum(puberty == "Post puberty"), " (",
             100*round(sum(puberty== "Post puberty")/length(puberty), digits = 2),
             "%)"),
    "Ethnicity [n (%)]" = "",
    "Hispanic" = paste(sum(hispanic =="Hispanic"), 
                       " (", 
                       100*round(sum(hispanic =="Hispanic")/length(hispanic), digits = 2),
                       "%)", sep = ""),
    "Non-Hispanic" = paste(sum(hispanic !="Hispanic"), 
                           " (", 
                           100*round(sum(hispanic !="Hispanic")/length(hispanic), digits = 2),
                           "%)", sep = ""),
    
    "Socioeconomic Status" = "",
    "Modified Hollingshead Four-Factor Index" = avg_sd_fxn(modified2factor_score, n.digits = 1),
    # Education
    "Household education level [n (%)]" = "",
    "Did not graduate High School" = npct(edu_house, "Did not graduate high school", 2),
    "High School Graduate" = npct(edu_house, "High school graduate", 2),
    "Partial college (at least one year) or specialized training" = 
      npct(edu_house, "Partial college (at least one year) or specialized training", 2),
    "Completed college/university" = 
      npct(edu_house, "Standard college or university graduation", 2),
    "Graduate professional training" = 
      npct(edu_house, "Graduate professional training (graduate degree)", 2),
    "Missing" = 
      npct(edu_house, "Missing", 2),
    "Study Wave [n (%)]" = "",
    "Wave 1 (2001-2003)" = npct(wave, "first wave", 2),
    "Wave 2 (2010-2012)" = npct(wave, "second wave", 2)) %>% 
  pivot_longer(cols = "General Characteristics":"Wave 2 (2010-2012)")


# ## Table 1 Sample Size -----------------------------------------------------
# table1_n <- t1_data %>% 
#   group_by(cohort, sex) %>% 
#   summarise(
#     "Sample Size" = length(sex), 
#     "Age, years" = sum(!is.na(age)), 
#     "BMI, kg/m2" = sum(!is.na(bmi)), 
#     "Modified Hollingshead Four-Factor Index" = sum(!is.na(modified2factor_score)),
#     "Puberty status [n (%)]" = sum(!is.na(puberty))
#     ) %>% 
#   pivot_longer(cols = "Puberty status [n (%)]", values_to = "N")

## Table 1 p-value  -----------------------------------------------------
# table1_p <- t1_data %>% 
#   # group_by(sex) %>% 
#   nest() %>% 
#   summarise(
#     "BMI, kg/m2" = map_chr(
#       data, ~t.test(bmi ~ cohort,data = ., paired = FALSE)$p.value %>% 
#         signif(., 2) %>% as.character(.)) , 
#     "Household education level [n (%)]" = 
#       map_chr(data, ~chisq.test(.$edu_house , .$cohort)$p.value %>%
#                 signif(., 2) %>% 
#                 as.character(.))) %>%
#   pivot_longer(cols = "BMI, kg/m2":"Household education level [n (%)]",
#                values_to = "p")


## Table 1 p-value overall  ------------------------------------------------
table1_p_overall <- t1_data %>% 
  nest(data = everything()) %>% 
  summarise(
    "BMI, kg/m2" = map_chr(
      data, ~t.test(bmi ~ cohort,data = ., paired = FALSE)$p.value %>% 
        signif(., 2) %>% as.character(.)) , 
    "Household education level [n (%)]" = 
      map_chr(data, ~chisq.test(.$edu_house , .$cohort)$p.value %>%
                signif(., 2) %>% 
                as.character(.)), 
    "Sex, Female [n (%)]" = 
      map_chr(data, ~chisq.test(.$sex , .$cohort)$p.value %>%
              signif(., 2) %>% 
              as.character(.))) %>%
  pivot_longer(cols = "BMI, kg/m2":"Sex, Female [n (%)]", 
               values_to = "p_overall")


# Join and format Table 1 ------------

# Merge with sample size data
# table1_b <- left_join(table1_a, table1_n)

## Subset and format CHS -----------------
chs_t1 <- table1_a %>% 
  filter(cohort == "CHS") %>% 
  mutate(value = case_when(str_detect(name, "Tanner") ~ "",
                           str_detect(name, "Modified Hollingshead") ~ "", 
                           str_detect(name, "Wave") ~ "",
                           TRUE ~ value)) %>% 
  rename(CHS = value)



## Subset and format SOLAR ----------------
solar_t1 <- table1_a %>% 
  filter(cohort == "SOLAR")  %>% 
  rename(SOLAR = value)

# Join SOLAR, CHS, and overall p-values ----------------------
table1_w <- full_join(solar_t1, 
                      chs_t1, 
                      by = "name") %>% 
  full_join(., table1_p_overall) %>% 
  mutate(p_overall = replace_na(p_overall, "")) %>% 
  select(name, 
         "SOLAR",
         "CHS",
         p_overall)


# Save File
writexl::write_xlsx(table1_w, 
                    here::here(dir_reports, 
                               "Table 1 Wide.xlsx"))



# rm(avg_sd_fxn,
#    fungm,
#    summary_data,
#    chs_data,
#    table1_a,
#    table1_b,
#    chs_t1,
#    solar_t1)



