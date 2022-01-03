# SOLAR --------------------------------
solar_edu <- read_rds(
  file = "R:/SOLAR EDCs/11.Jesse EDCs T2D 2021/Datasets/Outcome Data/Longitudinal outcomes and PFAS_v01.rds") %>% 
  select(id, visit, edu_house, edumom, edudad) %>% 
  group_by(id) %>% 
  filter(visit == min(visit)) %>% 
  ungroup() %>% 
  select(-visit)

# To Match CHS, truncate data so that "less than highschool" is the lowest option.
solar_edu <- solar_edu %>% 
  mutate(sch_mom = if_else(edumom < 4, 3, edumom), 
         sch_dad = if_else(edudad < 4, 3, edudad)) %>% 
  rowwise(id) %>% 
  mutate(edu_house = mean(c(sch_mom, sch_dad), na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(edu_house = case_when(edu_house == 3 ~ "Did not graduate high school", 
                               edu_house <= 4 ~ "High school graduate",  
                               edu_house <= 5 ~ "Partial college (at least one year) or specialized training",
                               edu_house <= 6 ~ "Completed college/university",
                               edu_house <= 7 ~ "Graduate professional training (graduate degree)"
                               ) %>% 
           replace_na(., "Missing")) %>% 
  select(id, edu_house)


# CHS -----------------------------
chs_edu <- exposure_outcome$chs %>% 
  mutate(across(c(sch_mom, sch_dad), 
                ~if_else(. == 9, 
                         NA_real_, 
                         as.numeric(.))), 
         sch_mom = sch_mom+2, 
         sch_dad = sch_dad+2) %>% 
  rowwise(id) %>% 
  mutate(edu_house = mean(c(sch_mom, sch_dad), na.rm = TRUE)) %>%
  ungroup() %>% 
  select(id, edu_house) %>% 
  mutate(edu_house = case_when(edu_house == 3 ~ "Did not graduate high school", 
                               edu_house <= 4 ~ "High school graduate",  
                               edu_house <= 5 ~ "Partial college (at least one year) or specialized training",
                               edu_house <= 6 ~ "Standard college or university graduation",
                               edu_house <= 7 ~ "Graduate professional training (graduate degree)") %>% 
           replace_na(., "Missing")) %>% 
  select(id, edu_house)

