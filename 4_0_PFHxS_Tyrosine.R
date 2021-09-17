# Identify main pathways of interest: Tyrosine and PFHxS
library(MetaboAnalystR)
library(tidylog)
library(RColorBrewer)
library(dendextend)
library(gplots)


# Read in MWAS Beta Coefficients --------------------
mwas_results_long <- read_rds(
  here::here("Temporary results", 
             "PFAS", 
             "SOL CHS all MWAS results long.rds")) 


# Scale estimates 
mwas_results_long <- mwas_results_long %>% 
  modify(~.x %>% 
           group_by(exposure) %>% 
           mutate(estimate_scaled = scale(estimate, center = FALSE))) 


# Pivot wider, change non-sig estimates to 0
mwas_results_wide <- mwas_results_long %>% 
  modify(~mutate(.x, estimate_scaled = if_else(p_value < 0.05, 
                                               estimate_scaled, 
                                               0)) %>%
           # filter(str_detect(exposure, "lg2"), 
           #        str_detect(exposure, "pfpes", negate = TRUE)) %>%
           select(exposure, name, estimate_scaled) %>% 
           pivot_wider(id_cols = name, 
                       names_from = exposure, 
                       values_from = estimate_scaled))



# Read in Dougs annotations ------------------------------------------------
annotations <- read_rds(
  fs::path(dir_data, 
           "4_Common_Metabolites_Annotation", 
           "Common_Metabolites_SOLAR_CHS_V1.RDS")) %>% 
  bind_rows(.id = "mode")


tem <- read_csv(here::here("Temporary results", 
                    "Pathway Results Annotated Metabolites", 
                    "All Annotated Metabolites Hand Curated.csv"))

tyrosine_pw <- readxl::read_xlsx(here::here("Temporary results", 
                                "Pathway Results Annotated Metabolites", 
                                "Pathway Metabolites",
                                "Tyrosine Pathway.xlsx"), col_names = FALSE) %>% 
  clean_names()

# Get list of all metabolites
cmpd.vec <- annotations$refmet_name %>% 
  str_split("; ") %>% 
  unlist() %>% 
  unique()


# Create name variable to match mwas results
annotations2 <- annotations %>% 
  mutate(name = str_c(mz, time, sep = "_"),
         across(refmet_name:number_studies_reported, ~str_split_fixed(., ";", n = Inf)[,1]))





# Tyrosine pathway
cpdnms <- c("2,5-Dihydroxybenzoate",
            "3-(4-Hydroxyphenyl)pyruvate",
            "3,4-Dihydroxy-L-phenylalanine",
            "3-Methoxy-4-hydroxymandelate",
            "3-Methoxytyramine",
            "4-Hydroxyphenylacetate",
            "Acetoacetate",
            "Dopamine",
            "Fumarate",
            "Homogentisate",
            "Homovanillate",
            "L-Adrenaline",
            "L-Normetanephrine",
            "L-Tyrosine",
            "Pyruvate",
            "Thyroxine",
            "Tyramine")
