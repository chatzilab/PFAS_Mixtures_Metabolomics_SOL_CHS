# Load annotated data ----------------------------------------------

  ## C18 Neg ------------
c18neg_annotate <- read.delim(
  fs::path(dir_data,
           "4_Common_Metabolites_Annotation",
           "HRE0026_C18neg_CommonMetabolites_SOLAR-CHS_V1.txt")) %>%
  janitor::clean_names()

## C18 Pos ------------
c18pos_annotate <- read.delim(
  fs::path(dir_data,
           "4_Common_Metabolites_Annotation",
           "HRE0026_C18pos_CommonMetabolites_SOLAR-CHS_V1.txt")) %>%
  janitor::clean_names()

## hilic Neg ------------
hilicneg_annotate <- read.delim(
  fs::path(dir_data,
           "4_Common_Metabolites_Annotation",
           "HRE0026_HILICneg_CommonMetabolites_SOLAR-CHS_V1.txt")) %>%
  janitor::clean_names()

## hilic Pos ------------
hilicpos_annotate <- read.delim(
  fs::path(dir_data,
           "4_Common_Metabolites_Annotation",
           "HRE0026_HILICpos_CommonMetabolites_SOLAR-CHS_V1.txt")) %>%
  janitor::clean_names()
