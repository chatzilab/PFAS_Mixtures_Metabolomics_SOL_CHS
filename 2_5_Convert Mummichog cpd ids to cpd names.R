# Compound ID to Compound Names
library(janitor)

# read in key linking -----------------------------------
file_location <- fs::path(
  dir_data_local,
  "Supporting Files", 
  "Cpd id to name keys",
  "Mummichog cpd id to cpd name key with db matches.xlsx")


# Modify annotated data to get common compound names ------------------------
library(tidylog)
# read annotated cmpd data
annotated_fts_from_mum <- read_rds(fs::path(dir_data_local,
                                            "Supporting Files", 
                                            "mummichog_pw_ec_feature_key_cpd_id_only.rds"))


# Read in metaboanalyst key 
nm_conv_metaboanalyst <- readxl::read_xlsx(file_location, 
                                           sheet = "metaboanalyst", 
                                           na = "NA") %>% 
  clean_names() %>% 
  rename(chem_id=query) %>% 
  filter(comment == 1) # Comment = 1 for compounds with a match

# Read in metanet key
nm_conv_metanet <- readxl::read_xlsx(file_location, 
                                     sheet = "MetaNetX") %>% 
  clean_names() %>% 
  rename(chem_id = number_query)

# Read in mbrole name conversions
nm_conv_mbrole <- readxl::read_xlsx(file_location, 
                                    sheet = "mbrole") %>% 
  clean_names()  %>% 
  rename(chem_id=input) 


# Modify mbrole 
mbrole_w <- nm_conv_mbrole %>% 
  group_by(chem_id, output_source) %>% 
  summarise(input_source = paste(unique(input_source),collapse = "; "), 
            output = str_c(unique(output), collapse = "; ")) %>% 
  pivot_wider(names_from = "output_source", 
              values_from = c("output") , 
              values_fn = list) %>% 
  janitor::clean_names() %>%
  unnest(cas:lipid_maps) %>% 
  select(-ymdb)


# Combine all annotations
name_cpd_key <- full_join(nm_conv_metaboanalyst, 
                          mbrole_w, by = "chem_id", 
                          suffix = c("_metab", "_mbrole")
) %>%
  full_join(nm_conv_metanet,  by = "chem_id", suffix = c("", "_metanet"))

# Merge matching compounds across databases 
name_cpd_key2 <- name_cpd_key %>% 
  select(-input_source) %>%  
  mutate(hmdb = if_else(is.na(hmdb_metab), hmdb_mbrole, hmdb_metab)) %>% 
  select(chem_id, match, name, hmdb, everything(), -hmdb_mbrole, -hmdb_metab) %>% 
  mutate(name_2 = if_else(is.na(match), name, match)) %>% 
  select(chem_id, match, hmdb, everything(), -name) %>% 
  janitor::remove_empty(which = c("rows", "cols")) %>%
  rename(matched_compound = chem_id, 
         met_name = name_2)

colnames(name_cpd_key2)

# Join Annotated fts from mumichog with names of compounds 

final_annotated_cpds <- left_join(annotated_fts_from_mum, name_cpd_key2)

# Change any remaining missing values to the compound id
final_annotated_cpds <- final_annotated_cpds %>% 
  mutate(met_name = if_else(is.na(met_name), 
                            matched_compound, 
                            met_name))

# Save Data
write_rds(final_annotated_cpds,
          fs::path(dir_data_local,
                   "Supporting Files", 
                   "mummichog_pw_ec_feature_key_with_cpd_names.rds"))


# 
# # Expand annotations from original dataset ---------------------------------------
# annotated_cpds <- annotated_fts_from_mum$solar %>% 
#   select(name:mz_max) %>% 
#   group_by(name) %>%
#   nest() %>%
#   mutate(
#     refmet_name = map(data,  ~ str_split(.x$refmet_name, "; ")), 
#     super_class = map(data,  ~ str_split(.x$super_class, "; ")), 
#     main_class  = map(data,  ~ str_split(.x$main_class,  "; ")), 
#     sub_class   = map(data,  ~ str_split(.x$sub_class,   "; ")), 
#     pubchem_cid = map(data,  ~ str_split(.x$pubchem_cid, "; ")), 
#     number_studies_reported = map(data, ~str_split(.x$number_studies_reported, "; ")))
# 
# # Expand across column lists
# ac_2 <- annotated_cpds %>% 
#   select(-data) %>% 
#   tidyr::unnest(refmet_name:number_studies_reported) %>% 
#   tidyr::unnest(refmet_name:number_studies_reported)
# 
# 
# # Combined with metaboanalyst cpds 
# ac_3 <- ac_2 %>% 
#   left_join(hc_met, by = "refmet_name")
# 
# # Combine with mbrole
# ac_4 <- ac_3 %>% 
#   left_join(mbrole, by = "refmet_name")
# 
# 
# 
# ## Join names from different dbs 
# ac_5 <- ac_4 %>% 
#   rename(number_studies_reported = number_studies_reported.x) %>%
#   mutate(kegg = if_else(is.na(kegg.x),
#                         as.character(kegg.y), 
#                         as.character(kegg.x)), 
#          pubchem_cid = if_else(is.na(pubchem_cid.x), 
#                                as.character(pubchem_cid.y), 
#                                as.character(pubchem_cid.x)),
#          hmdb = if_else(is.na(hmdb.x),
#                         as.character(hmdb.y), 
#                         as.character(hmdb.x)),  
#   ) %>% 
#   select(-c(kegg.y, kegg.x, pubchem_cid.y, pubchem_cid.x, 
#             hmdb.x, hmdb.y, number_studies_reported.y))
# 
# # Join with og data: 
# annotations_only <- left_join(ac_5, 
#                               annotated_fts_from_mum$solar %>% 
#                                 select(name,exactmass:median_cv))
# 
# 
# 
# # Join with original data
# annotated_fts_from_mum_final <-  annotated_fts_from_mum %>% 
#   modify(~select(.x, 
#                  -all_of(colnames(ac_5[-1])[
#                    colnames(ac_5[-1]) %in% colnames(annotated_fts_from_mum$solar)
#                  ] 
#                  )) %>% 
#            left_join(ac_5,.) %>% 
#            janitor::remove_empty())