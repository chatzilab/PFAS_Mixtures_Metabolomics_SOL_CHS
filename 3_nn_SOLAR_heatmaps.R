# All Metabolites
library(RColorBrewer)
library(dendextend)
library(gplots)


# Read in SOLAR mz/cpd key and mwas results

composite_exposure_vars = c("ocs", "pbde_num_detect", "pcb_num_detect")
continuous_exposure_vars = c("hexachlorobenzene_impute", "dde_impute",
                             "pbde_47_impute", "pbde_154_impute")


solar_single_matches_w_mwas <- readRDS(here::here("Temporary results", 
                                                  "solar_mz_cpd_pathway_key_w_mwas.rds"))

# Select top 50 features
top_50 <- solar_single_matches_w_mwas %>% 
   group_by(metabolite) %>% 
   filter(p_value == min(p_value)) %>% 
   ungroup() %>% 
   slice_max(n = 50, order_by = p_value)

# Pivot mwas data wider
solar_mwas_data_wide = solar_single_matches_w_mwas %>% 
   # filter(str_detect(exposure, pattern = "num_detect") | 
   #           str_detect(exposure, "dde") | 
   #           str_detect(exposure, "hexachlorobenzene")) %>%
   # filter(!(exposure %in% composite_exposure_vars)) %>%
   filter(exposure %in% continuous_exposure_vars) %>%
   select(metabolite, exposure2, beta) %>% 
   pivot_wider(id_cols = "metabolite", names_from = exposure2, values_from = beta) %>% 
   as.data.frame() %>% 
   filter(metabolite %in% top_50$metabolite)

# Data frame to matrix
rownames(solar_mwas_data_wide) <- solar_mwas_data_wide$metabolite 

# Create Dataset for the heatmap 
solar_mwas_matrix <- solar_mwas_data_wide %>% 
   select(-metabolite) %>% 
   as.matrix()


# Hierarchical clustering of Exposures  ----------------------------------


# Get matrix of exposures 
solar_exposure <- solar_exposure_outcome %>% 
   filter(!is.na(dde)) %>%
   # select(contains("num_detect"), 
   #        contains("dde_impute"), 
   #        contains("hexachlorobenzene_impute")) %>%
   dplyr::select(all_of(exposures),
                 -c(ocs, contains("num_detect"))) %>%
   # select(contains("impute")) %>%
   rename_all(~.x %>% 
                 str_remove("_impute") %>% 
                 str_replace("_ngml_detect", "*") %>% 
                 str_remove("_ngml") %>% 
                 str_replace("_", " ") %>% 
                 toupper()) %>%
   rename("HCB"= "HEXACHLOROBENZENE") %>% 
   select("HCB", "DDE",
          "PBDE 154", "PBDE 47", 
          "PBDE 85*", "PBDE 100*",
          "PBDE 153*",
          "PCB 118*", "PCB 138*", "PCB 153*",
          "PCB 180*"
          ) %>% 
   mutate(across(contains("*"), 
                 ~as.ordered(.)))
   

# Run Correlation Matrix
solar_exposure_cor_matrix <- polycor::hetcor(solar_exposure)

# Get pearson dissimilarity 
temp_dist <- as.dist(1 - solar_exposure_cor_matrix$correlations)

# Run Clustering algorhitm
temp_tree <- hclust(temp_dist, method = "complete")
# plot(temp_tree, cex=1 )

# create dendrogram object
temp_dend_exposure <- as.dendrogram(temp_tree) 

# # Get 10 major clusters
# clusters <- cutree(temp_dend_exposure, k=10)

# Plot 
# plot(color_branches(temp_dend, k=10),leaflab="none")
# clusters.df_exposure <- data.frame(metabolite = names(clusters), cluster = clusters)




# Metabolites -------------------------------------------------------------

# Hierarchical clustering of metabolites ----------------------------------

solar_c18_metab <- solar_single_matches_w_mwas %>% 
   group_by(metabolite) %>% 
   slice_head() %>% 
   filter(mode == "negative", 
          metabolite %in% top_50$metabolite) %>%
   ungroup()
solar_hilic_metab <- solar_single_matches_w_mwas %>% 
   group_by(metabolite) %>% 
   slice_head() %>% 
   filter(mode == "positive", 
          metabolite %in% top_50$metabolite) %>%
   ungroup()


solar_data <- left_join(solar_c18 %>% select(id, all_of(solar_c18_metab$metabolite)), 
                        solar_hilic %>% select(id, all_of(solar_hilic_metab$metabolite))) %>% 
   select(-id) %>%
   scale()

# solar_data <- solar_hilic %>% 
#    select(id, all_of(solar_hilic_metab$metabolite))%>% 
#    select(-id) %>%
#    scale()

# met_cor_plot <- function(data, file.location){ 
# Calculate Correlation Matrix
temp_cor_matrix <- solar_data %>%
   cor(use="pairwise.complete.obs")

# Get distance 
temp_dist <- as.dist(1 - temp_cor_matrix)

# Run Clustering algorhitm
temp_tree <- hclust(temp_dist, method="complete")
# plot(temp_tree, cex=0.000000001)

# create dendrogram object
temp_dend_solar <- as.dendrogram(temp_tree) 

# # Get 10 major clusters
# clusters <- cutree(temp_dend, k=10)

# Plot 
# plot(color_branches(temp_dend, k=10),leaflab="none")
clusters.df <- data.frame(metabolite = names(clusters), cluster = clusters)

# Set color Scheme
color.scheme <- rev(brewer.pal(10,"RdBu")) # generate the color scheme to use

# Create Plot
# jpeg(file=file.location)
# out <-
heatmap.2(solar_mwas_matrix,# reorderfun = 
          # Rowv = ladderize(temp_dend_solar), 
          # Colv = FALSE,
          dendrogram = "both", #c("both","row","column","none"),
          revC = TRUE,  # rev column order of dendrogram so conforms to natural representation
          trace = "none", 
          density.info = "none",
          col = color.scheme, 
          srtCol = 20,
          key = TRUE,
          labRow = FALSE, 
          # labCol = TRUE
)
# dev.off()

# }

# str(solar_mwas_matrix)

