# All Metabolites
library(RColorBrewer)
library(dendextend)
library(gplots)

# Read in chs mz/cpd key and mwas results

composite_exposure_vars = c("ocs", "pbde_num_detect", "pcb_num_detect")
continuous_exposure_vars = c("hexachlorobenzene_impute", "dde_impute",
                             "pbde_47_impute", "pbde_154_impute")

table(chs_single_matches_w_mwas$exposure)
chs_single_matches_w_mwas <- readRDS(here::here("Temporary results", 
                                                  "chs_mz_cpd_pathway_key_w_mwas.rds"))

top_50 <- chs_single_matches_w_mwas %>% 
  group_by(metabolite) %>% 
  filter(p_value == min(p_value)) %>% 
  ungroup() %>% 
  slice_max(n = 50, order_by = p_value)


chs_mwas_data = chs_single_matches_w_mwas %>% 
  # filter(str_detect(exposure, pattern = "num_detect") | 
  #           str_detect(exposure, "dde") | 
  #           str_detect(exposure, "hexachlorobenzene")) %>%
  # filter(!(exposure %in% composite_exposure_vars)) %>%
  filter(exposure %in% continuous_exposure_vars) %>%
  select(metabolite, exposure, beta) %>% 
  pivot_wider(id_cols = "metabolite", names_from = exposure, values_from = beta) %>% 
  as.data.frame() %>% 
  filter(metabolite %in% top_50$metabolite) %>% 
  rename_all(~str_remove(., "_impute") %>% toupper(.)) %>% 
  rename("HCB" = "HEXACHLOROBENZENE")


rownames(chs_mwas_data) <- chs_mwas_data$metabolite 

chs_mwas_matrix <- chs_mwas_data %>% 
  select(-METABOLITE) %>% 
  as.matrix()


## chs HCB new analysis-------------------
chs_c18_metab <- chs_single_matches_w_mwas %>% 
  group_by(metabolite) %>% 
  slice_head() %>% 
  filter(mode == "negative", 
         metabolite %in% top_50$metabolite) %>%
  ungroup()
chs_hilic_metab <- chs_single_matches_w_mwas %>% 
  group_by(metabolite) %>% 
  slice_head() %>% 
  filter(mode == "positive", 
         metabolite %in% top_50$metabolite) %>%
  ungroup()


chs_data <- left_join(chs_c18 %>% select(id, all_of(chs_c18_metab$metabolite)), 
                        chs_hilic %>% select(id, all_of(chs_hilic_metab$metabolite))) %>% 
  select(-id) %>%
  scale()


###################
## Exposures 
###################
imputed_ocs <- c("hexachlorobenzene_impute",
                 "dde_impute",
                 "pbde_154_impute",
                 "pbde_47_impute")


# First, select id and exposures and pivot longer
chs_exposure <- chs_exposure_outcome %>% 
  filter(!is.na(dde)) %>%
  # select(contains("num_detect"), 
  #        contains("dde_impute"), 
  #        contains("hexachlorobenzene_impute")) %>%
  # dplyr::select(all_of(exposures),
  #               -c(ocs, contains("num_detect"))) %>%
  select(contains("impute")) %>%
  rename_all(~.x %>% 
               str_remove("_impute") %>% 
               str_replace("_ngml_detect", "*") %>% 
               str_remove("_ngml") %>% 
               str_replace("_", " ") %>% 
               toupper()) %>%
  rename("HCB"= "HEXACHLOROBENZENE") %>% 
  select("HCB", "DDE",
         "PBDE 154", "PBDE 47", 
         # "PBDE 85*", "PBDE 100*",
         # "PBDE 153*", 
         # "PCB 118*", "PCB 138*", "PCB 153*",
         # "PCB 180*"
  )  %>%
  as.data.frame() %>% 
  as.matrix()


# Run Correlation Matrix
chs_exposure_cor_matrix <- polycor::hetcor(chs_exposure)

# Get distance 
temp_dist <- as.dist(1 - chs_exposure_cor_matrix$correlations)

# Run Clustering algorithm
temp_tree <- hclust(temp_dist, method="complete")
# plot(temp_tree, cex=0.000000001)

# create dendrogram object
temp_dend_exposure <- as.dendrogram(temp_tree) 

# # Get 10 major clusters
# clusters <- cutree(temp_dend_exposure, k=10)

# Plot 
# plot(color_branches(temp_dend, k=10),leaflab="none")
# clusters.df_exposure <- data.frame(metabolite = names(clusters), cluster = clusters)




# Metabolites -------------------------------------------------------------

# met_cor_plot <- function(data, file.location){ 
# Calculate Correlation Matrix
temp_cor_matrix <- chs_data %>%
  cor(use="pairwise.complete.obs")

# Get distance 
temp_dist <- as.dist(1 - temp_cor_matrix)

# Run Clustering algorhitm
temp_tree <- hclust(temp_dist, method="complete")
# plot(temp_tree, cex=0.000000001)

# create dendrogram object
temp_dend <- as.dendrogram(temp_tree) 

# # Get 10 major clusters
# clusters <- cutree(temp_dend, k=10)

# Plot 
# plot(color_branches(temp_dend, k=10),leaflab="none")
clusters.df <- data.frame(metabolite = names(clusters), cluster = clusters)

# Set color Scheme
color.scheme <- rev(brewer.pal(10,"RdBu")) # generate the color scheme to use

# Create Plot
# jpeg(file="chs_heatmap.jpg")
# out <-
heatmap.2(chs_mwas_matrix,
          Rowv = ladderize(temp_dend),
          Colv = ladderize(temp_dend_exposure),
          dendrogram = "both", #c("both","row","column","none"),
          revC = TRUE,  # rev column order of dendrogram so conforms to natural representation
          trace = "none", 
          density.info = "none",
          col = color.scheme, 
          srtCol = 20,
          key = TRUE,
          labRow = FALSE, 
)


# dev.off()