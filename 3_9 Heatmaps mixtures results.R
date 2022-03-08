library(colorspace)
library(tidylog)
library(RColorBrewer)
library(dendextend)
library(gplots)

# Read in MWAS Beta Coefficients --------------------
# read in mwas beta coefs 
mwas_results_long  <- read_rds(
  file = fs::path(dir_results_mixtures, 
                  "SOL CHS all Mixtures MWAS results long hyper_g.rds"))

# # select all features p < 0.05
p_05 <- mwas_results_long %>%
  modify(., ~filter(., p_value < 0.05))
 
# met_any_pfas_05 <- p_05 %>%
#   modify(., ~filter(., feature %in% p_05$feature))

# Scale estimates 
# mwas_results_long <- mwas_results_long %>% 
#   modify(~.x %>% 
#            group_by(exposure) %>% 
#            mutate(estimate_scaled = estimate_beta)) #scale(estimate, center = FALSE)))  


# Pivot wider, change non-sig estimates to 0
mwas_results_wide <- mwas_results_long %>% 
  modify(~.x %>% 
           # mutate(.x, 
           #       estimate_scaled = if_else(p < 0.05, 
           #                                 estimate_scaled, 
           #                                 estimate_scaled)) %>%
           select(exposure, feature, estimate_beta, estimate_pip) %>% 
           pivot_wider(id_cols = feature, 
                       names_from = exposure, 
                       values_from = c(estimate_beta, estimate_pip)) %>% 
           select(-contains("mixture")))


# Read in Dougs annotations ------------------------------------------------
annotations <- read_rds(
  fs::path(dir_data, 
           "Common_Metabolites_annotated_SOLAR_CHS_V1_from_doug.RDS")) %>% 
  bind_rows(.id = "mode")


# Get list of all metabolites
temp <- annotations$refmet_name %>% 
  str_split("; ") %>% 
  unlist() %>% 
  unique()

# write_csv(data.frame(met = temp), 
#           fs::path(dir_results, "All Annotated Metabolites.csv"))


# Create name variable to match mwas results
annotations2 <- annotations %>% 
  mutate(name = str_c(mz, time, sep = "_"),
         across(refmet_name:number_studies_reported,
                ~str_split_fixed(., ";", n = Inf)[,1]))

# Annotate MWAS and select top n metabolites  ----------------------------------

# Join MWAS with annotations
mwas_results_annotated_wide <- mwas_results_wide %>% 
  modify(~inner_join(.x, annotations2, by = c("feature" = "name")))

# Filter top 50 metabolites
top_50 <- mwas_results_annotated_wide %>%
  modify(~ .x %>%
           rowwise() %>%
           mutate(sum_abs_est = sum(abs(c_across(
             estimate_beta_pfda:estimate_beta_pfos))),
                  # max_mix_effect = abs(`Mixture effect`),
                  max_est = max(abs(c_across(
                    estimate_beta_pfda:estimate_beta_pfos)))) %>%
           ungroup() %>%
           group_by(refmet_name) %>%
           filter( # max_mix_effect == abs(`Mixture effect`),
                    sum_abs_est == max(sum_abs_est),
                    sum_abs_est != 0) %>%
           filter(max_est == max(max_est)) %>%
           filter(mass_error_ppm == min(mass_error_ppm)) %>%
           ungroup() %>%
           select(sum_abs_est, everything()) %>%
           slice_max(order_by = max_est, 
                     n = 100))

# Save results for temp analysis of metabolites 
# write_csv(top_50$solar, fs::path(dir_results, "top 50.csv"))

# reorder super class results
# top_50 <- mwas_results_annotated_wide %>% 
top_50 <- top_50 %>%
  modify(~mutate(.x, 
                 super_class = fct_infreq(super_class)) %>%
           arrange(super_class))



# Change to matrix, extract PIPs
pips <- top_50 %>% 
  modify(
    ~select(.x, 
            refmet_name, 
            all_of(str_c("estimate_pip_",
                         unique(mwas_results_long$chs$exposure)[-7]))) %>% 
      rename_all(~str_remove(., "estimate_pip_")) %>% 
      column_to_rownames(var = "refmet_name") %>% 
      as.matrix())

# Change to matrix, extract PIPs
eff_est <- top_50 %>% 
  modify(
    ~select(.x, 
            refmet_name, 
            all_of(str_c("estimate_beta_",
                         unique(mwas_results_long$chs$exposure)[-7]))) %>% 
      rename_all(~str_remove(., "estimate_beta_")) %>% 
      column_to_rownames(var = "refmet_name") %>% 
      as.matrix())


# Hierarchical clustering of Exposures  ----------------------------------

# Run Correlation Matrix
exposure_cor_matrix <- map(eff_est, ~cor(.x, method = "spearman") )

# Get pearson dissimilarity & run clustering
temp_dend_exposure <- map(exposure_cor_matrix, 
                          ~as.dist(1 - .x) %>% 
                            hclust(method = "complete") %>%
                            as.dendrogram() )



cols_branches <- c("darkred", "forestgreen", "blue")
# Set the colors of 4 branches
dend1 <- temp_dend_exposure %>% 
  map(~color_branches(.x, k = 3, col = cols_branches))

col_labels <- dend1 %>% 
  map(~get_leaves_branches_col(.x))

col_labels <- map2(col_labels,dend1, 
                   ~.x[order(order.dendrogram(.y))])


# Clustering of metabolites based on chemical class ---------------------------

# Create string var
# top_50 <- top_50 %>% 
#   modify(~mutate(.x, 
#                  pathstring = str_c(super_class, 
#                                     main_class, 
#                                     sub_class,
#                                     refmet_name,
#                                     sep = "_") %>% tolower()))
# 
# # Subset solar only
# solar <- top_50$solar
# 
# # Create dataframe to calcualte distance
# txtdat <- full_join(tibble(ps = solar$pathstring, 
#                            s = solar$super_class, 
#                            m = solar$main_class, 
#                            sb = solar$sub_class, 
#                            met = solar$refmet_name), 
#                     tibble(ps = solar$pathstring, 
#                            s = solar$super_class, 
#                            m = solar$main_class, 
#                            sb = solar$sub_class, 
#                            met = solar$refmet_name), 
#                     by = character()) %>% 
#   mutate(across(everything(), as.character))
# 
# 
# # Calculate distance 
# txtdat2 <- txtdat %>% 
#   mutate(dst = 0,  
#          dst1 = case_when(s.x  != s.y ~ dst + .2,
#                           m.x  != m.y ~ dst + .1,
#                           sb.x != sb.y ~ dst + .05,
#                           met.x  != met.y ~ dst + .01,
#                           TRUE ~ dst)) %>% 
#   select(ps.x, ps.y, dst1) %>% 
#   pivot_wider(id_cols = ps.x, names_from = ps.y, values_from = dst1) %>% 
#   column_to_rownames(var = "ps.x") %>% 
#   as.matrix()
# 
# # Convert to distance matrix
# txtdat2[upper.tri(txtdat2)] <- c(0)
# diag(txtdat2) <- 0
# txtdat2 <- as.dist(txtdat2)  
# 
# 
# met_cluster <- as.dendrogram(hclust(txtdat2))



# Create color scheme for metabolite names --------------------------------
out <- data.frame(super_class = unique(top_50$solar$super_class), 
                  color = viridis::turbo(length(unique(top_50$solar$super_class))))


top_50 <- top_50 %>%
  modify(~.x %>% 
           left_join(out))

# Plot Heatmap ------------------------------------------------------------

# Set color Scheme
color.scheme <- rev(diverging_hcl(palette = "Cork",n = 100))

# Create Plot
tiff(file=fs::path(dir_reports, "Heatmaps", "SOLAR PFAS Mixtures hyper_g.tiff"),
     width=5, height=5, res=300, units="in")  
out <-
heatmap.2(eff_est$solar,
          Colv = dend1$solar,
          # Rowv = met_cluster,
          dendrogram = "both", #c("both","row","column","none"),
          revC = TRUE,  # rev col order of dendrogram 
          trace = "none", 
          
          # Change Color of row names
          colRow = top_50$solar$color,
          
          # color key + density info
          key = TRUE,
          density.info = "none",
          # densadj = 0.25,
          key.title = "Effect Estimate",
          key.xlab = "Effect Estimate",
          # ColSideColors = col_labels,
          
          margins = c(5,1), 
          
          col = color.scheme,
          srtCol = 20,
          
          # Plot Labels
          xlab = "PFAS",
          # ylab = "Metabolite",
          main = "SOLAR",
          labRow = "",
          # labRow = top_50$solar$refmet_name
          # labCol = FALSE
)
dev.off()

# CHS Heatmap -------------------------------------------------------------
tiff(file=fs::path(dir_reports, "Heatmaps", 
                   "CHS PFAS Mixtures Heatmap_hyper_g.tiff"),
     width=5, height=5, res=300, units="in")  
out <-
heatmap.2(eff_est$chs,# reorderfun = 
          Colv = dend1$chs,
          # Rowv = met_cluster,
          dendrogram = "both", #c("both","row","column","none"),
          revC = TRUE,  # rev col order of dendrogram 
          trace = "none", 
          
          
          # color key + density info
          key = TRUE,
          density.info = "none",
          densadj = 0.25,
          key.title = "Effect Estimate",
          key.xlab = "Effect Estimate",
          # ColSideColors = col_labels,
          
          # Change Color of row names
          colRow = top_50$chs$color,
          
          margins = c(5,1), 
          
          col = color.scheme,
          srtCol = 20,
          
          # Plot Labels
          xlab = "PFAS",
          # ylab = "Metabolite",
          main = "CHS",
          labRow = ""
          # labRow = top_50$chs$refmet_name
          # labCol = NULL
)
dev.off()




















