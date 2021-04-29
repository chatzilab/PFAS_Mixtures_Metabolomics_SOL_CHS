library(RColorBrewer)
library(dendextend)
library(gplots)

# Get Correlation Matrix
data <- solar_c18

met_cor_plot <- function(data, file.location){ 
  # Calculate Correlation Matrix
  temp_cor_matrix <- data %>% 
    select(-id) %>% 
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
  # clusters.df <- data.frame(metabolite = names(clusters), cluster = clusters)
  
  # Set color Scheme
  color.scheme <- rev(brewer.pal(10,"RdBu")) # generate the color scheme to use
  
  # Create Plot
  jpeg(file=file.location)
  out <- heatmap.2(temp_cor_matrix, 
                   Rowv = ladderize(temp_dend), 
                   Colv = ladderize(temp_dend), 
                   dendrogram = "both", 
                   revC = TRUE,  # rev column order of dendrogram so conforms to natural representation
                   trace = "none", 
                   density.info = "none",
                   col = color.scheme, 
                   key = FALSE,
                   labRow = FALSE, labCol = FALSE)
  dev.off()
  
}




met_cor_plot(solar_c18, 
             file.location = here::here("figures", 
                                        "solar_c18_correlation_matrix.jpeg"))
met_cor_plot(solar_hilic, 
             file.location = 
               here::here("figures", 
                          "solar_hilic_correlation_matrix.jpeg"))

# CHS
met_cor_plot(chs_c18, 
             file.location = here::here("figures", 
                                        "chs_c18_correlation_matrix.jpeg"))
met_cor_plot(chs_hilic, 
             file.location = 
               here::here("figures", 
                          "chs_hilic_correlation_matrix.jpeg"))
