# Install MetaboanalystR 

# Install Dependencies ----------------------------------------------------
install.packages("fitdistrplus")

metanr_packages <- function(){
  
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()


# SSPA is not available in the current version of Bioconductor.
# (https://bioconductor.org/about/removed-packages/)
# I had to download the source packates off of bioconductor, but you also need:
BiocManager::install("qvalue")

# Test to see if it can be installed:
library(SSPA)

# Run metanr_packages again to see if anything is still missing: 
metanr_packages()

# Also need OptiLCMS: 
devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = FALSE, build_manual =TRUE)

# Step 2: Install MetaboAnalystR without documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)


# Check to see if it will load:
library(MetaboAnalystR)






# Install xmwas -----------------------------------------------------------

BiocManager::install(c("GO.db","graph","RBGL","impute","preprocessCore","mixOmics"),
                     dependencies=TRUE)

library(devtools); install_github("kuppal2/xMWAS")


