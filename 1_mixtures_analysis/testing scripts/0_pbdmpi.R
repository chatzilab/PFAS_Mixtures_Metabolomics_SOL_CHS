# script for metabolomics mixtures analysis
# Jesse Goodrich 111221

# Load Libraries 
library(rjags)
library(R2jags)
library(pbdMPI)
library(dplyr)

init()

#Set working directory
setwd("/project/dconti_624/Users/jagoodri/sol_chs")

# Load Data
load("SOLAR_mixtures_datasets.Rdata")

# Get metabolite names
# metabolites <- colnames(Y)
n_met <- ncol(Y)

# Model
model <- function(i){
  output <- BHRMA.g(X=X.obs,
                    Y=Y[,i],
                    U=U,
                    LOD=LOD,
                    profiles=profiles)
  
  output$metabolite = colnames(Y)[i]
  
  output <- output %>% 
    dplyr::select(metabolite, everything())
  # 
  return(output)
}


# Run model
coefs <- pbdLapply(1:n_met, model, pbd.mode = "spmd")

# Save results
comm.write.csv(coefs,  file = "results/test.csv")

# message(paste("SUCCESS from rank", comm.rank()))

finalize()