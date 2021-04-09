library(tidyverse)
# https://github.com/stephenslab/susieR
# devtools::install_github("stephenslab/susieR")
library(susieR)
source("5.1_susie_utils.R")


# SuSiE minimal example

  set.seed(1)
  n    <- 1000 # number of subjects in study
  p    <- 1000 # number of SNPs
  beta <- rep(0,p)
  beta[c(1,2,300,400)] <- 1
  X   <- matrix(rnorm(n*p), nrow=n, ncol=p) # SNP matrix
  y   <- X %*% beta + rnorm(n) # outcome matrix
  res <- susie(X,y,L=10)
  plot(coef(res), pch = 20)
  
  # Plot the ground-truth outcomes vs. the predicted outcomes:
  plot(y,predict(res), pch = 20)


# data formatting

  # get outcome and snp dataframe
  df_susie <- liver_enz %>% 
    dplyr::select(helixid, alt_n) %>% # only need liver injury outcome variable
    dplyr::mutate(alt_n = scale(log(alt_n))) %>% # log transform and scale
    inner_join(as.data.frame(snp_data$snp_imp), by = "helixid")

  # get outcome and SNPs as matrices
  outcome <- df_susie %>% select(alt_n) %>% as.matrix
  snp_cov <- df_susie %>% select(starts_with("rs")) %>% as.matrix()

  
# Simple regression summary statistics
sumstats <- univariate_regression(snp_cov, outcome)
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z")
  
  
# Fine-mapping with SuSiE
fit_liv_snp <- susie(snp_cov, outcome, # matrices formatted above
                     L = 122, # number of possible causal SNPs
                     # coverage = 0.9, # change the lower limit for credible sets
                     estimate_residual_variance = TRUE,
                     estimate_prior_variance = FALSE,
                     scaled_prior_variance = 0.1,
                     verbose = TRUE
)
  
  # plot coefficients
  plot(coef(fit_liv_snp), pch = 20)
  
  # susie plot
  susie_plot(fit_liv_snp, y = "PIP")
  
  # credible sets
  print(fit_liv_snp$sets)
  
  # Plot the ground-truth outcomes vs. the predicted outcomes:
  plot(outcome, predict(fit_liv_snp), pch = 10)
  
  