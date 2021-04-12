# 0.0.2_data_cleaning_funs.R

vzdt <- function(data){
  return( lst("Dimensions" = dim(data),
              "Head" = head(data[,1:7])))
}


transpose_ft <- function(feature_table, pct_metab_na_allowed = 0.2, mzrt = TRUE){
  if(mzrt == TRUE){
    feature_table <- feature_table %>% 
      mutate(name = str_c(mz, time, sep = "_")) %>% 
      select(name, everything(), -mz, -time)
  }
  feature_table <- as.data.frame(feature_table)
  
  ft_t = setNames(data.frame(t(feature_table[,-1])), feature_table[,1])
  # Specify Fraction Metabolite Missing Allowed
  mv <- pct_metab_na_allowed
  # Drop metabolites with more missing than allowed
  ft_t <- ft_t[, which(colMeans(ft_t == 0) < mv)]
  ### Rownames to variable called "sample"-----------
  ft_t$sample = rownames(ft_t)
  ft_t <- select(ft_t, sample, everything())
  rownames(ft_t) <- NULL
  return(ft_t)
}


untranspose_ft <- function(feature_table_t){
  feature_table_t <- as.data.frame(feature_table_t)
  ft = setNames(data.frame(t(feature_table_t[,-1])), feature_table_t[,1])
  ### Rownames to variable called "names"-----------
  ft$name = rownames(ft)
  ft <- select(ft, name, everything())
  rownames(ft) <- NULL
  return(ft)
}


# Cuberoot function
cuberoot <- function(x)sign(x)*abs(x)^(1/3)

# Dicotomize PFAS function
pfas_dicot <- function(x){
  qnt <- quantile(x, probs = 2/3)
  return(if_else(x > qnt, "High levels", "Low levels"))
}
