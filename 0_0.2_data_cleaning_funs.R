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




## Quantile function
qntle_fxn <- function(x, pct){ 
  pct <- quantile(x, pct) %>% signif(., digits = 3)
  
  pct2 <- if_else(pct == 0, 
                  "<LOD", 
                  as.character(pct))
  return(pct2)
}



####################
avg_sd_fxn <- function(x, n.digits = 1, include.n = FALSE) {
  out <- paste(round(mean(x, na.rm = T), n.digits) %>% 
                 formatC(. , format="f", digits=n.digits),
               round(sd(x, na.rm = T), n.digits) %>% 
                 formatC(. , format="f", digits=n.digits), sep = " ± ") 
  if(include.n == TRUE) {
    out <- paste(out, 
                 " (",
                 paste(length(which(!is.na(x))), ")", sep = ""), 
                 sep = "") }
  return(out)
}

fungm <- function(x, digits = 3){
  gm <- Gmean(x, conf.level = .95, na.rm = TRUE) %>% 
    signif(., digits) %>% 
    formatC(. , format="g", digits=digits)
  
  gm2 <- paste(gm[1], " [", gm[2], ", ", gm[3], "]", sep = "")
  # gm3 <- if_else(is.finite(Gmean(x)), gm2, "*")
  return(gm2)
}

fungsd <- function(x, digits = 1){
  gm <- Gmean(x, na.rm = TRUE) %>% 
    round(., digits) %>% 
    formatC(. , format="f", digits=digits)
  gs <- Gsd(x, na.rm = TRUE) %>% 
    round(., digits) %>% 
    formatC(. , format="f", digits=digits)
  return(paste(gm[1], " (", gs,")", sep = ""))
}

npct <- function(x, level_of_interest, digits = 1){ 
  number_meeting_criteria = sum(x == level_of_interest)
  pct_of_pop = 100*round(number_meeting_criteria/length(x),
                         digits = digits)
  return(paste(number_meeting_criteria, 
               " (", 
               pct_of_pop,
               "%)", sep = ""))
}


