# Functions

# Create Functions for summarizing Data ---------------
transpose_ft <- function(ft) {
  dataout <- ft %>%
    mutate(name = str_c(mz, time, sep = "_")) %>%
    select(name, everything(), -mz, -time) %>%
    gather(file_name, val, 2:ncol(.)) %>%
    spread(name, val) %>% 
    mutate(file_name = str_remove(file_name, "_mz_xml"))
  
  return(dataout)
}

