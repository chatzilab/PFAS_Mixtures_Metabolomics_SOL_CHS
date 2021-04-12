# work flow for cleaning data
library(tidyverse)
source("!directories.R") # directories for file architecture
source("0.0.1_format_vars_funs.R") # load formats for variables


# load NAFLD outcome data
source("0.1_outcome_nafld.R", print.eval = TRUE)

# imputed exposome variables (dependent on previous outcome file)
source("0.2_exposome.R")

# load SNP data
source("0.3_genetic_snp.R")

# methylome
source("0.4_methylome.R")

# transcriptome
source("0.5_transcriptome.R")

# proteome
source("0.6_proteome.R")

# metabolome
source("0.7_metabolome.R")
