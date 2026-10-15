# Load necessary libraries

library(data.table)
library(dplyr)
library(stringr)

########## Load data

# Load file containing as columns all the covariates that need to be included in the model and also all the phenotypes to be analysed, 
# and as rows all the sample eids, already filtered.  

cov_pheno_filter <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_all_pheno_cov_filter"))

print(paste0("Number of rows: ", nrow(cov_pheno_filter), "; Number of columns: ", ncol(cov_pheno_filter)))

eid <- cov_pheno_filter[[1]]

covariates <- names(cov_pheno_filter)[2:45]
traits     <- names(cov_pheno_filter)[46:98]

# Initialize output df with eid
cov_pheno_filter_new <- cov_pheno_filter[,c(1:45)]

# Loop over traits 
for(tr in traits) {

  print(tr)
  
  # Build formula: trait ~ cov1 + cov2 + ... + covN
  f <- as.formula(paste(tr, "~", paste(covariates, collapse = " + ")))
  
  # Fit with na.exclude so residuals() returns a vector the same length as data
  fit <- lm(f, data = cov_pheno_filter, na.action = na.exclude)
  
  # Extract residuals
  cov_pheno_filter_new[[tr]] <- residuals(fit)
  
}

write.csv(cov_pheno_filter_new, "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/cov_pheno_filter.csv", row.names = FALSE)







