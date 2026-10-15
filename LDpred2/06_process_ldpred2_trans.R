library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(qtl2)
library(readr)
library(data.table)
library(ggplot2)
library(tidyr)
library(scales)

rm(list=ls())

### NEW LDPRED2 DATA - NO CIS #################################################

################################################################################


ldpred2_no_cis<- read_excel("ldpred2_results_all_prot_combined_using_unstand_data_coef_final_no_cis.xlsx")

lower_boundary = -5
upper_boundary = 1

correct_quantiles <- function(median, low, up, lower_boundary = -5, upper_boundary = 1, tol = 1e-6) {
  
  # If the median is at either boundary, return NA for all quantiles.
  if(abs(median - lower_boundary) < tol || abs(median - upper_boundary) < tol) {
    return(c(NA, NA, NA))
  }
  
  # If both quantiles are stuck at the boundaries, there's no information to correct.
  if(abs(low - lower_boundary) < tol && abs(up - upper_boundary) < tol) {
    return(c(NA, NA, NA))
  }
  
  # If only the lower quantile is stuck, correct it using the symmetry assumption.
  if(abs(low - lower_boundary) < tol && abs(up - upper_boundary) > tol) {
    low <- 2 * median - up
  }
  
  # If only the upper quantile is stuck, correct it using the symmetry assumption.
  if(abs(up - upper_boundary) < tol && abs(low - lower_boundary) > tol) {
    up <- 2 * median - low
  }
  
  return(c(median, low, up))
}


# Using a for loop to update each row
ldpred2_no_cis_corrected <- ldpred2_no_cis  # copy original data

for(i in 1:nrow(ldpred2_no_cis)) {
  corrected <- correct_quantiles(ldpred2_no_cis$alpha_est_1[i],
                                 ldpred2_no_cis$alpha_est_2[i],
                                 ldpred2_no_cis$alpha_est_3[i],
                                 lower_boundary = -5,
                                 upper_boundary = 1)
  ldpred2_no_cis_corrected$alpha_est_1[i] <- corrected[1]
  ldpred2_no_cis_corrected$alpha_est_2[i] <- corrected[2]
  ldpred2_no_cis_corrected$alpha_est_3[i] <- corrected[3]
}


# RECOMPUTE NOW SE AND PVAL

# Step 1: Calculate Standard Error (SE) from 95% CI
ldpred2_no_cis_corrected$se <- (ldpred2_no_cis_corrected$alpha_est_3 - ldpred2_no_cis_corrected$alpha_est_2) / (2 * 1.96)

# Step 2: Compute Z-score for alpha
ldpred2_no_cis_corrected$z  <- ldpred2_no_cis_corrected$alpha_est_1 / ldpred2_no_cis_corrected$se

# Step 3: Compute p-value (two-sided)
ldpred2_no_cis_corrected$p_value <- 2 * (1 - pnorm(abs(ldpred2_no_cis_corrected$z)))

logp <- pnorm(-abs(ldpred2_no_cis_corrected$z ), lower.tail = TRUE, log.p = TRUE) + log(2)  # log( 2 * Phi(-|z|) )
ldpred2_no_cis_corrected$p_value_log <- logp

neglog10_p <- -(log(2) + pnorm(abs(ldpred2_no_cis_corrected$z ), lower.tail = FALSE, log.p = TRUE)) / log(10)
ldpred2_no_cis_corrected$neglog10_p <- neglog10_p 


ldpred2_no_cis_corrected <- ldpred2_no_cis_corrected[!is.na(ldpred2_no_cis_corrected$alpha_est_1),]
ldpred2_no_cis_corrected <- ldpred2_no_cis_corrected[!duplicated(ldpred2_no_cis_corrected$alpha_est_1),]

write_xlsx(ldpred2_no_cis_corrected, "~/selection_proteins/data_last/process/ldpred2-auto_proteins_coefficients_no_cis.xlsx")
