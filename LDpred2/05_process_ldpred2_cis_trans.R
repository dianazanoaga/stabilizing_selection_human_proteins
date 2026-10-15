library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(writexl)

rm(list= ls())

ldpred2_prot <- read_excel("~/selection_proteins/data_last/ldpred2_results_all_prot_combined_using_unstand_data_coef_final.xlsx")

ldpred2_prot <- ldpred2_prot[,-14] #removing fdr

ldpred2_prot[is.na(ldpred2_prot$alpha_est_1),] #0 

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
ldpred2_prot_corrected <- ldpred2_prot  # copy original data

for(i in 1:nrow(ldpred2_prot)) {
  corrected <- correct_quantiles(ldpred2_prot$alpha_est_1[i],
                                 ldpred2_prot$alpha_est_2[i],
                                 ldpred2_prot$alpha_est_3[i],
                                 lower_boundary = -5,
                                 upper_boundary = 1)
  ldpred2_prot_corrected$alpha_est_1[i] <- corrected[1]
  ldpred2_prot_corrected$alpha_est_2[i] <- corrected[2]
  ldpred2_prot_corrected$alpha_est_3[i] <- corrected[3]
}


# View the first few rows of the corrected data
head(ldpred2_prot_corrected)


# RECOMPUTE NOW SE AND PVAL

# Step 1: Calculate Standard Error (SE) from 95% CI
ldpred2_prot_corrected$se <- (ldpred2_prot_corrected$alpha_est_3 - ldpred2_prot_corrected$alpha_est_2) / (2 * 1.96)

# Step 2: Compute Z-score for alpha
ldpred2_prot_corrected$z  <- ldpred2_prot_corrected$alpha_est_1 / ldpred2_prot_corrected$se

# Step 3: Compute p-value (two-sided)
ldpred2_prot_corrected$p_value <- 2 * (1 - pnorm(abs(ldpred2_prot_corrected$z)))

logp <- pnorm(-abs(ldpred2_prot_corrected$z ), lower.tail = TRUE, log.p = TRUE) + log(2)  # log( 2 * Phi(-|z|) )
ldpred2_prot_corrected$p_value_log <- logp

neglog10_p <- -(log(2) + pnorm(abs(ldpred2_prot_corrected$z ), lower.tail = FALSE, log.p = TRUE)) / log(10)
ldpred2_prot_corrected$neglog10_p <- neglog10_p 

dim(ldpred2_prot_corrected[ldpred2_prot_corrected$neglog10_p > 1.3, ])
dim(ldpred2_prot_corrected[ldpred2_prot_corrected$p_value < 0.05, ])


p_a <- as.data.frame(read_excel("~/Desktop/codes/ldpred2_proteins_march/data/protein_annotation.xlsx"))
l_r <- ldpred2_prot_corrected

p_a_u <- p_a[,c(4,5,10)]
p_a_uu <- p_a_u[!duplicated(p_a_u$Assay),]

head(p_a)

colnames(l_r)[1] <- "Gene"

merge_lr <- merge(l_r, p_a_uu, by.x = "Gene", by.y = "Assay")

merge_lr_unique <- distinct(merge_lr, Gene, .keep_all = TRUE)

sum(is.na(merge_lr_unique$alpha_est_1)) # 4 are NAs

merge_lr_unique[is.na(merge_lr_unique$alpha_est_1),]

dim(merge_lr_unique[merge_lr_unique$alpha_est_1 < 0 & merge_lr_unique$p_value < 0.05/2922,])
dim(merge_lr_unique[merge_lr_unique$alpha_est_1 > 0 & merge_lr_unique$p_value < 0.05/2922,])

dim(merge_lr_unique[merge_lr_unique$alpha_est_1 == 0,])

write_xlsx(merge_lr_unique, "~/selection_proteins/data_last/process/ldpred2-auto_proteins_coefficients.xlsx")

