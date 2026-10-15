# load necessary library 

library(writexl)
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(readxl)
library(data.table)
library(magrittr)
library(dplyr)
library(arrow)

# List all files path - contain results of ldpred2 for protein 

all_prot_files <- list.files("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/proteins/data/unstand_prot_coef_no_cis/")

dt <- data.frame()

for(f in all_prot_files){

	dt_i <- as.data.frame(read_excel(paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/proteins/data/unstand_prot_coef_no_cis/", f)))
	
	dt <- rbind(dt, dt_i)
	
	}

dt[,c(2:10)] <- lapply(dt[,c(2:10)], as.numeric)

# Step 1: Calculate Standard Error (SE) from 95% CI
dt$se <- (dt$alpha_est_3 - dt$alpha_est_2) / (2 * 1.96)

# Step 2: Compute Z-score for alpha
dt$z  <- dt$alpha_est_1 / dt$se

# Step 3: Compute p-value (two-sided)
dt$p_value <- 2 * (1 - pnorm(abs(dt$z)))

dt$p_value_fdr <- p.adjust(dt$p_value, method = "fdr")

write_xlsx(dt, "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/proteins/data/ldpred2_results_all_prot_combined_using_unstand_data_coef_final_no_cis.xlsx")





