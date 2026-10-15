# Load Libraries 
library(data.table)
library(glmnet)
library(writexl)   
library(dplyr)


# ==========================
# 1) Load data
# ==========================
prot <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/prot_file_for_regression_filtered_Sample"))
prot[prot == ""] <- NA

rownames(prot) <- as.character(prot[[1]])
prot <- prot[, -1]  # proteins only (columns = proteins)

prot <- prot[,c(1:2923)]

# Covariates + complex traits
cov_pheno_filter_new <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/cov_pheno_filter.csv"))

# keep eid + the listed traits
trait_sub <- cov_pheno_filter_new %>%
  dplyr::select(
    eid,
    albumin, ALP, ALT, ApoA, ApoB,
    AST, bilirubin, BMI, BMR, body_fat_mass,
    Ca, cholesterol, Cr, CRP, cystatinC,
    diastolic_BP, fluid_intelligence, FVC, GGT, glucose,
    grip_strength_right, HbA1c, HDL, heart_rate,
    height, hip_circumference, IGF1, LDL, LPA,
    lymphocyte_count, MCH, monocyte_count, neuroticism, neutrophil_count,
    phosphate, platelet_count, RBC_count, reticulocyte_count, SHBG,
    systolic_BP, TG, total_protein, urate,
    urea, vitamin_D, waist_circumference, WBC_count, weight
  )


rownames(trait_sub) <- as.character(trait_sub$eid)
trait_sub <- trait_sub[,-1]

# Align individuals (by eid rownames)
common_ids <- intersect(rownames(prot), rownames(trait_sub))
common_ids <- sort(common_ids)
prot       <- prot[common_ids, , drop = FALSE]
trait_sub  <- trait_sub[common_ids, , drop = FALSE]
trait_sub <- as.data.frame(scale(trait_sub))


####### LASSO regression for each protein #########


# Initialize an empty data frame to store coefficients
total_coeffic <- data.frame()

# Loop over each protein 
for(prot_id in colnames(prot)){
  
  print(prot_id)
  
  X <- as.matrix(trait_sub)
  y <- prot[[prot_id]]
  
  y <- as.numeric(y)
  
  impute_mean <- function(x) {
  replace(x, is.na(x), mean(x, na.rm = TRUE))
	}

  trait_sub_imputed <- apply(trait_sub, 2, impute_mean)
  y_imputed <- ifelse(is.na(y), mean(y, na.rm = TRUE), y)

  # Fit the model 
  	cvfit <- cv.glmnet(trait_sub_imputed, y_imputed, 
                   alpha = 1, standardize = FALSE)

	# Extract optimal lambda and predict on the training data
	yhat <- predict(cvfit, newx = trait_sub_imputed, s = "lambda.min")

	yhat <- as.vector(yhat)
	
	# Compute R-squared
	r2_in <- 1 - sum((y_imputed - yhat)^2) / sum((y_imputed - mean(y_imputed))^2)

	# Compute R-squared oos
	r2_oof  <- 1 - min(cvfit$cvm) / var(y_imputed)

  # Extract the optimal lambda and corresponding coefficients
  best_lambda <- cvfit$lambda.min
  coef_fit <- coef(cvfit, s = "lambda.min")
  
  # Convert the coefficient matrix (a sparse dgCMatrix) to a data frame
  curr_coefficients <- data.frame(
    Trait = rownames(coef_fit),
    Coefficient = as.numeric(coef_fit)
  )
  
  # Add the protein ID to the results
  curr_coefficients$prot <- prot_id
  
  curr_coefficients$r2_in <- r2_in
  curr_coefficients$r2_oof <- r2_oof
  
  # filter out coefficients that are zero (except for the intercept)
  nonzero <- curr_coefficients$Coefficient != 0 | curr_coefficients$Trait == "(Intercept)"
  curr_coefficients <- curr_coefficients[nonzero, ]
  
  # Append current coefficients to the total results
  total_coeffic <- rbind(total_coeffic, curr_coefficients)
}



######################################## Save the coefficients #################################################

write_xlsx(total_coeffic, "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/03_lasso/data/coefficients_proteins_traits_lm_scaled_lasso_lambdamin.xlsx")  
