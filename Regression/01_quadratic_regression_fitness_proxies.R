# Load necessary libraries

library(data.table)
library(dplyr)
library(stringr)
library(writexl)
library(readxl)


# Load file containing protein levels
 
prot <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/prot_file_for_regression_filtered_Sample"))
prot <- prot[,c(1:2924)]

prot[prot == ""] <- NA # eid, 2923 proteins 

# Column 1 is 'eid', so convert columns 2 through the last column (ncol(prot)) to numeric 
prot[, 2:ncol(prot)] <- lapply(prot[, 2:ncol(prot)], as.numeric)

# Load file containing covariates and complex traits 

all_cov <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_all_pheno_cov_filter_removed_and_add"))
all_cov <- all_cov[,-1]
all_cov_sampled <- all_cov[, c(1:45)] # eid, age, age2, sex, batch, PC1-40

# Load file containing longevity info

load("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/longevity.RData")  #df_both 
df_both <- df_both[,c(1,14)] #eid, avg long 

# Load file where there are covariates and composite health score(last column)
df_z2 <- read_excel("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_composite_regression_and_correlation_zscore.xlsx")

# select eid and health score
df_to_use <- df_z2[,c(1,46)] #eid, composite score 


# Load file containing nchildren info 

load("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/pheno_n.children_M_and_F.RData") #"data_n_child_fathered"  "data_n_live_births_fem"
head(data_n_child_fathered)
head(data_n_live_births_fem)

# Prepare pheno : n.children 
merge_child <- merge(data_n_child_fathered[,c(1,6)], data_n_live_births_fem[,c(1,6)], by ="eid")
merge_child$n.children <- rowSums(merge_child[,c(2,3)], na.rm =T)

merge_child$childlessness <- ifelse(merge_child$n.children == 0 , 1, ifelse(merge_child$n.children > 0, 0, NA))
merge_child <- merge_child[,c(1,4,5)]

merge_child <- merge_child[merge_child$n.children < 10,] # Filtering out individuals with 10 or more children
merge_child <- merge_child[,c(1,2)] #eid, nchildren 


# Load file containing other covariates 
ea <- load("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/pheno_education_att.RData") #"data_qualif" "dq_save" 

head(data_qualif)
head(dq_save)

data_qualif <- data_qualif[,c(1, 30)]
data_qualif$avg_abs <- ifelse(is.na(data_qualif$avg_abs), NA_real_, 
							  ifelse(data_qualif$avg_abs > 19, 1, 0))
							  
data_qualif$avg_abs <- as.numeric(data_qualif$avg_abs) 


all_cov_sampled_edu <-  merge(all_cov_sampled, data_qualif, by = "eid")

# Create the interaction term while handling NA values
all_cov_sampled_edu$edu_sex <- ifelse(is.na(all_cov_sampled_edu$sex) | is.na(all_cov_sampled_edu$avg_abs),
                                     NA_real_,
                                    all_cov_sampled_edu$sex * all_cov_sampled_edu$avg_abs)


reg_fun <- function(dataset_output, id, cov_id){

	# initialize a list 
	
	cov_pheno_list <- list()
	
	prot_to_use <- prot
	
	if(id == "longevity" || (id == "n.children" & cov_id == "no_other_cov")){
	
	print(paste0("Looking at phenotype: ", id ))
	
	cov_pheno_filter <- merge(all_cov_sampled,dataset_output, by = "eid")
	
	prot_to_use <- prot_to_use[prot_to_use$eid %in% cov_pheno_filter$eid,]
	cov_pheno_filter <- cov_pheno_filter[cov_pheno_filter$eid %in% prot_to_use$eid,]

	# order individuals 
	prot_to_use <- prot_to_use[order(prot_to_use$eid),]
	cov_pheno_filter <- cov_pheno_filter[order(cov_pheno_filter$eid),]
	
	if( nrow(cov_pheno_filter) == nrow(prot_to_use) ) {

		print("Same number of sample in cov_pheno and prot files.")

	} else {

	print("Not same number of sample in cov_pheno and prot files.")
	
	}
	
	temp_cov_pheno_filter <- cov_pheno_filter %>% dplyr::select(c(1:45, ncol(cov_pheno_filter)))
	
	# rename the phenotype column as "pheno"
	colnames(temp_cov_pheno_filter)[ncol(temp_cov_pheno_filter)] <- "pheno"
	
	
	# regress out covariates pheno ~ age + age^2 + sex + batch + PC1-40
	
	resid <- residuals(lm(pheno ~ . , data = temp_cov_pheno_filter[, -c(1)], na.action = na.exclude))

	resid <- scale(resid)
	
	} else if(id == "n.children" & cov_id == "edu_eduxsex_cov"){
	
	print(paste0("Looking at phenotype: ", id ))
	
	cov_pheno_filter <- merge(all_cov_sampled_edu, dataset_output, by = "eid")
	
	# filter protein and covariate data to keep only common eids
	prot_to_use <- prot_to_use[prot_to_use$eid %in% cov_pheno_filter$eid,]
	cov_pheno_filter <- cov_pheno_filter[cov_pheno_filter$eid %in% prot_to_use$eid,]

	# order individuals 
	prot_to_use <- prot_to_use[order(prot_to_use$eid),]
	cov_pheno_filter <- cov_pheno_filter[order(cov_pheno_filter$eid),]

	if(nrow(cov_pheno_filter) == nrow(prot_to_use) ) {
		print("Same number of sample in cov_pheno and prot files.")
	} else {
	    stop("Not same number of sample in cov_pheno and prot files. Check merging.")
	}
	
	temp_cov_pheno_filter <- cov_pheno_filter %>% dplyr::select(c(1:47, ncol(cov_pheno_filter)))
	
	# rename the phenotype column as "pheno"
	colnames(temp_cov_pheno_filter)[ncol(temp_cov_pheno_filter)] <- "pheno"
	
	# regress out covariates
	resid <- residuals(lm(pheno ~ . , data = temp_cov_pheno_filter[, -c(1)], na.action = na.exclude))

	resid <- scale(resid)
	
	} else if (id == "health_score" & cov_id == "no_other_cov"){
	print(paste0("Looking at phenotype: ", id ))
	
	cov_pheno_filter <- merge(all_cov_sampled, dataset_output, by = "eid")
	
	# filter protein and covariate data to keep only common eids
	prot_to_use <- prot_to_use[prot_to_use$eid %in% cov_pheno_filter$eid,]
	cov_pheno_filter <- cov_pheno_filter[cov_pheno_filter$eid %in% prot_to_use$eid,]

	# order individuals 
	prot_to_use <- prot_to_use[order(prot_to_use$eid),]
	cov_pheno_filter <- cov_pheno_filter[order(cov_pheno_filter$eid),]

	if(nrow(cov_pheno_filter) == nrow(prot_to_use) ) {
		print("Same number of sample in cov_pheno and prot files.")
	} else {
	    stop("Not same number of sample in cov_pheno and prot files. Check merging.")
	}
	
	health_score_col <- cov_pheno_filter[, ncol(cov_pheno_filter)]
    
    resid <- scale(-health_score_col)
	
	}
	
for(gp in c(2:ncol(prot_to_use))){

	cov_pheno_each <- data.frame(Pheno = id)

	# create a temporary dt containing only eid, and the respective prot column
	temp_gene_prs_filter <- data.frame(eid = prot_to_use$eid)
	mom <- prot_to_use %>% dplyr::select(c(gp))

	# extract gene name 
	gene_name <- colnames(prot_to_use)[c(gp)]
	
	print(paste0("Looking at gene: ", gene_name, "; which is column number: ", gp))
	
	# rename 
	colnames(mom) <- "protein"
	temp_gene_prs_filter$protein <- mom$protein
	temp_gene_prs_filter$protein <- as.numeric(temp_gene_prs_filter$protein)
	temp_gene_prs_filter$protein <- scale(temp_gene_prs_filter$protein) 

	
	# add a new column when storing residuals from the first regression
	temp_gene_prs_filter$resid <- resid
	
	lm_prs <- lm(resid ~ protein + I(protein^2), temp_gene_prs_filter[, -c(1)], na.action = na.exclude)
	
	lm_summary <- summary(lm_prs)

	# store in cov_pheno_each the information from regression for the respective gene 
	cov_pheno_each$Gene <- gene_name

	cov_pheno_each$Coefficients_intercept <- lm_summary$coefficients[1, "Estimate"]
  	cov_pheno_each$StdErrors_intercept <-  lm_summary$coefficients[1, "Std. Error"]
  	cov_pheno_each$TValues_intercept <-  lm_summary$coefficients[1, "t value"]
  	cov_pheno_each$PValues_intercept <-  lm_summary$coefficients[1, "Pr(>|t|)"]
	
	cov_pheno_each$Coefficients_lin <- lm_summary$coefficients[2, "Estimate"]
  	cov_pheno_each$StdErrors_lin <-  lm_summary$coefficients[2, "Std. Error"]
  	cov_pheno_each$TValues_lin <-  lm_summary$coefficients[2, "t value"]
  	cov_pheno_each$PValues_lin <-  lm_summary$coefficients[2, "Pr(>|t|)"]


	cov_pheno_each$Coefficients_quad <-  lm_summary$coefficients[3, "Estimate"]
  	cov_pheno_each$StdErrors_quad <-  lm_summary$coefficients[3, "Std. Error"]
  	cov_pheno_each$TValues_quad <- lm_summary$coefficients[3, "t value"]
  	cov_pheno_each$PValues_quad <-  lm_summary$coefficients[3, "Pr(>|t|)"]

	cov_pheno_list[[gp-1]] <- cov_pheno_each
	

	}
	
	# bind the information related to the single gene to all the other genes already analysed for that phenotype
	cov_pheno_all <- do.call(rbind, cov_pheno_list)
	
	write_xlsx(cov_pheno_all, paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/01_quadratic_regression/data/protein_regression_", id, "_", cov_id, "_.xlsx"))
	
	}
	
	
reg_fun(df_both, "longevity", "no_other_cov")
reg_fun(df_to_use, "health_score", "no_other_cov")
reg_fun(merge_child, "n.children", "no_other_cov")
reg_fun(merge_child, "n.children", "edu_eduxsex_cov")




