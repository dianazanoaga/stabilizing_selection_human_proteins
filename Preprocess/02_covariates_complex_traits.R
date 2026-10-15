# Load necessary libraries

library(data.table)
library(dplyr)
library(stringr)
library(readxl)

########## CREATE A DT CONTAINING ALL COVARIATES AND PHENOTYPES INFO 

########## 1. Create dt containing all the phenotypes info 

# Read phenotypes and respective UKBiobank ID
#pheno_list <- read.table("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/eqtl_g/03_Regression/data/pheno_list_chiara.txt", head=T)

pheno_list <- read.table("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/pheno_list_chiara.txt", head=T)

# Modify the path adding "old" 
modified_path <- sub("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/", "/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/", pheno_list$File)
pheno_list$File <- modified_path

corresp <- as.data.frame(read_excel("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/correspondence_uk_pheno.xlsx"))
        
# Extract numeric portion from file column
value  <- as.numeric(gsub("\\D", "", pheno_list$File))
value_new <- vector()

for(i in c(1:length(value))){
	
	n <- corresp[(corresp$original == value[i]) | (corresp$"refresh 1" == value[i]),]$"refresh 2"
	value_new[i] <- n[!is.na(n)]
	pheno_list$File_diana[i] <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/pheno/ukb", value_new[i], ".csv")
	i = i+1

	} 

ID <- pheno_list[,1]

pheno_list$Col_diana <- NA

for(id in ID){

	print(id)

	df <- fread(pheno_list[pheno_list$FieldID == id,]$File_diana, nrows=0)
  		
	id_occurrences <- which(grepl(paste0("\\b", id ,"-", "\\b"), colnames(df)))

	if(length(id_occurrences) > 0) 	{
	
			pheno_list[pheno_list$FieldID == id,]$Col_diana <- paste0(id_occurrences, collapse = "_")
	
		}		

	}


fwrite(pheno_list, "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/pheno_list_update.txt", col.names = T, row.names = F, quote = F, sep = "\t")

########################################################################

########################################################################

# Load file containing all eid samples
samples_pheno <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/plink/_001_ukb_cal_chr1_v2.fam", header = F, select = c(1), col.names = c("eid")))
print(paste0("Total number of samples eid :", dim(samples_pheno)[1]))

# Loop over all phenotypes to extract the data and calculate the average over different measurement instances, using pheno_list from Chiara 

for(id in pheno_list$FieldID) {
	
	# Define the phenotype
	pheno <- pheno_list[pheno_list$FieldID == id,]$"Pheno" 
	ID <- pheno_list[pheno_list$FieldID == id, ]$"FieldID"
	file <- pheno_list[pheno_list$FieldID == id,]$"File_diana" 
	col <- pheno_list[pheno_list$FieldID == id,]$"Col_diana" 
	print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))

	# Extract columns corresponding to the defined phenotype
	temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))

	# Average (if more than two column present, one is eid and the other the pheno)
	if (ncol(temp_pheno) > 2) {
	temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))
	}

	# Merge
	colnames(temp_pheno) <- c("eid", pheno)

	# I do a full join, so I keep all the columns, because the eid changes based on the phenotypes, and put NA where there is no information  
	samples_pheno <- full_join(samples_pheno, temp_pheno, by = "eid")

} 

print(paste0("Dimensions of phenotype table: ", ncol(samples_pheno), " pheno x ", nrow(samples_pheno), " eids"))

# Replace entries by NA 
samples_pheno[samples_pheno < 0] <- NA

# Exclude sex specific traits 
samples_pheno <- samples_pheno[, !names(samples_pheno) %in% c("facial_hair", "balding", "menarche", "menopause","birth_weight_first_child")]

# Save the file 
fwrite(samples_pheno, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_all_pheno")


########### 2. Create dt containing all the phenotypes info + the other covariates info

# Retrieve all the covariates - load file about covariates created by Chiara 
covariates <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/covariates.txt")

# Do not considering 5 column which is "array", but only considering "batch" because is giving more information 
cov <- covariates[,-5]
colnames(cov)[1] <- "eid"

# Merge to have a dataframe with all the covariates and phenotypes information
cov_pheno <- as.data.frame(full_join(cov, samples_pheno))

# Load samples eid after filtering 
sample_filtered <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_filtering_sample")
eid_filter<- sample_filtered[,1]

# Take only filtered samples from cov_pheno
cov_pheno_filter <- cov_pheno[cov_pheno$eid %in% eid_filter$eid,]
print(paste0( "After taking only filtered samples of cov_pheno file, there are : ", nrow(cov_pheno_filter), " rows"))

# Save the file 
fwrite(cov_pheno_filter, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_all_pheno_cov_filter")

print("End!")

########## add two new phenotypes ######################

before_pheno <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_all_pheno_cov_filter"))

# create a new column for body fat mass percentage 
before_pheno$"bdfp" <- rep(NA, dim(before_pheno)[1])

before_pheno$"bdfp" <- ifelse(before_pheno$sex == 1, 1.20 * before_pheno$BMI + 0.23 * before_pheno$age - 16.2, 1.20 * before_pheno$BMI + 0.23 * before_pheno$age - 5.4)

# create a new column for eGFR
before_pheno$"gfr" <- rep(NA, dim(before_pheno)[1])

before_pheno$"gfr" <- ifelse((before_pheno$sex == 1 & before_pheno$Cr <= (0.9 * 0.0113)), (141 * (before_pheno$Cr * 0.0113 / 0.9)^(-0.411)) * ((0.993)^(before_pheno$age)),
                          ifelse((before_pheno$sex == 1 & before_pheno$Cr > (0.9 * 0.0113)), (141 * (before_pheno$Cr * 0.0113 / 0.9)^(-1.209)) * ((0.993)^(before_pheno$age)),
                                 ifelse((before_pheno$sex == 2 & before_pheno$Cr <= (0.7 * 0.0113)), (144 * (before_pheno$Cr * 0.0113/0.7)^(-0.329)) * ((0.993)^(before_pheno$age)), ifelse((before_pheno$sex == 2 & before_pheno$Cr > (0.7 * 0.0113)), (144 * (before_pheno$Cr * 0.0113/ 0.7)^(-1.209)) * ((0.993)^(before_pheno$age)), NA))))

before_pheno <- as.data.frame(before_pheno)
after_pheno <- before_pheno[, !(colnames(before_pheno) %in% c("Cr", "body_fat_mass"))]

write.table(after_pheno, "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_all_pheno_cov_filter_removed_and_add")


