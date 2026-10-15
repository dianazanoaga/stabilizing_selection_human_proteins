# load necessary libraries

library(data.table)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(rbgen)
library(genio)
library(Matrix)
library(mbend)
library(readxl)


# Set the directory containing the CSV files
folder_path <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/pheno"

# Get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

corresp <- as.data.frame(read_excel("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/correspondence_uk_pheno.xlsx"))
    
# Specify the column names you're interested in
column_of_interest <- c("6138", "845") # respectively : Qualifications, Age completed full time education


# Initialize a data frame to store the results

new_pheno <- data.frame(id = character(),
                        file = character(),
                        col = character(),
                        stringsAsFactors = FALSE)


# Loop over each column of interest
for (i in 1:length(column_of_interest)) {
  
  print(i)
  # Loop over each CSV file
  for (file in csv_files) {
    
    # Read only the header of the CSV file
    data <- fread(file, nrows = 0, header = TRUE)
    
    # Get the column names
    column_names <- colnames(data)
    
    matching_indices <- grep(paste0("^", column_of_interest[i], "-"), column_names)
    
    if (length(matching_indices) > 0) {
      
	  # Save the column of interest, file, and matching indices to the new_pheno data frame
	  mom_pheno <- data.frame(id = column_of_interest[i],
                              file = file,
                              col = paste(column_names[matching_indices], collapse = "_"),
                              stringsAsFactors = FALSE)
	
      new_pheno <- rbind(new_pheno, mom_pheno)

    }

  }
}


new_pheno$numeric_part <- as.numeric(gsub("\\D", "", new_pheno$file))


# Filter rows where the numeric part is in the `refresh 2` column of corresp
filtered_data <- new_pheno[new_pheno$numeric_part %in% corresp$`refresh 2`, ]


#################### Analysis of "Qualifications" #6138  ###############################################################################


# Some double check 
file_paths <- filtered_data[filtered_data$id == "6138", "file"]  #dim = 3, there are 3 files csv 

col_1 <- filtered_data[filtered_data$id == "6138" & filtered_data$"file" == file_paths[1], "col"]
temp_pheno_all_1 <- as.data.table(fread(file_paths[1], header = T, select = c("eid", unlist(strsplit(col_1, "_")))))

col_2 <- filtered_data[filtered_data$id == "6138" & filtered_data$"file" == file_paths[2], "col"]
temp_pheno_all_2 <- as.data.table(fread(file_paths[2], header = T, select = c("eid", unlist(strsplit(col_2, "_")))))

col_3 <- filtered_data[filtered_data$id == "6138" & filtered_data$"file" == file_paths[3], "col"]
temp_pheno_all_3 <- as.data.table(fread(file_paths[3], header = T, select = c("eid", unlist(strsplit(col_3, "_")))))

# Check dimensions
if (!identical(dim(temp_pheno_all_1), dim(temp_pheno_all_2)) || !identical(dim(temp_pheno_all_2), dim(temp_pheno_all_3))) {
  print("Datasets have different dimensions")
} else {
  # Check column names and types
  if (!identical(names(temp_pheno_all_1), names(temp_pheno_all_2)) || !identical(names(temp_pheno_all_2), names(temp_pheno_all_3)) ||
      !identical(sapply(temp_pheno_all_1, class), sapply(temp_pheno_all_2, class)) || !identical(sapply(temp_pheno_all_2, class), sapply(temp_pheno_all_3, class))) {
    print("Datasets have different column names or types")
  } else {
    # Check data equality
    if (all.equal(temp_pheno_all_1, temp_pheno_all_2) && all.equal(temp_pheno_all_2, temp_pheno_all_3)) {
      print("Datasets are identical")
    } else {
      print("Datasets have different data")
    }
  }
}

# Datasets are identical - just take one ! ##########################################################

col <- filtered_data[filtered_data$id == "6138" & filtered_data$"file" == file_paths[1], "col"]
data_qualif_init <- as.data.frame(fread(file_paths[1], header = T, select = c("eid", unlist(strsplit(col, "_")))))

# Transform column related to ISCED coding - mapping to UKBB 

UKB_coding= c(1,2,3,4,5,6,-7, -3)
ISCED_coding_years= c(20,13,10,10,19,15,7,NA)
mapping <- setNames(ISCED_coding_years, UKB_coding)

data_qualif <- data_qualif_init %>%
  mutate(across(2:25, ~ mapping[as.character(.)]))

# highest qualification (max ISCED years) per instance
columns_to_min_0 <- grep("^6138-0\\.", colnames(data_qualif), value = TRUE)
data_qualif$min_0 <- apply(data_qualif[, columns_to_min_0], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))

columns_to_min_1 <- grep("^6138-1\\.", colnames(data_qualif), value = TRUE)
data_qualif$min_1 <- apply(data_qualif[, columns_to_min_1], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))

columns_to_min_2 <- grep("^6138-2\\.", colnames(data_qualif), value = TRUE)
data_qualif$min_2 <- apply(data_qualif[, columns_to_min_2], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))

columns_to_min_3 <- grep("^6138-3\\.", colnames(data_qualif), value = TRUE)
data_qualif$min_3 <- apply(data_qualif[, columns_to_min_3], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))

# Across instances take the average answer

columns_to_avg_tot <- c("min_0", "min_1", "min_2", "min_3") 

# Calculate the row-wise mean, handling NAs properly
data_qualif$avg_abs <- apply(data_qualif[, columns_to_avg_tot], 1, function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(mean(x, na.rm = TRUE))
  }
})

#################### Analysis of "Age completed full education" #845 ###############################################################################


# 845 - "Age completed full education"

file_paths <- filtered_data[filtered_data$id == "845", "file"]
col <- filtered_data[filtered_data$id == "845" & filtered_data$"file" == file_paths[1], "col"]
age_edu <- as.data.frame(fread(file_paths[1], header = T, select = c("eid", unlist(strsplit(col, "_")))))


# Address special values 

age_edu[age_edu < 0 ] <- NA

# Compute averave across istances 
age_edu$age_edu_avg <- rowMeans(age_edu[,c(2,3,4)], na.rm=T)

# Age completed education (not used downstream) - Tabea way:

dq <- data_qualif[,c("eid", "avg_abs")]
ae <- age_edu

dq_merge_ae <- merge(dq, ae, by = "eid")
dq_merge_ae[is.nan(dq_merge_ae$"age_edu_avg"),"age_edu_avg"] <- NA

dq_merge_ae$age_edu_avg_transf <- dq_merge_ae$"age_edu_avg"

dq_merge_ae$age_edu_avg_transf[which(dq_merge_ae$"age_edu_avg" < 0)] <- NA 
dq_merge_ae$age_edu_avg_transf[which(dq_merge_ae$"age_edu_avg" < 14)] <- 14 
dq_merge_ae$age_edu_avg_transf[which(dq_merge_ae$"age_edu_avg" > 19)] <- 19 
dq_merge_ae$age_edu_avg_transf[which(dq_merge_ae$"avg_abs" == 1)] <- 19

print("Save:")

dq_save <- dq_merge_ae[,c(-2)]
#save(data_qualif, dq_save, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/input_output_files_creation/input_files/covariates/data/pheno_education_att.RData")

save(data_qualif, dq_save, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/pheno_education_att.RData")


