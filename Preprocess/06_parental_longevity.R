# Load necessary libraries

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(rbgen)
library(genio)
library(Matrix)
library(mbend)
library(gridExtra)
library(readxl)
library(ggplot2)

# Set the directory containing the CSV files
folder_path <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/pheno"

# Get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Load file with refresh informations about the files
#corresp <- as.data.frame(read_excel("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/eqtl_g/03_Regression/data/correspondence_uk_pheno.xlsx"))

corresp <- as.data.frame(read_excel("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/correspondence_uk_pheno.xlsx"))

# Column of interest
column_of_interest <- c("1807", "3526") # respectively : longevity father, longevity mother
													   
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

# Save the numeric part of the files and check which one correspond to refresh 2
new_pheno$numeric_part <- as.numeric(gsub("\\D", "", new_pheno$file))

# Filter rows where the numeric part is in the `refresh 2` column of corresp
filtered_data <- new_pheno[new_pheno$numeric_part %in% corresp$`refresh 2`, ]


#################### Analysis of longevity - father" # 1807  #####################################################

file_paths <- filtered_data[filtered_data$id == "1807", "file"] #take file path where I have found #1807
col <- filtered_data[filtered_data$id == "1807" & filtered_data$file == file_paths[1], "col"] #take all columns name I have found for the phenotype #1807

long_father <- as.data.frame(fread(file_paths[1], header = T, select = c("eid", unlist(strsplit(col, "_"))))) #access the file, select only eid and the columns related to #1807

# Put NA where values are less than 0 (-1 and -3, corresponding to "Do not know" and "Prefer not to answer") and take average across instances 


long_father[long_father < 0 ] <- NA
long_father$long_father_avg <- rowMeans(long_father[,c(2:ncol(long_father))], na.rm=T)


#################### Analysis of longevity mother  #3526   ##################################################


file_paths <- filtered_data[filtered_data$id == "3526", "file"]
col <- filtered_data[filtered_data$id == "3526" & filtered_data$"file" == file_paths[1], "col"]

long_mother <- as.data.frame(fread(file_paths[1], header = T, select = c("eid", unlist(strsplit(col, "_"))))) #access the file, select only eid and the columns related to #3526

# Put NA where values are less than 0 (-1 and -3, corresponding to "Do not know" and "Prefer not to answer") and take average across instances 


long_mother[long_mother < 0 ] <- NA
long_mother$long_mother_avg <- rowMeans(long_mother[,c(2:ncol(long_mother))], na.rm=T)


### inverse normal quantile transform both traits

inverse_normal_transform <- function(x) {
  # Count the non-missing values
  n <- sum(!is.na(x))
  # Rank the values, keeping NA values as NA
  ranks <- rank(x, na.last = "keep", ties.method = "average")
  # Compute the inverse normal transform using (rank - 0.5)/n
  transformed <- qnorm((ranks - 0.5) / n)
  return(transformed)
}

### merge both 

df_both <- merge(long_father, long_mother, by = "eid")

# Apply the transformation to both traits
df_both$maternal_trans <- inverse_normal_transform(df_both$long_mother_avg)
df_both$paternal_trans <- inverse_normal_transform(df_both$long_father_avg)

# now take average

df_both$parental_avg <- rowMeans(df_both[, c("maternal_trans", "paternal_trans")], na.rm = TRUE)

#save(df_both, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/input_output_files_creation/output_files/longevity /longevity.RData")
save(df_both, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/longevity.RData")

print("END!")


