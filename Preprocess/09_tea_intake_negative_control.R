# Load necessary libraries
library(data.table)
library(dplyr)
library(readxl)
library(stringr)

# Load file 

folder_path <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/pheno"
corresp_path <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/correspondence_uk_pheno.xlsx"

# Phenotype field IDs of interest
column_of_interest <- c("90012",  # accelerometer average acceleration
                        "864",    # days/week walked ≥10 min
                        "1498",   # coffee intake
                        "1488",   # tea intake
                        "1528",   # water intake
                        "1329",   # oily fish intake
                        "1408",   # cheese intake
                        "1478")   # salt added to food
# Output file

output_path <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/diet_physact.RData"


# List all CSVs in phenotype directory
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Load correspondence table
corresp <- as.data.frame(read_excel(corresp_path))

# Find phenotype columns 

new_pheno <- data.frame(id = character(),
                        file = character(),
                        col = character(),
                        stringsAsFactors = FALSE)

for (i in seq_along(column_of_interest)) {
  pheno_id <- column_of_interest[i]
  
  for (file in csv_files) {
    data_header <- fread(file, nrows = 0)
    colnames_file <- colnames(data_header)
    
    matching_indices <- grep(paste0("^", pheno_id, "-"), colnames_file)
    
    if (length(matching_indices) > 0) {
      mom_pheno <- data.frame(id = pheno_id,
                              file = file,
                              col = paste(colnames_file[matching_indices], collapse = "_"),
                              stringsAsFactors = FALSE)
      new_pheno <- rbind(new_pheno, mom_pheno)
    }
  }
}

# Keep only files corresponding to "refresh 2"
new_pheno$numeric_part <- as.numeric(gsub("\\D", "", new_pheno$file))
filtered_data <- new_pheno[new_pheno$numeric_part %in% corresp$`refresh 2`, ]


# Inverse normal quantile transformation
inverse_normal_transform <- function(x) {
  n <- sum(!is.na(x))
  ranks <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((ranks - 0.5) / n)
}


pheno_list <- list()

for (pheno_id in column_of_interest) {
  cat("Processing phenotype:", pheno_id, "\n")
  
  pheno_info <- filtered_data[filtered_data$id == pheno_id, ]
  if (nrow(pheno_info) == 0) {
    cat("  --> WARNING: Phenotype not found in filtered data.\n")
    next
  }
  
  file_path <- pheno_info$file[1]
  cols <- unlist(strsplit(pheno_info$col[1], "_"))
  
  tmp <- as.data.frame(fread(file_path, select = c("eid", cols)))
  tmp[tmp < 0] <- NA  # set invalid codes (-1, -3) to NA
  tmp$avg <- rowMeans(tmp[, -1, drop = FALSE], na.rm = TRUE)

  #tmp$avg_norm <- inverse_normal_transform(tmp$avg)

  pheno_list[[pheno_id]] <- tmp[, c("eid", "avg")]
  colnames(pheno_list[[pheno_id]]) <- c("eid", pheno_id)
}


final_df <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE), pheno_list)

# Save

save(final_df, file = output_path)

cat("\n All done! File saved to:\n", output_path, "\n")

