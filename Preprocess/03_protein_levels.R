# load necessary libraries

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

##########  Load data

cov_pheno_filter <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_all_pheno_cov_filter_removed_and_add"))
cov_pheno_filter <- cov_pheno_filter[,-1] #V1 column added in the loading

print(paste0("Number of rows: ", nrow(cov_pheno_filter), "; Number of columns: ", ncol(cov_pheno_filter)))

# Order eid 
cov_pheno_filter <- cov_pheno_filter[order(cov_pheno_filter$eid),]

# Load file containing individual eid, protein id and protein level 

prot <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/phenotypes/olink_data.txt.gz", header = T)
prot <- as.data.frame(prot)


# Load file with protein name to match to protein id 

names_prot <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/codings/coding143.tsv"))

names_prot_to_save <- names_prot
split_meaning <- strsplit(names_prot_to_save$meaning, ";")
names_prot_to_save$prot <- sapply(split_meaning, function(x) x[1])

# Merge to have also protein name not only id protein 
prot_merged <- merge(prot, names_prot, by.x = "protein_id", by.y = "coding")

# Split the strings in the "meaning" column by semicolon - e.g. A1BG;Alpha-1B-glycoprotein
split_meaning <- strsplit(prot_merged$meaning, ";")

# Extract only the first part from each split string which is the one with the protein name - e.g. A1BG
prot_merged$prot <- sapply(split_meaning, function(x) x[1])

prot_merged_1_ins <- prot_merged[prot_merged$ins_index == 0, ]

# Take only column eid, level and name of protein 
prot_final <- prot_merged_1_ins[,c(2,4,6)]

# Order by eid again
prot_final <- prot_final[order(prot_final$eid),]

# Make a wide format having each protein as column 
final_prot_wide <- prot_final %>% 
    pivot_wider(names_from = prot, values_from = result)

final_prot_wide <- as.data.frame(final_prot_wide)

# Replace NULL with NA
final_prot_wide[final_prot_wide == "NULL"] <- NA

final_prot_wide <- as.data.frame(final_prot_wide)

# Take only those sample filtered 
sample_filter <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_filtering_sample"))

final_prot_wide_fin <- final_prot_wide[final_prot_wide$eid %in% sample_filter$eid, ]

fwrite(final_prot_wide, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/prot_file_for_regression_all_Sample")
fwrite(final_prot_wide_fin, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/prot_file_for_regression_filtered_Sample")

f <- final_prot_wide_fin[,c(1:2924)] 
f[, 2:2924] <- lapply(f[, 2:2924], as.numeric)

fwrite(f, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/prot_file_for_regression_filtered_Sample_only_prot_ids")




