# Load necessary libraries
library(data.table)

# Load data for sample selection - sample quality control (QC); This corresponds to the "ukb_sqc_v2.txt" file available from the UKBB portal
sample_qc <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/pheno/ukb_sqc_v2.txt",
                                 header = FALSE, select = c(3:68)))

colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", 
  "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", 
  "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", 
  "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", 
  paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")

# Sample eid; from the .fam file
sample_eid <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/plink/_001_ukb_cal_chr1_v2.fam",
                                  header = FALSE, select = c(1,5), col.names = c("eid", "sex")))

# Merge these files
df <- cbind(sample_eid, sample_qc)
print(paste0("Start: ", nrow(df), " individuals"))

# Start the Filtering 

# 1. Exclude related samples (pca_calculation = 1)
df <- df[df$pca_calculation == 1, ]
print(paste0("STEP 1: Exclude related samples: ", nrow(df), " individuals"))

# 2. Exclude non-white, non-British ancestry samples
df <- df[df$white_british == 1, ]
print(paste0("STEP 2: Exclude non-white, British ancestry samples: ", nrow(df), " individuals"))

# 3. Exclude redacted samples
df <- df[!(df$eid < 0), ]
print(paste0("STEP 3: Exclude redacted samples: ", nrow(df), " individuals"))

# 4. Exclude retracted samples; file from UKBB portal
retracted <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/org/w16389_2023-08-21.csv",
                                 header = FALSE, col.names = "eid"))
# Ensure same type for matching
df$eid <- as.numeric(df$eid)
retracted$eid <- as.numeric(retracted$eid)

df <- df[!df$eid %in% retracted$eid, ]
print(paste0("STEP 4: Exclude retracted samples: ", nrow(df), " individuals"))

out <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_filtering_sample"
fwrite(df, out)
