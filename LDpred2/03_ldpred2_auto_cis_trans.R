# load necessary library 

library(writexl)
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(readxl)
library(data.table)
library(magrittr)
library(ggplot2)
library(readxl)
library(dplyr)
library(stringr)
library(arrow)
library(writexl)
library(ieugwasr)



# file_with_all_prot <- list.files("/data/FAC/FBM/DBC/zkutalik/default_sensitive/_avdg/2025_02_11_ukb_ppp_non_normalized/harmonized/")
# "ZPR1:O75312:OID30801:v1:Neurology_II.parq" 
#write.table(file_with_all_prot, "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/linear_regression_likelihood/data/file_with_all_prot.txt", row.names = F)


args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

# Remove unwanted quotes 
filename <- gsub('^"|"$', '', filename)  # Removes quotes at the start and end

print(paste("Filename after cleaning:", filename))

#filename <- "ZPR1:O75312:OID30801:v1:Neurology_II.parq" 
# Construct the full file path
file <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/_avdg/2025_02_11_ukb_ppp_non_normalized/harmonized/", filename)

print(paste("Full file path:", file))

# Extract the first part before the first colon
first_part <- strsplit(filename, ":")[[1]][1]
pheno <- first_part
print(paste("Phenotype extracted:", pheno))

trait <- sub(".*:(.*)\\.parq", "\\1", filename)
print(paste("Trait extracted:", trait))


# Check if the file exists
if (!file.exists(file)) {
    cat("File does not exist:", file, "\nSkipping this file...\n")
} else {
    # Attempt to read the file
    sumstats_p <- tryCatch(
        {
            read_parquet(file)  # Read the Parquet file
        },
        error = function(e) {
            cat("Error reading file:", file, "\n", e$message, "\nSkipping this file...\n")
            return(NULL)  # Return NULL and skip
        }
    )
    
    if (!is.null(sumstats_p)) {
        sumstats_p <- as.data.frame(sumstats_p)
        print("Successfully read the file!")
    }
}


###################### LOAD PRE-COMPUTED LD MATRICES ############


# Directory containing precomputed LD matrices
ld_dir <-  "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/data_ld_reference/ld_blocks/"

# Temporary file for the combined sparse matrix
tmp <- tempfile(tmpdir = "tmp-data2", fileext = ".sbk")

# Initialize variables
ld <- NULL
corr <- NULL


# Loop through chromosomes
for (chr in 1:22) {
  print(paste("Processing chromosome", chr))
  
  # Load the precomputed LD matrix for this chromosome
  ld_file <- paste0(ld_dir, "LD_with_blocks_chr", chr, ".rds")
  corr0 <- readRDS(ld_file)  # Load precomputed LD matrix
  
 
  # Combine into a single SFBM object
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)  # LD scores for first chromosome
    corr <- as_SFBM(corr0, compact = TRUE)  # Initialize the SFBM
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))  # Append LD scores
    corr$add_columns(corr0, nrow(corr))  # Append columns to the SFBM
  }
}


##################################################################


##################################################################

info <- readRDS("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/data_ld_reference/map_hm3_plus.rds")



  # Skip to the next iteration if sumstats_p is NULL
  if (is.null(sumstats_p)) next
  
   	var <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/sumstat/GWAS/neale/variants.tsv.gz")

	var <- as.data.frame(var)

	sumstats_s <- merge(sumstats_p, var[,c(1:6,10, 17)], by.y = c("chr","pos"), by.x = c("chromosome", "base_pair_location"))

	sumstats <- sumstats_s[,c(10, 1,2, 3, 4, 5, 6, 11, 8, 7, 16)]

# LDpred 2 require the header to follow the exact naming
	names(sumstats) <-
    c("rsid",
    "chr",
    "pos",
    "a1",
    "a0",
    "beta",
    "beta_se",
    "N",
    "p", 
    "af1",
    "info")
    

	sumstats[,2] <- as.numeric(sumstats[,2])
	
	sumstats$n_eff <- sumstats$N


# QC
# SD REFERENCE

matched_variants <- intersect(sumstats$rsid, info$rsid)

# Get valid indices
ind.df_beta <- match(matched_variants, sumstats$rsid)
ind.ref <- match(matched_variants, info$rsid)

# Subset correctly
sums <- sumstats[ind.df_beta,]
ld_s <- ld[ind.ref]
sums_n <- merge(sums, info, by = "rsid")


# Process the dataset - check for flipped alleles between info and my sumstat 

dt_cleaned <- sums_n %>%
  mutate(
    flip = (a1.x == a0.y & a0.x == a1.y),  # Identify flipped alleles
    valid = (a1.x == a1.y & a0.x == a0.y) | flip  # Keep only valid rows
  ) %>%
  filter(valid) %>%  # Remove rows where alleles don’t match
  mutate(
    beta = ifelse(flip, -beta, beta),  # Flip beta if alleles are flipped
    temp = ifelse(flip, a1.x, a1.x),  # Temp variable to swap alleles
    a1.x = ifelse(flip, a0.x, a1.x),  
    a0.x = ifelse(flip, temp, a0.x),  
    af1 = ifelse(flip, 1 - af1, af1)  # Flip allele frequency if needed
  ) %>%
  select(-flip, -valid, -temp)  # Remove helper columns


dt_cleaned <- dt_cleaned %>%
  mutate(maf = pmin(af1, 1 - af1)) %>%
  filter(maf > 0.01) 



protein_annotation <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/sumstat/QTL/protein/ukb_ppp_cohort/protein_annotation_3k.tsv_orig"))
protein_annotation_short <- protein_annotation[,c("UniProt", "Assay", "chr", "gene_start", "gene_end")]
prot_ann_one <- protein_annotation_short[protein_annotation_short$Assay == pheno,]
gene_start_fin <- prot_ann_one$gene_start - 1000000
gene_end_fin <- prot_ann_one$gene_end + 1000000
chr_gene <- prot_ann_one$chr 

dt_cis <- dt_cleaned[dt_cleaned$chr.x == chr_gene & dt_cleaned$pos.x >= gene_start_fin & dt_cleaned$pos.x <= gene_end_fin,]

dt_no_cis <- dt_cleaned[!(dt_cleaned$chr.x == chr_gene & dt_cleaned$pos.x >= gene_start_fin & dt_cleaned$pos.x <= gene_end_fin),]

matched_variants <- intersect(dt_cleaned$rsid, info$rsid)

# Get valid indices
ind.df_beta <- match(matched_variants, dt_cleaned$rsid)
ind.ref <- match(matched_variants, info$rsid)

# Subset correctly
sums_new <- dt_cleaned[ind.df_beta,]
ld_s_new <- ld[ind.ref]

sums_new_new <- sums_new

# Define output file path
output_file <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/proteins/data/info_snps/results_", pheno,"_", trait, ".gz")

# Write and gzip the file
gz_out <- gzfile(output_file, "wt")  # "wt" = write text

write.table(sums_new_new, gz_out, row.names = FALSE, quote = FALSE, sep = "\t")
close(gz_out)  # Close the connection to finalize compression


# Estimate of h2 from LD Score regression
(ldsc <- with(sums_new, snp_ldsc(ld_s_new, 1444196 , chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
                                

ldsc_h2_est <- ldsc[["h2"]]

NCORES = 4

set.seed(1)

coef_shrink <- 0.95 

start_time <- Sys.time()  # Start the timer

if (ldsc_h2_est < 0) {
  ldsc_h2_est <- 0.05
}


multi_auto <- snp_ldpred2_auto(
  corr, sums_new, h2_init = ldsc_h2_est, ind.corr = ind.ref,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 50), ncores = NCORES,
  burn_in = 500, num_iter = 500, report_step = 20, alpha_bounds= c(-5,1),
  allow_jump_sign = FALSE, shrink_corr = coef_shrink
)

end_time <- Sys.time()  # End the timer

# Calculate the duration
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))


(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
alp <- quantile(all_alpha, c(0.5, 0.025, 0.975), na.rm = T)

all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
h2 <- quantile(all_h2, c(0.5, 0.025, 0.975))

all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
p <- quantile(all_p, c(0.5, 0.025, 0.975))

df_t <- cbind(pheno, h2[1],h2[2],h2[3], p[1], p[2], p[3], alp[1], alp[2], alp[3])

colnames(df_t) <- c("pheno", "h2_est_1", "h2_est_2", "h2_est_3", "p_est_1", "p_est_2", "p_est_3", "alpha_est_1", "alpha_est_2", "alpha_est_3")

df_t <- as.data.frame(df_t)

write_xlsx(df_t, paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/proteins/data/unstand_prot_coef/results_", pheno,"_", trait, ".xlsx"))

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
beta_est_df <- cbind(matched_variants, beta_auto)

beta_est_df <- as.data.frame(beta_est_df)
beta_est_df[,2] <- as.numeric(beta_est_df[,2])

write.table(beta_est_df, paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/02_ldpred2/proteins/data/unstand_prot_beta/results_", pheno,"_", trait, ".xlsx"))


