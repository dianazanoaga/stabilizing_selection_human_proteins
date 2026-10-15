# ##############################################################################
# ##############################################################################
# BiocManager::install(c(
#   "Biostrings",
#   "GenomicRanges",
#   "GenomicFeatures",
#   "Rsamtools",
#   "rtracklayer",
#   "pwalign"
# ))
# 
# # install CRAN dependencies
# install.packages(c("doParallel", "foreach", "ape", "Rdpack", "benchmarkme", "devtools"))
# devtools::install_github("drostlab/orthologr")
# 
# 
# # install BLAST dependency metablastr from GitHub
# devtools::install_github("drostlab/metablastr")
# 
# # install DIAMOND dependency rdiamond from GitHub
# devtools::install_github("drostlab/rdiamond")
# ##############################################################################
# ##############################################################################

library("doParallel")
library("foreach")
library("ape") 
library("Rdpack")
library("benchmarkme")
library("devtools")
library(ape)
library(orthologr)
library(rdiamond)
library(metablastr)
library(biomartr)
library(Biostrings)
library(data.table)

# download all coding sequences for Mus musculus
Mmusculus_file <- biomartr::getCDS(organism = "Mus musculus", path = getwd())
# download all coding sequences for Homo sapiens
Hsapiens_file <- biomartr::getCDS(organism = "Homo sapiens", path = getwd())
# compute dN/dS values for Homo sapiens versus Mus musculus
Hs_vs_Mm_dNdS <- 
  dNdS(query_file      = Hsapiens_file,
       subject_file    = Mmusculus_file,
       delete_corrupt_cds = FALSE,
       ortho_detection = "RBH", 
       aa_aln_type     = "pairwise",
       aa_aln_tool     = "NW", 
       codon_aln_tool  = "pal2nal", 
       dnds_est.method = "Comeron", 
       comp_cores      = 1 )

fasta_file <- Hsapiens_file
x <- readDNAStringSet(fasta_file)

hdr <- names(x)
hdr2 <- as.character(hdr)

query_id <- sub(" .*", "", hdr2)
gene   <- sub(".*\\[gene=([^]]+)\\].*", "\\1", hdr2)

result <- data.frame(
 query_id = query_id,
  gene = gene,
  stringsAsFactors = FALSE
)

head(result)

merged <- merge(result, Hs_vs_Mm_dNdS , by = "query_id")

fwrite(merged, "~/selection_proteins/data_manuscript_2026_april/dn_ds_per_gene_human_mouse")


