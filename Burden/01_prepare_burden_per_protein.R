### Load libraries

library(data.table)

base_dir <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/tcavinat/Gene_burden_all_genes_vs_all_proteins"

# all chromosome zip files
zip_files <- list.files(
  base_dir,
  pattern = "step2_gene_burden_prot.all_vs_all.chr[0-9]+\\.parallel\\.zip$",
  full.names = TRUE
)

prot_list <- list()

#zf <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/tcavinat/Gene_burden_all_genes_vs_all_proteins/step2_gene_burden_prot.all_vs_all.chr9.parallel.zip" 

for (zf in zip_files) {

  cat("Processing zip:", zf, "\n")
  
  # list files inside zip
  z_contents <- unzip(zf, list = TRUE)$Name
  # keep only *.regenie.gz
  regenies <- z_contents[grepl("\\.regenie\\.gz$", z_contents)]
  
  for (rf in regenies) {
    print(rf)
  	#rf <- "step2_gene_burden_prot.all_vs_all.batch_1.chr9_olink_A1BG.regenie.gz"   

    # protein name = bit before ".regenie.gz", after "_olink_"
    prot <- sub(".*_olink_([^.]*)\\.regenie\\.gz$", "\\1", basename(rf))
    
    # unzip this single file to tempdir
    tmpdir <- tempdir()
    gz_path <- unzip(zf, files = rf, exdir = tmpdir, overwrite = TRUE)
    
    # read compressed regenie file
    dt <- fread(cmd = paste("zcat", shQuote(gz_path)))
    
    # keep only TEST == "ADD" and ALLELE1 == "M1.0.001"
    dt_sub <- dt[TEST == "ADD" & ALLELE1 == "M1.0.001",]
    
    if (nrow(dt_sub) == 0) next
    
    # tag with protein name
    dt_sub[, PROT := prot]
    
    # rbind with existing entries for this protein
    if (is.null(prot_list[[prot]])) {
      prot_list[[prot]] <- dt_sub
    } else {
      prot_list[[prot]] <- rbindlist(list(prot_list[[prot]], dt_sub), use.names = TRUE)
    }
  }
}


# prot_list is a list: it contains one data.table per protein, across all chromosomes.
# Example: 
a1bg_dt <- prot_list[["A1BG"]]

base_dir_f <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/likelihood_zoltan_for_gene_burden/data/"

# save 
out_dir <- file.path(base_dir_f, "per_protein_0.001")
dir.create(out_dir, showWarnings = FALSE)

for (prot in names(prot_list)) {
  fwrite(prot_list[[prot]],
         file = file.path(out_dir, paste0("gene_burden_", prot, ".tsv")),
         sep = "\t")
}

