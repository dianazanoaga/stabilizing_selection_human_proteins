rm(list = ls())

# ---- Libraries ----
library(readxl)     
library(data.table) 

# Regression results from the cluster (regression/01 and 02)

paths <- list(
  nchild_edu   = "~/selection_proteins/data_last/protein_regression_n.children_edu_eduxsex_cov_.xlsx",
  nchil_noedu  = "~/selection_proteins/data_last/protein_regression_n.children_no_other_cov_.xlsx",
  longevity    = "~/selection_proteins/data_last/protein_regression_longevity_no_other_cov_.xlsx",
  health_score = "~/selection_proteins/data_last/protein_regression_health_score_no_other_cov_.xlsx",
  tea_intake   = "~/selection_proteins/data_last/tea_intake_all_no_edu_regression_results_proteins.xlsx"
)


.make_dist_name <- function(th) paste0("dist_", gsub("\\.", "_", as.character(th)))

preprocess_trait_data <- function(filepath,
                                  output_prefix,
                                  pval_quad_col = "PValues_quad",  # name of quadratic p-value column
                                  pval_lin_col  = "PValues_lin",   # name of linear p-value column 
                                  lin_col       = "Coefficients_lin",
                                  quad_col      = "Coefficients_quad",
                                  thresh_vals   = c(0.5, 1),       # compute dist_* for each of these
                                  save_tsv      = TRUE) {          # whether to write a TSV snapshot
  
  # Read Excel as a plain data.frame 
  df <- as.data.frame(read_excel(filepath))
  
  # Optional cleanup: PSCA gene naming fix (if present) 
  if ("Gene" %in% names(df)) {
    df$Gene[df$Gene == "PSCA.x"] <- "PSCA"
  }
  
  # Sanity-check: require quad p-values and both coefficients 
  required_cols <- c(pval_quad_col, lin_col, quad_col)
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", filepath, ": ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Multiple-testing correction (FDR) 
  df$fdr_quad <- p.adjust(df[[pval_quad_col]], method = "fdr")
  if (pval_lin_col %in% names(df)) {
    df$fdr_lin <- p.adjust(df[[pval_lin_col]], method = "fdr")
  } else {
    df$fdr_lin <- NA_real_
  }
  
  # Optimum of the quadratic curve: s = -lin/(2*quad)
  df$optimum <- ifelse(df[[quad_col]] == 0, NA_real_,
                       -df[[lin_col]] / (2 * df[[quad_col]]))
  df$s <- df$optimum
  
  for (th in thresh_vals) {
    dist_col <- .make_dist_name(th)  # "dist_0_5" or "dist_1"
    df[[dist_col]] <- ifelse(!is.na(df$s) & df$s > -th & df$s < th, 1, 0)
  }
  
  if (isTRUE(save_tsv)) {
    out_path <- paste0("~/selection_proteins/data_last/process/p_", output_prefix, ".tsv")
    fwrite(df, out_path, sep = "\t")
    message("Saved processed data to: ", out_path)
  }
  
  return(df)
}


preprocess_all <- function(paths, thresh_vals = c(0.5, 1), save_tsv = TRUE) {
  mapply(
    function(fp, nm) {
      preprocess_trait_data(
        filepath      = fp,
        output_prefix = nm,           # used only for TSV name
        thresh_vals   = thresh_vals,  # compute both dist_0_5 and dist_1
        save_tsv      = save_tsv
      )
    },
    fp = paths,
    nm = names(paths),
    SIMPLIFY = FALSE
  )
}

# Writes p_<proxy>.tsv to data_last/process/
pre_list <- preprocess_all(paths, thresh_vals = c(0.5, 1), save_tsv = TRUE)
