library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(writexl)

cov_pheno_filter_new <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/cov_pheno_filter.csv"))

pheno_subset <- cov_pheno_filter_new %>%
  select(
    albumin, ALP, ALT, ApoA, ApoB,
    AST, bilirubin, BMI, BMR, body_fat_mass,
    Ca, cholesterol, Cr, CRP, cystatinC,
    diastolic_BP, fluid_intelligence, FVC, GGT, glucose,
    grip_strength_right, HbA1c, HDL, heart_rate,
    height, hip_circumference, IGF1, LDL, LPA,
    lymphocyte_count, MCH, monocyte_count, neuroticism, neutrophil_count,
    phosphate, platelet_count, RBC_count, reticulocyte_count, SHBG,
    systolic_BP, TG, total_protein, urate,
    urea, vitamin_D, waist_circumference, WBC_count, weight
  )
  

trait_categories <- tribble(
  ~trait,               ~category,
  "albumin",            "Liver Function",
  "ALP",                "Liver Function",
  "ALT",                "Liver Function",
  "AST",                "Liver Function",
  "bilirubin",          "Liver Function",
  "GGT",                "Liver Function",
  "total_protein",      "Liver Function",
  "ApoA",               "Lipid Metabolism",
  "ApoB",               "Lipid Metabolism",
  "cholesterol",        "Lipid Metabolism",
  "HDL",                "Lipid Metabolism",
  "LDL",                "Lipid Metabolism",
  "TG",                 "Lipid Metabolism",
  "LPA",                "Lipid Metabolism",
  "glucose",            "Glycaemic Control",
  "HbA1c",              "Glycaemic Control",
  "Cr",                 "Renal Function",
  "cystatinC",          "Renal Function",
  "urea",               "Renal Function",
  "urate",              "Renal Function",
  "systolic_BP",        "Cardiovascular Vitals",
  "diastolic_BP",       "Cardiovascular Vitals",
  "heart_rate",         "Cardiovascular Vitals",
  "BMI",                "Anthropometry",
  "BMR",                "Anthropometry",
  "body_fat_mass",      "Anthropometry",
  "height",             "Anthropometry",
  "weight",             "Anthropometry",
  "waist_circumference","Anthropometry",
  "hip_circumference",  "Anthropometry",
  "Ca",                 "Mineral Metabolism",
  "phosphate",          "Mineral Metabolism",
  "vitamin_D",          "Mineral Metabolism",
  "IGF1",               "Endocrine Function",
  "SHBG",               "Endocrine Function",
  "CRP",                "Inflammation",
  "RBC_count",          "Haematology",
  "WBC_count",          "Haematology",
  "platelet_count",     "Haematology",
  "lymphocyte_count",   "Haematology",
  "monocyte_count",     "Haematology",
  "neutrophil_count",   "Haematology",
  "reticulocyte_count", "Haematology",
  "MCH",                "Haematology",
  "FVC",                "Respiratory Function",
  "grip_strength_right","Muscle Strength",
  "fluid_intelligence", "Cognitive Function",
  "neuroticism",        "Cognitive Function"
)


desired_order <- trait_categories$trait

pheno_subset <- pheno_subset[, desired_order]

cor_matrix <- cor(pheno_subset, use = "pairwise.complete.obs", method = "spearman")

tri <- cor_matrix
tri[lower.tri(tri, diag = TRUE)] <- NA

# find all pairs with |ρ| > 0.4
high_pairs <- which(abs(tri) > 0.4, arr.ind = TRUE)

# decide which variable(s) to drop —  drop the second column in each pair
vars_to_drop <- unique(colnames(tri)[high_pairs[, "col"]])

# subset
pheno_reduced <- pheno_subset[, !(colnames(pheno_subset) %in% vars_to_drop) ]

length(vars_to_drop)              # how many removed
colnames(pheno_reduced)    # remaining variables

# z score and avg

z_mat <- scale(pheno_reduced, 
               center = TRUE, 
               scale  = TRUE)

z2_mat <- z_mat^2

mean_zscore <- rowMeans(z2_mat, na.rm = TRUE)
pheno_reduced$composite_z <- mean_zscore

df_z_to_save <- cbind(cov_pheno_filter_new[,c(1:45)], pheno_reduced$composite_z)

write_xlsx(df_z_to_save, "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/df_composite_regression_and_correlation_zscore.xlsx")



