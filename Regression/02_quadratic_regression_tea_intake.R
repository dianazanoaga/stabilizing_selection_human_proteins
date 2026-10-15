# Load necessary libraries

library(data.table)
library(dplyr)
library(stringr)
library(writexl)

# Load file containing protein levels
 
prot <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/prot_file_for_regression_filtered_Sample"))
prot <- prot[,c(1:2924)]

prot[prot == ""] <- NA

prot[, -1] <- lapply(prot[, -1], function(x) as.numeric(as.character(x)))

# Load covariates and phenotypes
cov_pheno_filter_new <- as.data.frame(fread(
  "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/cov_pheno_filter.csv"
))
# Only covariates
cov <- cov_pheno_filter_new[, 1:45]

# Load diet/physical activity
load("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/diet_physact.RData")
colnames(final_df) <- c("eid",
                        "overall_acceleration",
                        "days_week_walked",
                        "coffee_intake",
                        "tea_intake",
                        "water_intake",
                        "oily_fish_intake",
                        "cheese_intake",
                        "salt_added_to_food")


# Load education
load("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/pheno_education_att.RData")
data_qualif <- data_qualif[, c(1, 30)]
data_qualif$avg_abs <- ifelse(data_qualif$avg_abs > 19, 1,
                              ifelse(is.na(data_qualif$avg_abs), NA, 0))
data_qualif$avg_abs <- as.numeric(data_qualif$avg_abs)

# Merge covariates + education
cov_pheno_filter_1 <- merge(cov, data_qualif, by = "eid")

# Create interaction term
cov_pheno_filter_1$edu_sex <- ifelse(is.na(cov_pheno_filter_1$sex) | is.na(cov_pheno_filter_1$avg_abs),
                                     NA,
                                     cov_pheno_filter_1$sex * cov_pheno_filter_1$avg_abs)

# Merge tea intake
cov_pheno_filter <- merge(cov_pheno_filter_1, final_df[, c("eid", "tea_intake")], by = "eid")


# Take only individuals in both datasets

prot <- prot[prot$eid %in% cov_pheno_filter$eid,]
cov_pheno_filter <- cov_pheno_filter[cov_pheno_filter$eid %in% prot$eid,]


# Order individuals 
prot <- prot[order(prot$eid),]
cov_pheno_filter <- cov_pheno_filter[order(cov_pheno_filter$eid),]

if( nrow(cov_pheno_filter) == nrow(prot) ) {

	print("Same number of sample in cov_pheno and prot files.")

} else {

	print("Not same number of sample in cov_pheno and prot files.")
}

# Create 6 dataset variants 

cov_all_with_edu    <- cov_pheno_filter
cov_all_no_edu      <- cov_pheno_filter %>% select(-avg_abs, -edu_sex)
cov_male_with_edu   <- cov_pheno_filter %>% filter(sex == 1)
cov_male_no_edu     <- cov_male_with_edu %>% select(-avg_abs, -edu_sex)
cov_female_with_edu <- cov_pheno_filter %>% filter(sex == 2)
cov_female_no_edu   <- cov_female_with_edu %>% select(-avg_abs, -edu_sex)

datasets <- list(
  all_with_edu    = cov_all_with_edu,
  all_no_edu      = cov_all_no_edu,
  male_with_edu   = cov_male_with_edu,
  male_no_edu     = cov_male_no_edu,
  female_with_edu = cov_female_with_edu,
  female_no_edu   = cov_female_no_edu
)

# Regression

id <- "tea_intake"
print(paste0("Looking at phenotype: ", id))


for (set_name in names(datasets)) {
  cov_data <- datasets[[set_name]]
  print(paste0("Running model for dataset: ", set_name))

  # Align phenotype and covariate datasets
  pheno_subset_sub <- prot[prot$eid %in% cov_data$eid, ]
  pheno_subset_sub <- pheno_subset_sub[order(pheno_subset_sub$eid), ]
  cov_data <- cov_data[order(cov_data$eid), ]

  # Dynamically select the phenotype column (tea_intake)
  tea_col <- which(colnames(cov_data) == "tea_intake")
  temp_cov_pheno_filter <- cov_data %>% dplyr::select(all_of(1:tea_col))
  colnames(temp_cov_pheno_filter)[ncol(temp_cov_pheno_filter)] <- "pheno"

  # Regress out covariates
  resid <- residuals(lm(pheno ~ ., data = temp_cov_pheno_filter[, -1], na.action = na.exclude))
  resid <- scale(resid)

  # Initialize results
  cov_pheno_all <- data.frame()

# loop over proteins
for (gp in 2:ncol(pheno_subset_sub)) {
    cov_pheno_each <- data.frame(Pheno = id)
    mom <- pheno_subset_sub %>% dplyr::select(all_of(gp))
    trait_name <- colnames(pheno_subset_sub)[gp]

    print(paste0("Looking at trait: ", trait_name, "; column: ", gp))

    temp_gene_prs_filter <- data.frame(
      eid = pheno_subset_sub$eid,
      trait = scale(mom[[1]]),
      resid = resid
    )

    # Quadratic regression
    lm_prs <- lm(resid ~ trait + I(trait^2), data = temp_gene_prs_filter[, -1], na.action = na.exclude)
    lm_summary <- summary(lm_prs)

    cov_pheno_each$Gene <- trait_name
    cov_pheno_each$Coefficients_intercept <- lm_summary$coefficients[1, "Estimate"]
    cov_pheno_each$StdErrors_intercept    <- lm_summary$coefficients[1, "Std. Error"]
    cov_pheno_each$TValues_intercept      <- lm_summary$coefficients[1, "t value"]
    cov_pheno_each$PValues_intercept      <- lm_summary$coefficients[1, "Pr(>|t|)"]

    cov_pheno_each$Coefficients_lin <- lm_summary$coefficients[2, "Estimate"]
    cov_pheno_each$StdErrors_lin    <- lm_summary$coefficients[2, "Std. Error"]
    cov_pheno_each$TValues_lin      <- lm_summary$coefficients[2, "t value"]
    cov_pheno_each$PValues_lin      <- lm_summary$coefficients[2, "Pr(>|t|)"]

    cov_pheno_each$Coefficients_quad <- lm_summary$coefficients[3, "Estimate"]
    cov_pheno_each$StdErrors_quad    <- lm_summary$coefficients[3, "Std. Error"]
    cov_pheno_each$TValues_quad      <- lm_summary$coefficients[3, "t value"]
    cov_pheno_each$PValues_quad      <- lm_summary$coefficients[3, "Pr(>|t|)"]

    cov_pheno_all <- rbind(cov_pheno_all, cov_pheno_each)
  }

  # Save results
  out_file <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/01_quadratic_regression/data/",
                     id, "_", set_name, "_regression_results_proteins.xlsx")
  write_xlsx(cov_pheno_all, out_file)
  cat("Saved results for", set_name, "→", out_file, "\n")
}

cat("All 6 regressions completed successfully.\n")

