# Load libraries
library(dplyr)
library(readr)
library(readxl)
library(writexl)

# Load transcriptome data
transcriptome <- read_table("~/selection_proteins/data_last/transcriptome_liu_BM_tbl_filt80.tsv")
#transc <- transcriptome[,c("gene_name", "evolRate")]

p_a <- as.data.frame(read_excel("~/Desktop/codes/ldpred2_proteins_march/data/protein_annotation.xlsx"))
transc_with_assay_name <- merge(transcriptome, p_a[,c("ensembl_id", "Assay")], by.y = "ensembl_id", by.x = "gene_name")

bra <- transc_with_assay_name[transc_with_assay_name$tissue == "Brain",]
liv <- transc_with_assay_name[transc_with_assay_name$tissue == "Liver",]
kid<- transc_with_assay_name[transc_with_assay_name$tissue == "Kidney",]

colnames(bra)[7] <- "Pheno"
colnames(liv)[7] <- "Pheno"
colnames(kid)[7] <- "Pheno"

# ============================================================
# ANOVA: is mean transcript evolutionary rate the same across
# the 3 organs (Brain, Liver, Kidney)?
# ============================================================

anova_df <- rbind(
  data.frame(tissue = "Brain",  evolRate = bra$evolRate),
  data.frame(tissue = "Liver",  evolRate = liv$evolRate),
  data.frame(tissue = "Kidney", evolRate = kid$evolRate)
)

anova_fit <- aov(evolRate ~ tissue, data = anova_df)
summary(anova_fit)

# pairwise comparisons, in case the overall test is significant
TukeyHSD(anova_fit)

# ============================================================
# Boxplot of evolRate by organ, with group means annotated
# ============================================================

library(ggplot2)

anova_df$tissue <- factor(anova_df$tissue, levels = c("Brain", "Kidney", "Liver"))

group_means <- anova_df %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(mean_evolRate = mean(evolRate, na.rm = TRUE), .groups = "drop")

av    <- summary(anova_fit)[[1]]
f_val <- av$`F value`[1]
p_val <- av$`Pr(>F)`[1]
df1   <- av$Df[1]
df2   <- av$Df[2]

p_box <- ggplot(anova_df, aes(x = tissue, y = evolRate, fill = tissue)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  geom_text(
    data = group_means,
    aes(x = tissue, y = mean_evolRate, label = sprintf("mean = %.4f", mean_evolRate)),
    inherit.aes = FALSE,
    vjust = -1.4,
    size = 3.6,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("Brain" = "#4C72B0", "Kidney" = "#55A868", "Liver" = "#DD8452")) +
  labs(
    x = NULL,
    y = "Transcript evolutionary rate",
    title = "Transcript evolutionary rate by organ",
    subtitle = sprintf("One-way ANOVA: F(%d, %d) = %.2f, p = %.2e", df1, df2, f_val, p_val)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10.5, margin = margin(b = 8)),
    axis.text = element_text(size = 11, color = "black")
  )

p_box

### meta analysis
bra <- bra %>%
  mutate(
    se = (CI_2 - CI_1) / (2 * 1.96),
    var = se^2
  )

liv <- liv %>%
  mutate(
    se = (CI_2 - CI_1) / (2 * 1.96),
    var = se^2
  )

kid <- kid %>%
  mutate(
    se = (CI_2 - CI_1) / (2 * 1.96),
    var = se^2
  )

bra <- bra[,c(3,7,9)]
liv <- liv[,c(3,7,9)]
kid <- kid[,c(3,7,9)]

merged <- merge(
  bra, liv,
  by = "Pheno",
  all = TRUE,
  suffixes = c("_bra", "_liv")
)
colnames(kid)[colnames(kid) != "Pheno"] <-
  paste0(colnames(kid)[colnames(kid) != "Pheno"], "_kid")

merged <- merge(
  merged, kid,
  by = "Pheno",
  all = TRUE
)

merged$w_bra <- 1 / merged$var_bra
merged$w_liv <- 1 / merged$var_liv
merged$w_kid <- 1 / merged$var_kid

merged$w_bra[is.na(merged$w_bra)] <- 0
merged$w_liv[is.na(merged$w_liv)] <- 0
merged$w_kid[is.na(merged$w_kid)] <- 0

y_bra <- ifelse(is.na(merged$evolRate_bra), 0, merged$evolRate_bra)
y_liv <- ifelse(is.na(merged$evolRate_liv), 0, merged$evolRate_liv)
y_kid <- ifelse(is.na(merged$evolRate_kid), 0, merged$evolRate_kid)

w_sum <- merged$w_bra + merged$w_liv + merged$w_kid
wy_sum <- merged$w_bra * y_bra + merged$w_liv * y_liv + merged$w_kid * y_kid

merged$meta_evolRate <- wy_sum / w_sum

z <- 1.96
# meta SE and 95% CI
merged$meta_var <- 1 / w_sum
merged$meta_se  <- sqrt(merged$meta_var)

merged$meta_CI_1 <- merged$meta_evolRate - z * merged$meta_se
merged$meta_CI_2 <- merged$meta_evolRate + z * merged$meta_se

merged <- merged[,c(1,11)]

library(writexl)
write_xlsx(merged,"~/meta_analysis_tran_rate.xlsx")

### ============================================================
### META-ANALYSIS AFTER SCALING COEFFICIENTS WITHIN EACH TISSUE
### ============================================================
# Load transcriptome data
transcriptome <- read_table("~/selection_proteins/data_last/transcriptome_liu_BM_tbl_filt80.tsv")
#transc <- transcriptome[,c("gene_name", "evolRate")]

p_a <- as.data.frame(read_excel("~/Desktop/codes/ldpred2_proteins_march/data/protein_annotation.xlsx"))
transc_with_assay_name <- merge(transcriptome, p_a[,c("ensembl_id", "Assay")], by.y = "ensembl_id", by.x = "gene_name")

bra <- transc_with_assay_name[transc_with_assay_name$tissue == "Brain",]
liv <- transc_with_assay_name[transc_with_assay_name$tissue == "Liver",]
kid<- transc_with_assay_name[transc_with_assay_name$tissue == "Kidney",]

colnames(bra)[7] <- "Pheno"
colnames(liv)[7] <- "Pheno"
colnames(kid)[7] <- "Pheno"

# First calculate SE from the original 95% CI
bra <- bra %>%
  mutate(
    se_raw = (CI_2 - CI_1) / (2 * 1.96)
  )

liv <- liv %>%
  mutate(
    se_raw = (CI_2 - CI_1) / (2 * 1.96)
  )

kid <- kid %>%
  mutate(
    se_raw = (CI_2 - CI_1) / (2 * 1.96)
  )


# ============================================================
# Scale evolutionary-rate coefficients within each tissue
# scale() = mean 0, SD 1
# ============================================================

sd_bra <- sd(bra$evolRate, na.rm = TRUE)
sd_liv <- sd(liv$evolRate, na.rm = TRUE)
sd_kid <- sd(kid$evolRate, na.rm = TRUE)

mean_bra <- mean(bra$evolRate, na.rm = TRUE)
mean_liv <- mean(liv$evolRate, na.rm = TRUE)
mean_kid <- mean(kid$evolRate, na.rm = TRUE)


bra <- bra %>%
  mutate(
    evolRate_scaled = (evolRate - mean_bra) / sd_bra,
    se_scaled       = se_raw / sd_bra,
    var_scaled      = se_scaled^2
  )

liv <- liv %>%
  mutate(
    evolRate_scaled = (evolRate - mean_liv) / sd_liv,
    se_scaled       = se_raw / sd_liv,
    var_scaled      = se_scaled^2
  )

kid <- kid %>%
  mutate(
    evolRate_scaled = (evolRate - mean_kid) / sd_kid,
    se_scaled       = se_raw / sd_kid,
    var_scaled      = se_scaled^2
  )


bra_meta <- bra %>%
  dplyr::select(
    Pheno,
    evolRate = evolRate_scaled,
    var = var_scaled
  )

liv_meta <- liv %>%
  dplyr::select(
    Pheno,
    evolRate = evolRate_scaled,
    var = var_scaled
  )

kid_meta <- kid %>%
  dplyr::select(
    Pheno,
    evolRate = evolRate_scaled,
    var = var_scaled
  )


# ============================================================
# Merge tissues
# ============================================================

merged <- merge(
  bra_meta,
  liv_meta,
  by = "Pheno",
  all = TRUE,
  suffixes = c("_bra", "_liv")
)

colnames(kid_meta)[colnames(kid_meta) != "Pheno"] <-
  paste0(colnames(kid_meta)[colnames(kid_meta) != "Pheno"], "_kid")

merged <- merge(
  merged,
  kid_meta,
  by = "Pheno",
  all = TRUE
)


# ============================================================
# Inverse-variance weights
# ============================================================

merged$w_bra <- 1 / merged$var_bra
merged$w_liv <- 1 / merged$var_liv
merged$w_kid <- 1 / merged$var_kid

merged$w_bra[is.na(merged$w_bra)] <- 0
merged$w_liv[is.na(merged$w_liv)] <- 0
merged$w_kid[is.na(merged$w_kid)] <- 0


# Missing effect estimates contribute zero because their weight = 0
y_bra <- ifelse(is.na(merged$evolRate_bra), 0, merged$evolRate_bra)
y_liv <- ifelse(is.na(merged$evolRate_liv), 0, merged$evolRate_liv)
y_kid <- ifelse(is.na(merged$evolRate_kid), 0, merged$evolRate_kid)


# ============================================================
# Fixed-effect meta-analysis
# ============================================================

w_sum <- merged$w_bra +
  merged$w_liv +
  merged$w_kid

wy_sum <- merged$w_bra * y_bra +
  merged$w_liv * y_liv +
  merged$w_kid * y_kid

merged$meta_evolRate_scaled <- wy_sum / w_sum


# Meta-analysis variance and SE
merged$meta_var <- 1 / w_sum
merged$meta_se  <- sqrt(merged$meta_var)


# 95% CI
z <- 1.96

merged$meta_CI_1 <- merged$meta_evolRate_scaled -
  z * merged$meta_se

merged$meta_CI_2 <- merged$meta_evolRate_scaled +
  z * merged$meta_se


# ============================================================
# Final output
# ============================================================

meta_final <- merged %>%
  dplyr::select(
    Pheno,
    meta_evolRate_scaled,
    meta_se,
    meta_CI_1,
    meta_CI_2
  )

write_xlsx(
  meta_final,
  "~/selection_proteins/scripts_figure_may_26/data/meta_analysis_tran_rate_scaled.xlsx"
)
