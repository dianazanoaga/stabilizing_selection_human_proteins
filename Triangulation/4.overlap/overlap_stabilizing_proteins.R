# Load libraries
library(data.table)
library(dplyr)
library(purrr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(UpSetR)

ldpred2_prot_no_cis <- read_excel("~/selection_proteins/data_last/process/ldpred2-auto_proteins_coefficients_no_cis.xlsx")
ldpred2_prot_no_cis <- ldpred2_prot_no_cis[!is.na(ldpred2_prot_no_cis$alpha_est_1), ]
colnames(ldpred2_prot_no_cis)[1] <- "Pheno"

pheno_to_use <- unique(ldpred2_prot_no_cis$Pheno)

se <- (ldpred2_prot_no_cis$p_est_3 - ldpred2_prot_no_cis$p_est_2) / (2 * 1.96)
z <- (ldpred2_prot_no_cis$p_est_1 - 0) / se
ldpred2_prot_no_cis$pval_p <- 2 * (1 - pnorm(abs(z)))

se <- (ldpred2_prot_no_cis$h2_est_3 - ldpred2_prot_no_cis$h2_est_2) / (2 * 1.96)
z <- (ldpred2_prot_no_cis$h2_est_1 - 0) / se
ldpred2_prot_no_cis$pval_h2 <- 2 * (1 - pnorm(abs(z)))

ldpred2_prot_no_cis$pval_h2_fdr <- p.adjust(ldpred2_prot_no_cis$pval_h2, method = "fdr")
ldpred2_prot_no_cis$pval_p_fdr  <- p.adjust(ldpred2_prot_no_cis$pval_p, method = "fdr")

colnames(ldpred2_prot_no_cis)[8] <- "alpha_lpred2"
colnames(ldpred2_prot_no_cis)[13] <- "pval_lpred2"
colnames(ldpred2_prot_no_cis)[14] <- "pval_lpred2_adj"

ldpred2_prot_no_cis$dist <- rep(1, nrow(ldpred2_prot_no_cis))

ldpred2_prot_no_cis <- ldpred2_prot_no_cis[!duplicated(ldpred2_prot_no_cis$Pheno),]
ldpred2_prot_no_cis <- ldpred2_prot_no_cis[!duplicated(ldpred2_prot_no_cis$Pheno),]


ldpred2_prot_no_cis_to_use_low <- ldpred2_prot_no_cis[
  ldpred2_prot_no_cis$h2_est_1 > 0.1 &
    ldpred2_prot_no_cis$p_est_1 > 0.00025 &
    ldpred2_prot_no_cis$pval_h2 < 0.05 &
    ldpred2_prot_no_cis$pval_p < 0.05 &
    ldpred2_prot_no_cis$pval_lpred2 < 0.05, ]

ldpred2_prot_no_cis_to_use_low <- ldpred2_prot_no_cis_to_use_low[, c(1, 2, 8)]

ldpred2_prot_no_cis_to_use_high <- ldpred2_prot_no_cis[
  ldpred2_prot_no_cis$h2_est_1 > 0.1 &
    ldpred2_prot_no_cis$p_est_1 > 0.00025 &
    ldpred2_prot_no_cis$pval_h2 < 0.05 &
    ldpred2_prot_no_cis$pval_p < 0.05 &
    ldpred2_prot_no_cis$pval_lpred2_adj < 0.05, ]

ldpred2_prot_no_cis_to_use_high <- ldpred2_prot_no_cis_to_use_high[, c(1, 2, 8)]

### gene burden 
res_all <- fread("~/alpha_sigmaG_results_all_proteins_0.001_correction_optim.tsv")

colnames(res_all)[3] <- "alpha_burden"
colnames(res_all)[6] <- "pval_burden"
colnames(res_all)[11] <- "pval_burden_adj"

res_all$dist <- rep(1, nrow(res_all))

res_all_to_use_low <- res_all[res_all$pval_burden < 0.05, ]
res_all_to_use_low <- res_all_to_use_low[, c(2, 3)]

res_all_to_use_high <- res_all[res_all$pval_burden_adj < 0.05, ]
res_all_to_use_high <- res_all_to_use_high[, c(2, 3)]


### prot_lev

prot_lev_nchild <- fread("~/selection_proteins/data_last/process/p_nchild_edu.tsv")
prot_lev_long <- fread("~/selection_proteins/data_last/process/p_longevity.tsv")
prot_lev_compscore <- fread("~/selection_proteins/data_last/process/p_health_score.tsv")

colnames(prot_lev_nchild)[11] <- "alpha_nchild"
colnames(prot_lev_long)[11] <- "alpha_long"
colnames(prot_lev_compscore)[11] <- "alpha_health_score"

colnames(prot_lev_nchild)[14] <- "pval_nchild"
colnames(prot_lev_long)[14] <- "pval_long"
colnames(prot_lev_compscore)[14] <- "pval_health_score"

colnames(prot_lev_nchild)[15] <- "pval_nchild_adj"
colnames(prot_lev_long)[15] <- "pval_long_adj"
colnames(prot_lev_compscore)[15] <- "pval_health_score_adj"

colnames(prot_lev_nchild)[20] <- "dist"
colnames(prot_lev_long)[20] <- "dist"
colnames(prot_lev_compscore)[20] <- "dist"

prot_lev_nchild <- prot_lev_nchild[, -1]
prot_lev_long <- prot_lev_long[, -1]
prot_lev_compscore <- prot_lev_compscore[, -1]

colnames(prot_lev_nchild)[1] <- "Pheno"
colnames(prot_lev_long)[1] <- "Pheno"
colnames(prot_lev_compscore)[1] <- "Pheno"

prot_lev_nchild_to_use <- merge(prot_lev_nchild, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_nchild_to_use_low <- prot_lev_nchild_to_use[
  prot_lev_nchild_to_use$dist.x == 1 &
    prot_lev_nchild_to_use$h2_est_1 > 0.1 &
    prot_lev_nchild_to_use$pval_nchild < 0.05]
prot_lev_nchild_to_use_low <- prot_lev_nchild_to_use_low[, c(1, 10)]

prot_lev_nchild_to_use <- merge(prot_lev_nchild, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_nchild_to_use_high <- prot_lev_nchild_to_use[
  prot_lev_nchild_to_use$dist.x == 1 &
    prot_lev_nchild_to_use$h2_est_1 > 0.1 &
    prot_lev_nchild_to_use$pval_nchild_adj < 0.05]
prot_lev_nchild_to_use_high <- prot_lev_nchild_to_use_high[, c(1, 10)]

prot_lev_long_to_use <- merge(prot_lev_long, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_long_to_use_low <- prot_lev_long_to_use[
  prot_lev_long_to_use$dist.x == 1 &
    prot_lev_long_to_use$h2_est_1 > 0.1 &
    prot_lev_long_to_use$pval_long < 0.05]
prot_lev_long_to_use_low <- prot_lev_long_to_use_low[, c(1, 10)]

prot_lev_long_to_use <- merge(prot_lev_long, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_long_to_use_high <- prot_lev_long_to_use[
  prot_lev_long_to_use$dist.x == 1 &
    prot_lev_long_to_use$h2_est_1 > 0.1 &
    prot_lev_long_to_use$pval_long_adj < 0.05]
prot_lev_long_to_use_high <- prot_lev_long_to_use_high[, c(1, 10)]

prot_lev_compscore_to_use <- merge(prot_lev_compscore, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_compscore_to_use_low <- prot_lev_compscore_to_use[
  prot_lev_compscore_to_use$dist.x == 1 &
    prot_lev_compscore_to_use$h2_est_1 > 0.1 &
    prot_lev_compscore_to_use$pval_health_score < 0.05]
prot_lev_compscore_to_use_low <- prot_lev_compscore_to_use_low[, c(1, 10)]

prot_lev_compscore_to_use <- merge(prot_lev_compscore, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_compscore_to_use_high <- prot_lev_compscore_to_use[
  prot_lev_compscore_to_use$dist.x == 1 &
    prot_lev_compscore_to_use$h2_est_1 > 0.1 &
    prot_lev_compscore_to_use$pval_health_score_adj < 0.05]
prot_lev_compscore_to_use_high <- prot_lev_compscore_to_use_high[, c(1, 10)]

ldpred2_prot_no_cis_to_use_low_neg <- ldpred2_prot_no_cis_to_use_low[ldpred2_prot_no_cis_to_use_low$alpha_lpred2 <0,]
res_all_to_use_low_neg <- res_all_to_use_low[res_all_to_use_low$alpha_burden < 0,]
prot_lev_nchild_to_use_low_neg <- prot_lev_nchild_to_use_low[prot_lev_nchild_to_use_low$alpha_nchild <0,]
prot_lev_long_to_use_low_neg <- prot_lev_long_to_use_low[prot_lev_long_to_use_low$alpha_long <0,]
prot_lev_compscore_to_use_low_neg <- prot_lev_compscore_to_use_low[prot_lev_compscore_to_use_low$alpha_health_score <0,]

ldpred_sig <- unique(ldpred2_prot_no_cis_to_use_low_neg$Pheno)

burden_sig <- unique(res_all_to_use_low_neg$Pheno)

nchild_sig <- unique(prot_lev_nchild_to_use_low_neg$Pheno)

long_sig <- unique(prot_lev_long_to_use_low_neg$Pheno)

health_sig <- unique(prot_lev_compscore_to_use_low_neg$Pheno)


all_proteins <- unique(c(
  ldpred_sig,
  burden_sig,
  nchild_sig,
  long_sig,
  health_sig
))


final_table <- data.frame(
  Pheno = all_proteins
)

final_table$ldpred2 <- final_table$Pheno %in% ldpred_sig
final_table$burden <- final_table$Pheno %in% burden_sig
final_table$nchild <- final_table$Pheno %in% nchild_sig
final_table$longevity <- final_table$Pheno %in% long_sig
final_table$health_score <- final_table$Pheno %in% health_sig

final_table$n_methods <- rowSums(
  final_table[, c(
    "ldpred2",
    "burden",
    "nchild",
    "longevity",
    "health_score"
  )]
)

final_table$status <- ifelse(final_table$n_methods >= 3, "OK", "NO")

final_table <- final_table %>%
  arrange(desc(n_methods))

head(final_table)

methods <- c("ldpred2", "burden", "nchild", "longevity", "health_score")

plot_df <- final_table %>%
  filter(n_methods >= 1) %>%
  arrange(n_methods, Pheno) %>%
  mutate(
    protein_rank = row_number(),
    label = ifelse(n_methods >= 4, Pheno, NA_character_)
  )

count_df <- plot_df %>%
  count(n_methods) %>%
  mutate(
    x_pos = -30,
    label = paste0("n = ", n)
  )

p <- ggplot(plot_df, aes(x = protein_rank, y = n_methods)) +
  geom_segment(
    aes(xend = protein_rank, y = 0, yend = n_methods),
    linewidth = 0.25,
    alpha = 0.35
  ) +
  geom_point(
    aes(size = n_methods, fill = factor(n_methods)),
    shape = 21,
    color = "black",
    stroke = 0.25
  ) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  geom_text(
    data = count_df,
    aes(
      x = x_pos,
      y = n_methods,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0.7,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_y_continuous(
    breaks = 1:5,
    limits = c(0, 5.4)
  ) +
  scale_size_continuous(
    range = c(1.5, 5),
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(-60, max(plot_df$protein_rank) + 15),
    clip = "off"
  )+
  labs(
    x = "Proteins ranked by cross-method support",
    y = "Number of supporting methods",
    fill = "# methods",
    title = "Cross-method support for stabilizing selection in proteins",
    subtitle = "Proteins passing method-specific filters and nominal significance threshold (P < 0.05)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    plot.margin = margin(8, 20, 8, 80)
  )

p


methods <- c("ldpred2", "burden", "nchild", "longevity", "health_score")

upset_df <- final_table %>%
  filter(n_methods >= 1) %>%
  mutate(across(all_of(methods), ~ as.integer(as.logical(.x)))) %>%
  dplyr::select(all_of(methods))

methods <- c("ldpred2", "burden", "nchild", "longevity", "health_score")

upset_df <- final_table %>%
  filter(n_methods >= 1) %>%
  mutate(across(all_of(methods), ~ as.integer(as.logical(.x)))) %>%
  dplyr::select(all_of(methods))

# rename methods for figure
colnames(upset_df) <- c(
  "LDpred2",
  "Burden",
  "N.children",
  "Longevity",
  "Health Score"
)

method_cols <- c(
  "LDpred2" = "#1b9e77",
  "Burden" = "#d95f02",
  "N.children" = "#7570b3",
  "Longevity" = "#e7298a",
  "Health Score" = "#66a61e"
)

pdf("stabilizing_selection_upset_manuscript.pdf", width = 8.5, height = 5.5)

upset(
  as.data.frame(upset_df),
  
  sets = c(
    "LDpred2",
    "Burden",
    "Health Score",
    "Longevity",
    "N.children"
  ),
  
  nsets = ncol(upset_df),
  order.by = "freq",
  keep.order = TRUE,
  
  main.bar.color = "grey25",
  
  sets.bar.color = c(
    "#1b9e77",
    "#d95f02",
    "#66a61e",
    "#e7298a",
    "#7570b3"
  ),
  
  matrix.color = "grey20",
  shade.color = "grey90",
  shade.alpha = 0.35,
  
  point.size = 3.2,
  line.size = 0.8,
  
  mb.ratio = c(0.65, 0.35),
  
  mainbar.y.label = "Protein intersections",
  sets.x.label = "Proteins per method",
  
  text.scale = c(1.4, 1.4, 1.2, 1.2, 1.4, 1.2)
)
