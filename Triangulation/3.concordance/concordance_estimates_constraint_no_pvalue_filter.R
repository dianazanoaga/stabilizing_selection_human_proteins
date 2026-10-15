# Load libraries 

library(data.table)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())
gc()

# ============================================================
# Helper functions
# ============================================================

dedup_by_pheno_mean <- function(dt, name) {
  dt <- as.data.table(dt)
  if (!"Pheno" %in% names(dt)) stop(name, " is missing 'Pheno'")
  
  numcols <- setdiff(names(dt)[sapply(dt, is.numeric)], "Pheno")
  if (length(numcols) == 0) stop(name, " has no numeric columns")
  
  dt[, lapply(.SD, function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }), by = Pheno, .SDcols = numcols]
}

prep_var <- function(dt, value_col, var_name, name_for_checks) {
  dt <- as.data.table(dt)[, .(Pheno, value = get(value_col))]
  setnames(dt, "value", var_name)
  dedup_by_pheno_mean(dt, name_for_checks)
}

spearman_xy <- function(x, y, expected_sign = NULL) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  
  if (n < 3) {
    return(list(rho = NA_real_, p = NA_real_, n = n))
  }
  
  ct <- suppressWarnings(
    cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  )
  
  rho <- unname(ct$estimate)
  p_two <- ct$p.value
  
  if (is.null(expected_sign)) {
    p <- p_two
  } else {
    if (sign(rho) == expected_sign) {
      p <- p_two / 2
    } else {
      p <- 1 - p_two / 2
    }
  }
  
  list(rho = rho, p = p, n = n)
}

pairwise_spearman_dt <- function(dt1, v1, dt2, v2, expected_sign = NULL) {
  dt1 <- as.data.table(dt1)
  dt2 <- as.data.table(dt2)
  
  tmp1 <- dt1[, c("Pheno", v1), with = FALSE]
  setnames(tmp1, v1, "x")
  
  tmp2 <- dt2[, c("Pheno", v2), with = FALSE]
  setnames(tmp2, v2, "y")
  
  tmp <- merge(tmp1, tmp2, by = "Pheno", all = FALSE)
  
  spearman_xy(tmp$x, tmp$y, expected_sign)
}

# ============================================================
# Load LDpred2
# ============================================================

ldpred2_prot_no_cis <- read_excel(
  "~/selection_proteins/data_last/process/ldpred2-auto_proteins_coefficients_no_cis.xlsx"
)

ldpred2_prot_no_cis <- ldpred2_prot_no_cis[
  !is.na(ldpred2_prot_no_cis$alpha_est_1), ]

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

ldpred2_prot_no_cis$dist <- 1
ldpred2_prot_no_cis <- ldpred2_prot_no_cis[
  !duplicated(ldpred2_prot_no_cis$Pheno), ]

# p < 1
ldpred2_to_use <- ldpred2_prot_no_cis[
  ldpred2_prot_no_cis$h2_est_1 > 0.1 &
    ldpred2_prot_no_cis$p_est_1 > 0.00025 &
    ldpred2_prot_no_cis$pval_h2 < 0.05 &
    ldpred2_prot_no_cis$pval_p < 0.05 &
    ldpred2_prot_no_cis$pval_lpred2 < 1,
]

ldpred2_to_use <- ldpred2_to_use[, c(1, 2, 8)]

# ============================================================
# Load burden
# ============================================================

res_all <- fread("~/alpha_sigmaG_results_all_proteins_0.001_correction_optim.tsv")

colnames(res_all)[3] <- "alpha_burden"
colnames(res_all)[6] <- "pval_burden"
colnames(res_all)[11] <- "pval_burden_adj"

res_all$dist <- 1

res_all_to_use <- res_all[res_all$pval_burden < 1, ]
res_all_to_use <- res_all_to_use[, c(2, 3)]

# ============================================================
# Load observed fitness proxy regressions
# ============================================================

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

# p < 1 
prot_lev_nchild_to_use <- merge(prot_lev_nchild, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_nchild_to_use <- prot_lev_nchild_to_use[
  prot_lev_nchild_to_use$dist.x == 1 &
    prot_lev_nchild_to_use$h2_est_1 > 0.1 &
    prot_lev_nchild_to_use$pval_nchild < 1, ]
prot_lev_nchild_to_use <- prot_lev_nchild_to_use[, c(1, 10)]

prot_lev_long_to_use <- merge(prot_lev_long, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_long_to_use <- prot_lev_long_to_use[
  prot_lev_long_to_use$dist.x == 1 &
    prot_lev_long_to_use$h2_est_1 > 0.1 &
    prot_lev_long_to_use$pval_long < 1, ]
prot_lev_long_to_use <- prot_lev_long_to_use[, c(1, 10)]

prot_lev_compscore_to_use <- merge(prot_lev_compscore, ldpred2_prot_no_cis, by = "Pheno")
prot_lev_compscore_to_use <- prot_lev_compscore_to_use[
  prot_lev_compscore_to_use$dist.x == 1 &
    prot_lev_compscore_to_use$h2_est_1 > 0.1 &
    prot_lev_compscore_to_use$pval_health_score < 1, ]
prot_lev_compscore_to_use <- prot_lev_compscore_to_use[, c(1, 10)]



# ============================================================
# Load constraint metrics
# ============================================================

gnomad_to_use <- fread(
  "~/selection_proteins/data_last/gnomad.v4.1.constraint_metrics.tsv",
  select = c("gene", "mane_select", "canonical", "lof_hc_lc.pLI", "lof.oe_ci.upper")
)

gnomad_to_use <- gnomad_to_use[
  mane_select == TRUE | canonical == TRUE,
  .(Pheno = gene, pLI = `lof_hc_lc.pLI`, LOEUF = `lof.oe_ci.upper`)
]

gnomad_to_use <- unique(gnomad_to_use, by = "Pheno")
gnomad_to_use<- gnomad_to_use[gnomad_to_use$Pheno %in% pheno_to_use,]


ev_tr <- read_xlsx("~/selection_proteins/scripts_figure_may_26/data/meta_analysis_tran_rate_scaled.xlsx")

colnames(ev_tr) <- c("Pheno", "meta_transcript_ev_rate")

coef_reg <- read_excel(
  "~/selection_proteins/data/lasso/coefficients_proteins_traits_lm_scaled_lasso_lambdamin.xlsx"
)

coef_reg[coef_reg$prot == "PSCA.x", ]$prot <- "PSCA"

coef_reg_u <- coef_reg[!duplicated(coef_reg$prot), ]
colnames(coef_reg_u)[3] <- "Pheno"
coef_reg_u <- coef_reg_u[, c(3, 5)]
colnames(coef_reg_u) <- c("Pheno", "Lasso R2")

coef_reg_u <- coef_reg_u[coef_reg_u$Pheno %in% pheno_to_use,]


dnds <- as.data.frame(
  fread("~/selection_proteins/data_manuscript_2026_april/dn_ds_per_gene_human_mouse")
)
dnds <- dnds[, c("gene", "dNdS")]
colnames(dnds)[1] <- "Pheno"
dnds <- dnds[dnds$Pheno %in% pheno_to_use,]

moran <- read_excel("~/combined_Results_morgan_tabea.xlsx")
colnames(moran)[1] <- "Pheno"

moran <- moran[moran$Pheno %in% pheno_to_use,]

gdi <- read_excel("~/Downloads/pnas.1518646112.sd01.xlsx")
colnames(gdi)[1] <- "Pheno"
colnames(gdi)[2] <- "gdi"
gdi <- gdi[gdi$Pheno %in% pheno_to_use,]


# Extract gene symbol after the "_"

omega <- fread("~/selection_proteins/scripts_figure_may_26/data/GeneTable_w_thibault.tsv")

# Extract gene symbol after the "_"
omega$Pheno <- sub(".*_", "", omega$ENSG)
omega <- as.data.frame(omega)

omega <- omega[, c("Pheno", "ωA_phy")]
colnames(omega) <- c("Pheno", "wA")

omega <- omega[omega$Pheno %in% pheno_to_use,]


# ============================================================
# Build variable lists
# ============================================================

alpha_vars <- list(
  ldpred2 = prep_var(
    ldpred2_to_use,
    "alpha_lpred2",
    "ldpred2",
    "LDpred2"
  ),
  
  burden = prep_var(
    res_all_to_use,
    "alpha_burden",
    "burden",
    "Burden"
  ),
  
  health_score = prep_var(
    prot_lev_compscore_to_use,
    "alpha_health_score",
    "health_score",
    "Health score"
  ),
  
  long = prep_var(
    prot_lev_long_to_use,
    "alpha_long",
    "long",
    "Longevity"
  ),
  
  nchild = prep_var(
    prot_lev_nchild_to_use,
    "alpha_nchild",
    "nchild",
    "N.children"
  )
)

metric_vars <- list(
  h2_est_1 = prep_var(
    ldpred2_prot_no_cis,
    "h2_est_1",
    "h2_est_1",
    "h2"
  ),
  
  pLI = prep_var(
    gnomad_to_use,
    "pLI",
    "pLI",
    "pLI"
  ),
  
  LOEUF = prep_var(
    gnomad_to_use,
    "LOEUF",
    "LOEUF",
    "LOEUF"
  ),
  
  meta_transcript_ev_rate = prep_var(
    ev_tr,
    "meta_transcript_ev_rate",
    "meta_transcript_ev_rate",
    "Transcript evolutionary rate"
  ),
  
  `Lasso R2` = prep_var(
    coef_reg_u,
    "Lasso R2",
    "Lasso R2",
    "Lasso R2"
  ),
  
  dnds = prep_var(
    dnds,
    "dNdS",
    "dnds",
    "dN/dS"
  ),
  wA = prep_var(
    omega,
    "wA",
    "wA",
    "wA"
  ),
  
  gdi = prep_var(
    gdi,
    "gdi",
    "gdi",
    "GDI"
  ),
  
  moran = prep_var(
    moran,
    "moran",
    "moran",
    "Moran"
  )
)

expected_direction <- list(
  pLI = -1,
  LOEUF = 1,
  meta_transcript_ev_rate = 1,
  `Lasso R2` = -1,
  dnds = 1,
  gdi = 1,
  moran = -1, 
  wA = 1
)


# ============================================================
# Compute full correlation matrix
# ============================================================

all_vars <- c(alpha_vars, metric_vars)
var_order <- names(all_vars)
alpha_names <- names(alpha_vars)

res <- rbindlist(lapply(seq_along(var_order), function(i) {
  rbindlist(lapply(seq_along(var_order), function(j) {
    
    v_i <- var_order[i]
    v_j <- var_order[j]
    
    dt_i <- all_vars[[v_i]]
    dt_j <- all_vars[[v_j]]
    
    exp_sign <- NULL
    
    if (v_i %in% alpha_names && v_j %in% names(expected_direction)) {
      exp_sign <- expected_direction[[v_j]]
    } else if (v_j %in% alpha_names && v_i %in% names(expected_direction)) {
      exp_sign <- expected_direction[[v_i]]
    }
    
    s <- pairwise_spearman_dt(dt_i, v_i, dt_j, v_j, exp_sign)
    
    data.table(
      var_i = v_i,
      var_j = v_j,
      rho = s$rho,
      p = s$p,
      n = s$n
    )
  }))
}))

# ============================================================
# Pretty labels
# ============================================================

pretty_names <- c(
  ldpred2 = "LDpred2",
  burden = "Burden",
  nchild = "N.children",
  long = "Longevity",
  health_score = "Health Score",
  h2_est_1 = "h2",
  pLI = "pLI",
  LOEUF = "LOEUF",
  meta_transcript_ev_rate = "Transcript evol. rate",
  `Lasso R2` = "Lasso R2",
  dnds = "dN/dS",
  wA = "omega[A]", 
  gdi = "GDI",
  moran = "Moran's I"
)

plot_order <- pretty_names[var_order]

res[, var_i_plot := pretty_names[as.character(var_i)]]
res[, var_j_plot := pretty_names[as.character(var_j)]]

res[, var_i_plot := factor(var_i_plot, levels = plot_order)]
res[, var_j_plot := factor(var_j_plot, levels = plot_order)]

res[, sig := !is.na(p) & p < 0.05]

res[, label := fifelse(
  is.na(rho),
  "",
  sprintf("%.2f\np=%.1e", rho, p)
)]

parse_axis_labels <- function(x) {
  txt <- vapply(x, function(z) {
    if (z == "omega[A]") {
      "omega[A]"
    } else {
      encodeString(z, quote = "'")
    }
  }, character(1))
  
  parse(text = txt)
}
# ============================================================
# Manuscript heatmap
# ============================================================

p_all <- ggplot(res, aes(x = var_j_plot, y = var_i_plot, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.75) +
  geom_tile(
    data = res[sig == TRUE],
    aes(x = var_j_plot, y = var_i_plot),
    fill = NA,
    color = "black",
    linewidth = 0.75,
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = label),
    size = 2.2,
    lineheight = 0.9,
    color = "black"
  )+
  scale_fill_gradient2(
    low = "#B23A48",
    mid = "white",
    high =  "#3B6EA8",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey92",
    name = "Spearman\nrho"
  ) +
  scale_x_discrete(labels = parse_axis_labels) +
  scale_y_discrete(labels = parse_axis_labels) + 
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL,
    title = "Concordance between inferred selection estimates and independent metrics",
    subtitle = "No p-value filtering was applied to selection estimates. Black boxes indicate Spearman's correlation p < 0.05."
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(4, 4, 4, 4)
  )

p_all

# ============================================================
# OVERLAP SAMPLE SIZE TABLE
# ============================================================

# Create overlap labels
res[, overlap_label := fifelse(
  is.na(n),
  "",
  as.character(n)
)]

# Dynamic label color: white for dark tiles, black for light tiles
threshold <- median(res$n, na.rm = TRUE)

res[, label_color := fifelse(
  is.na(n),
  "black",
  fifelse(n > threshold, "white", "black")
)]

# Plot
p_overlap <- ggplot(
  res,
  aes(
    x = var_j_plot,
    y = var_i_plot,
    fill = n
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.75
  ) +
  
  geom_text(
    aes(
      label = overlap_label,
      color = label_color
    ),
    size = 3,
    fontface = "bold"
  ) +
  
  scale_fill_gradient(
    low = "grey95",
    high = "grey30",
    na.value = "grey95",
    name = "Overlap\nN"
  ) +
  
  scale_color_identity() +
  scale_x_discrete(labels = parse_axis_labels) +
  scale_y_discrete(labels = parse_axis_labels) +
  coord_fixed() +
  
  labs(
    x = NULL,
    y = NULL,
    title = "Protein overlap used for pairwise correlations",
    subtitle = "Number of overlapping proteins shared between each pair of analyses \n(no p value filtering applied on selection estimates)"
  ) +
  
  theme_classic(base_size = 11, base_family = "Arial") +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 10)
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(8, 8, 8, 8)
  )

p_overlap


# ============================================================
# SAVE
# ============================================================



