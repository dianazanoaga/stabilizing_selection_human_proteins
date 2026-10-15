# =========================
# Load packages
# =========================

library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)
library(gridExtra)
library(grid)
library(tidyr)
library(patchwork)

rm(list = ls())

# Load data 

ldpred2_prot_no_cis <- read_excel("~/selection_proteins/data_last/process/ldpred2-auto_proteins_coefficients_no_cis.xlsx")
ldpred2_prot_no_cis <- ldpred2_prot_no_cis[!is.na(ldpred2_prot_no_cis$alpha_est_1), ]
colnames(ldpred2_prot_no_cis)[1] <- "Pheno"

colnames(ldpred2_prot_no_cis)[8] <- "alpha_lpred2"
colnames(ldpred2_prot_no_cis)[13] <- "pval_lpred2"
colnames(ldpred2_prot_no_cis)[14] <- "pval_lpred2_adj"

ldpred2_prot_no_cis$dist <- rep(1, nrow(ldpred2_prot_no_cis))

se <- (ldpred2_prot_no_cis$p_est_3 - ldpred2_prot_no_cis$p_est_2) / (2 * 1.96)
z  <- (ldpred2_prot_no_cis$p_est_1 - 0) / se
ldpred2_prot_no_cis$pval_p <- 2 * (1 - pnorm(abs(z)))

se <- (ldpred2_prot_no_cis$h2_est_3 - ldpred2_prot_no_cis$h2_est_2) / (2 * 1.96)
z  <- (ldpred2_prot_no_cis$h2_est_1 - 0) / se
ldpred2_prot_no_cis$pval_h2 <- 2 * (1 - pnorm(abs(z)))

ldpred2_prot_no_cis <- ldpred2_prot_no_cis[!duplicated(ldpred2_prot_no_cis$Pheno), ]

h2_lookup <- as.data.frame(ldpred2_prot_no_cis)[, c("Pheno", "h2_est_1")]

# QC filter 
ldpred2_prot_no_cis <- ldpred2_prot_no_cis[
  ldpred2_prot_no_cis$h2_est_1 > 0.1 &
    ldpred2_prot_no_cis$p_est_1 > 0.00025 &
    ldpred2_prot_no_cis$pval_h2 < 0.05 &
    ldpred2_prot_no_cis$pval_p < 0.05, ]

### gene burden 
res_all <- fread("~/alpha_sigmaG_results_all_proteins_0.001_correction_optim.tsv")

colnames(res_all)[3] <- "alpha_burden"
colnames(res_all)[6] <- "pval_burden"
colnames(res_all)[11] <- "pval_burden_adj"

res_all <- res_all[res_all$pval_burden < 1, ]
res_all$dist <- rep(1, nrow(res_all))

### prot_lev

prot_lev_nchild     <- fread("~/selection_proteins/data_last/process/p_nchild_edu.tsv")
prot_lev_long       <- fread("~/selection_proteins/data_last/process/p_longevity.tsv")
prot_lev_compscore  <- fread("~/selection_proteins/data_last/process/p_health_score.tsv")

colnames(prot_lev_nchild)[11]    <- "alpha_nchild"
colnames(prot_lev_long)[11]      <- "alpha_long"
colnames(prot_lev_compscore)[11] <- "alpha_health_score"

colnames(prot_lev_nchild)[14]    <- "pval_nchild"
colnames(prot_lev_long)[14]      <- "pval_long"
colnames(prot_lev_compscore)[14] <- "pval_health_score"

colnames(prot_lev_nchild)[15]    <- "pval_nchild_adj"

colnames(prot_lev_long)[15]      <- "pval_long_adj"
colnames(prot_lev_compscore)[15] <- "pval_health_score_adj"

colnames(prot_lev_nchild)[20]    <- "dist"
colnames(prot_lev_long)[20]      <- "dist"
colnames(prot_lev_compscore)[20] <- "dist"

prot_lev_nchild    <- prot_lev_nchild[, -1]
prot_lev_long      <- prot_lev_long[, -1]
prot_lev_compscore <- prot_lev_compscore[, -1]

colnames(prot_lev_nchild)[1]    <- "Pheno"
colnames(prot_lev_long)[1]      <- "Pheno"
colnames(prot_lev_compscore)[1] <- "Pheno"

prot_lev_nchild <- merge(prot_lev_nchild, h2_lookup, by = "Pheno")
prot_lev_nchild <- prot_lev_nchild[prot_lev_nchild$dist == 1 & prot_lev_nchild$h2_est_1 > 0.1, ]

prot_lev_long <- merge(prot_lev_long, h2_lookup, by = "Pheno")
prot_lev_long <- prot_lev_long[prot_lev_long$dist == 1 & prot_lev_long$h2_est_1 > 0.1, ]

prot_lev_compscore <- merge(prot_lev_compscore, h2_lookup, by = "Pheno")
prot_lev_compscore <- prot_lev_compscore[prot_lev_compscore$dist == 1 & prot_lev_compscore$h2_est_1 > 0.1, ]

all_dt <- list(ldpred2_prot_no_cis, res_all, prot_lev_nchild, prot_lev_long, prot_lev_compscore)


# Names for the 5 methods
method_names <- c(
  "LDPred2",
  "Burden",
  "N.children (dist ==1)",
  "Longevity (dist ==1)",
  "Health_Score (dist ==1)"
)


# For each dataset: which alpha column + which pval column to use
spec <- data.table(

  method = method_names,

  alpha_col = c("alpha_lpred2", "alpha_burden", "alpha_nchild", "alpha_long", "alpha_health_score"),
  pval_col  = c("pval_lpred2", "pval_burden", "pval_nchild", "pval_long", "pval_health_score"),
  dist_col = rep("dist", 5)
)

spec <- data.table(

  method = method_names,

  alpha_col = c("alpha_lpred2", "alpha_burden", "alpha_nchild", "alpha_long", "alpha_health_score"),
  pval_col  = c("pval_lpred2_adj", "pval_burden_adj", "pval_nchild_adj", "pval_long_adj", "pval_health_score_adj"),
  dist_col = rep("dist", 5)
)

# helper that returns the 3 stats for one data.table/data.frame
summ_one <- function(dt, alpha_col, pval_col, dist_col, pthr = 0.05) {

  dt <- as.data.table(dt)

  a <- dt[[alpha_col]]
  p <- dt[[pval_col]]
  d <- dt[[dist_col]]

  ok <- !is.na(a) & !is.na(p) & !is.na(d)
  a <- a[ok]; p <- p[ok]; d <- d[ok]

  neg_sig <- sum(a < 0 & p < pthr & d == 1)
  pos_sig <- sum(a > 0 & p < pthr & d == 1)

  pct_neg <- if ((neg_sig + pos_sig) == 0) NA_real_ else 100 * neg_sig / (neg_sig + pos_sig)

  list(neg_sig = neg_sig, pos_sig = pos_sig, pct_neg = pct_neg)
}

# build the 3-row x 5-column table
out <- rbindlist(lapply(seq_along(all_dt), function(i) {
  s <- summ_one(all_dt[[i]], spec$alpha_col[i], spec$pval_col[i], spec$dist_col[i], pthr = 0.05)
  data.table(method = spec$method[i], neg_sig = s$neg_sig, pos_sig = s$pos_sig, pct_neg = s$pct_neg)
}))


tab <- as.data.frame(t(out))
tab$"category" <- rownames(tab)

colnames(tab) <- tab[1, ]
tab <- tab[-1, ]

tab_df <- as.data.frame(tab)

tab_df <- tab_df[, c(1, 2, 5, 4, 3, 6)]


tt <- ttheme_default(
  core = list(fg_params = list(cex = 1)),
  colhead = list(fg_params = list(cex = 1))
)


g <- tableGrob(tab_df, rows = NULL, theme = tt)

summarize_alpha_counts <- function(df, alpha_col, pval_col,
                                   method_name = "method",
                                   filter_dist1 = FALSE,
                                   dist_col = "dist",
                                   p_thresh = 0.05) {
  df <- as.data.frame(df)

  if (filter_dist1) {
    if (!dist_col %in% colnames(df)) {
      stop(sprintf("Column '%s' not found, but filter_dist1=TRUE was set for %s.",
                   dist_col, method_name))
    }
    df <- df[df[[dist_col]] == 1, , drop = FALSE]
  }

  df <- df[!is.na(df[[alpha_col]]), , drop = FALSE]

  df$cat <- "Not significant / NA"
  sig <- !is.na(df[[pval_col]]) & df[[pval_col]] < p_thresh

  df$cat[sig & df[[alpha_col]] < 0] <- "coefficient < 0 (p<0.05)"
  df$cat[sig & df[[alpha_col]] > 0] <- "coefficient > 0 (p<0.05)"

  out <- df %>%
    dplyr::count(cat) %>%
    dplyr::mutate(method = method_name)

  out$cat <- factor(out$cat,
                    levels = c("coefficient < 0 (p<0.05)",
                               "coefficient > 0 (p<0.05)",
                               "Not significant / NA"))
  out
}


counts_ldpred2 <- summarize_alpha_counts(
  ldpred2_prot_no_cis,
  alpha_col = "alpha_lpred2",
  pval_col  = "pval_lpred2",
  method_name = "LDpred2",
  filter_dist1 = FALSE
)

counts_burden <- summarize_alpha_counts(
  res_all,
  alpha_col = "alpha_burden",
  pval_col  = "pval_burden",
  method_name = "Burden",
  filter_dist1 = FALSE
)

counts_obs_nchild <- summarize_alpha_counts(
  prot_lev_nchild,
  alpha_col = "alpha_nchild",
  pval_col  = "pval_nchild",
  method_name = "N.children (dist == 1)",
  filter_dist1 = TRUE,
  dist_col = "dist"
)

counts_obs_long <- summarize_alpha_counts(
  prot_lev_long,
  alpha_col = "alpha_long",
  pval_col  = "pval_long",
  method_name = "Longevity (dist == 1)",
  filter_dist1 = TRUE,
  dist_col = "dist"
)

counts_obs_health <- summarize_alpha_counts(
  prot_lev_compscore,
  alpha_col = "alpha_health_score",
  pval_col  = "pval_health_score",
  method_name = "Health_Score (dist == 1)",
  filter_dist1 = TRUE,
  dist_col = "dist"
)


###
all_counts <- dplyr::bind_rows(
  counts_ldpred2,
  counts_burden,
  counts_obs_nchild,
  counts_obs_long,
  counts_obs_health
)

# Optional: set panel order
all_counts$method <- factor(all_counts$method, levels = c(
  "LDpred2", "Burden", "Health_Score (dist == 1)", "Longevity (dist == 1)", "N.children (dist == 1)"
))



### NICEST PLOT ##########


# nicer labels
all_counts <- all_counts %>%
  dplyr::mutate(
    method_clean = dplyr::recode(as.character(method),
                          "LDpred2" = "LDpred2",
                          "Burden" = "Burden",
                          "Health_Score (dist == 1)" = "Health Score",
                          "Longevity (dist == 1)" = "Longevity",
                          "N.children (dist == 1)" = "N.children"
    ),
    method_clean = factor(method_clean,
                          levels = c("LDpred2", "Burden", "Health Score",
                                     "Longevity", "N.children")
    ),
    cat_clean = dplyr::recode(as.character(cat),
                       "coefficient < 0 (p<0.05)" = "Stabilizing",
                       "coefficient > 0 (p<0.05)" = "Disruptive",
                       "Not significant / NA" = "Not significant"
    ),
    cat_clean = factor(cat_clean,
                       levels = c("Stabilizing", "Disruptive", "Not significant")
    )
  )

manuscript_cols <- c(
  "Stabilizing" = "aquamarine4",
  "Disruptive" = "tomato3",
  "Not significant" = "#BDBDBD"
)


p_counts <- ggplot(all_counts, aes(x = cat_clean, y = n, fill = cat_clean)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = n),
            size = 3.2,
            fontface = "plain",
            color = "black",
            vjust = -0.3) +
  facet_wrap(~ method_clean, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = manuscript_cols) +
  labs(
    x = NULL,
    y = "Number of proteins",
    title = "Distribution and direction of inferred selection coefficients",
    subtitle = paste0(
      "Proteins classified by coefficient sign and nominal significance (p < 0.05)"
    )
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13.5, hjust = 0),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing.x = unit(0.8, "lines"), plot.subtitle = element_text(
      size = 11.5,
      margin = margin(b = 12)
    ),
  )


# =========================
# Pie data: significant only
# =========================

pie_data <- all_counts %>%
  dplyr::filter(cat_clean %in% c("Stabilizing", "Disruptive")) %>%
  dplyr::group_by(method_clean) %>%
  dplyr::mutate(
    pct = round(100 * n / sum(n), 1),
    label = paste0(pct, "%")
  ) %>%
  dplyr::ungroup()

# =========================
# Function to make one full pie
# =========================

make_pie <- function(method_name) {

  df <- pie_data %>%
    dplyr::filter(method_clean == method_name)

  ggplot(df, aes(x = "", y = n, fill = cat_clean)) +
    geom_col(
      width = 1,
      color = "white",
      linewidth = 0.6
    ) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4,
      fontface = "bold",
      color = "white"
    ) +
    scale_fill_manual(
      values = c(
        "Stabilizing" = "aquamarine4",
        "Disruptive" = "tomato3"
      )
    ) +
    theme_void(base_size = 11, base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 35, 0, 35)
    )
}

# =========================
# Make pies
# =========================

pie1 <- make_pie("LDpred2")
pie2 <- make_pie("Burden")
pie3 <- make_pie("Health Score")
pie4 <- make_pie("Longevity")
pie5 <- make_pie("N.children")

pies <- pie1 + pie2 + pie3 + pie4 + pie5 +
  plot_layout(ncol = 5)

# =========================
# Method-name row between pies and barplots
# =========================

method_labels <- data.frame(
  method_clean = factor(
    c("LDpred2", "Burden", "Health Score", "Longevity", "N.children"),
    levels = c("LDpred2", "Burden", "Health Score", "Longevity", "N.children")
  ),
  x = 1,
  y = 1
)

p_method_labels <- ggplot(method_labels, aes(x = x, y = y, label = method_clean)) +
  geom_text(
    fontface = "bold",
    size = 4,
    family = "Arial"
  ) +
  facet_wrap(~ method_clean, nrow = 1) +
  theme_void(base_size = 11, base_family = "Arial") +
  theme(
    strip.text = element_blank(),
    plot.margin = margin(0, 35, 0, 35),
    panel.spacing.x = unit(0.8, "lines")
  )

# =========================
# Clean barplot
# =========================

p_counts_clean <- p_counts +
  labs(
    title = NULL,
    subtitle = NULL
  ) +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 11,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 11,
      face = "bold"
    )
  )

# =========================
# Final combined plot (nominal p < 0.05)
# =========================

combined_plot <- (pies / p_method_labels / p_counts_clean) +
  plot_layout(
    heights = c(1.2, 0.18, 2.5)
  ) +
  plot_annotation(
    title = "Distribution and direction of inferred selection coefficients",
    subtitle = paste0(
      "Proteins classified by coefficient sign and nominal significance (p < 0.05).\n",
      "Pie charts show the percentage of significant proteins under stabilizing versus disruptive selection.\n"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5,
        margin = margin(b = 10)
      )
    )
  )

combined_plot


##
### with adj (FDR) ###

counts_ldpred2 <- summarize_alpha_counts(
  ldpred2_prot_no_cis,
  alpha_col = "alpha_lpred2",
  pval_col  = "pval_lpred2_adj",
  method_name = "LDpred2",
  filter_dist1 = FALSE
)

counts_burden <- summarize_alpha_counts(
  res_all,
  alpha_col = "alpha_burden",
  pval_col  = "pval_burden_adj",
  method_name = "Burden",
  filter_dist1 = FALSE
)

counts_obs_nchild <- summarize_alpha_counts(
  prot_lev_nchild,
  alpha_col = "alpha_nchild",
  pval_col  = "pval_nchild_adj",
  method_name = "N.children (dist == 1)",
  filter_dist1 = TRUE,
  dist_col = "dist"
)

counts_obs_long <- summarize_alpha_counts(
  prot_lev_long,
  alpha_col = "alpha_long",
  pval_col  = "pval_long_adj",
  method_name = "Longevity (dist == 1)",
  filter_dist1 = TRUE,
  dist_col = "dist"
)


counts_obs_health <- summarize_alpha_counts(
  prot_lev_compscore,
  alpha_col = "alpha_health_score",
  pval_col  = "pval_health_score_adj",
  method_name = "Health_Score (dist == 1)",
  filter_dist1 = TRUE,
  dist_col = "dist"
)
###
all_counts <- dplyr::bind_rows(
  counts_ldpred2,
  counts_burden,
  counts_obs_nchild,
  counts_obs_long,
  counts_obs_health
)

all_counts$method <- factor(all_counts$method, levels = c(
  "LDpred2", "Burden", "Health_Score (dist == 1)", "Longevity (dist == 1)", "N.children (dist == 1)"
))


### NICEST PLOT ##########


all_counts <- all_counts %>%
  dplyr::mutate(
    method_clean = dplyr::recode(as.character(method),
                          "LDpred2" = "LDpred2",
                          "Burden" = "Burden",
                          "Health_Score (dist == 1)" = "Health Score",
                          "Longevity (dist == 1)" = "Longevity",
                          "N.children (dist == 1)" = "N.children"
    ),
    method_clean = factor(method_clean,
                          levels = c("LDpred2", "Burden", "Health Score",
                                     "Longevity", "N.children")
    ),
    cat_clean = dplyr::recode(as.character(cat),
                       "coefficient < 0 (p<0.05)" = "Stabilizing",
                       "coefficient > 0 (p<0.05)" = "Disruptive",
                       "Not significant / NA" = "Not significant"
    ),
    cat_clean = factor(cat_clean,
                       levels = c("Stabilizing", "Disruptive", "Not significant")
    )
  )

manuscript_cols <- c(
  "Stabilizing" = "aquamarine4",
  "Disruptive" = "tomato3",
  "Not significant" = "#BDBDBD"
)


p_counts <- ggplot(all_counts, aes(x = cat_clean, y = n, fill = cat_clean)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = n),
            size = 3.2,
            fontface = "plain",
            color = "black",
            vjust = -0.3) +
  facet_wrap(~ method_clean, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = manuscript_cols) +
  labs(
    x = NULL,
    y = "Number of proteins",
    title = "Distribution and direction of inferred selection coefficients",
    subtitle = paste0(
      "Proteins classified by coefficient sign and FDR significance (p < 0.05)"
      
    )
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13.5, hjust = 0),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing.x = unit(0.8, "lines"), plot.subtitle = element_text(
      size = 11.5,
      margin = margin(b = 12)
    ),
  )


# =========================
# Pie data: significant only
# =========================

pie_data <- all_counts %>%
  dplyr::filter(cat_clean %in% c("Stabilizing", "Disruptive")) %>%
  dplyr::group_by(method_clean) %>%
  dplyr::mutate(
    pct = round(100 * n / sum(n), 1),
    label = paste0(pct, "%")
  ) %>%
  dplyr::ungroup()

make_pie <- function(method_name) {

  df <- pie_data %>%
    dplyr::filter(method_clean == method_name)

  ggplot(df, aes(x = "", y = n, fill = cat_clean)) +
    geom_col(
      width = 1,
      color = "white",
      linewidth = 0.6
    ) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4,
      fontface = "bold",
      color = "white"
    ) +
    scale_fill_manual(
      values = c(
        "Stabilizing" = "aquamarine4",
        "Disruptive" = "tomato3"
      )
    ) +
    theme_void(base_size = 11, base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 35, 0, 35)
    )
}

pie1 <- make_pie("LDpred2")
pie2 <- make_pie("Burden")
pie3 <- make_pie("Health Score")
pie4 <- make_pie("Longevity")
pie5 <- make_pie("N.children")

pies <- pie1 + pie2 + pie3 + pie4 + pie5 +
  plot_layout(ncol = 5)

method_labels <- data.frame(
  method_clean = factor(
    c("LDpred2", "Burden", "Health Score", "Longevity", "N.children"),
    levels = c("LDpred2", "Burden", "Health Score", "Longevity", "N.children")
  ),
  x = 1,
  y = 1
)

p_method_labels <- ggplot(method_labels, aes(x = x, y = y, label = method_clean)) +
  geom_text(
    fontface = "bold",
    size = 5,
    family = "Arial"
  ) +
  facet_wrap(~ method_clean, nrow = 1) +
  theme_void(base_size = 11, base_family = "Arial") +
  theme(
    strip.text = element_blank(),
    plot.margin = margin(0, 35, 0, 35),
    panel.spacing.x = unit(0.8, "lines")
  )

p_counts_clean <- p_counts +
  labs(
    title = NULL,
    subtitle = NULL
  ) +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 11,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 11,
      face = "bold"
    )
  )

# =========================
# Final combined plot (FDR)
# =========================

combined_plot_fdr <- (pies / p_method_labels / p_counts_clean) +
  plot_layout(
    heights = c(1.2, 0.18, 2.5)
  ) +
  plot_annotation(
    title = "Distribution and direction of inferred selection coefficients",
    subtitle = paste0(
      "Proteins classified by coefficient sign and FDR significance (p < 0.05).\n",
      "Pie charts show the percentage of significant proteins under stabilizing versus disruptive selection.\n"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5,
        margin = margin(b = 10)
      )
    )
  )

combined_plot_fdr


# ============================================================
# Extra quadratic-regression datasets (N.children no edu, Tea intake),
# restricted to dist == 1 & h2 (LDpred2) > 0.1, same as above
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

prot_noedu <- fread("~/selection_proteins/data_last/process/p_nchil_noedu.tsv")
prot_tea   <- fread("~/selection_proteins/data_last/process/p_tea_intake.tsv")

# Use quadratic coefficient + raw p-value + dist_1
prot_noedu$alpha_noedu <- prot_noedu$Coefficients_quad
prot_tea$alpha_tea     <- prot_tea$Coefficients_quad

prot_noedu$pval_noedu <- prot_noedu$PValues_quad
prot_tea$pval_tea     <- prot_tea$PValues_quad

prot_noedu$dist <- prot_noedu$dist_1
prot_tea$dist   <- prot_tea$dist_1

# Same Gene -> Pheno trick as above, to merge in LDpred2's h2_est_1
prot_noedu_pheno <- prot_noedu$Gene
prot_tea_pheno   <- prot_tea$Gene

prot_noedu$Pheno <- prot_noedu_pheno
prot_tea$Pheno   <- prot_tea_pheno

prot_noedu <- merge(prot_noedu, h2_lookup, by = "Pheno")
prot_noedu <- prot_noedu[prot_noedu$dist == 1 & prot_noedu$h2_est_1 > 0.1, ]

prot_tea <- merge(prot_tea, h2_lookup, by = "Pheno")
prot_tea <- prot_tea[prot_tea$dist == 1 & prot_tea$h2_est_1 > 0.1, ]

# =========================
# Helper function
# =========================

summarize_alpha_counts <- function(df, alpha_col, pval_col,
                                   method_name = "method",
                                   filter_dist1 = TRUE,
                                   dist_col = "dist",
                                   p_thresh = 0.05) {
  df <- as.data.frame(df)

  if (filter_dist1) {
    df <- df[df[[dist_col]] == 1, , drop = FALSE]
  }

  df <- df[!is.na(df[[alpha_col]]), , drop = FALSE]

  df$cat <- "Not significant / NA"
  sig <- !is.na(df[[pval_col]]) & df[[pval_col]] < p_thresh

  df$cat[sig & df[[alpha_col]] < 0] <- "coefficient < 0 (p<0.05)"
  df$cat[sig & df[[alpha_col]] > 0] <- "coefficient > 0 (p<0.05)"

  out <- df %>%
    dplyr::count(cat) %>%
    dplyr::mutate(method = method_name)

  out$cat <- factor(
    out$cat,
    levels = c(
      "coefficient < 0 (p<0.05)",
      "coefficient > 0 (p<0.05)",
      "Not significant / NA"
    )
  )

  return(out)
}

# =========================
# Summaries
# =========================

counts_noedu <- summarize_alpha_counts(
  prot_noedu,
  alpha_col = "alpha_noedu",
  pval_col  = "pval_noedu",
  method_name = "N.children no edu",
  filter_dist1 = TRUE
)

counts_tea <- summarize_alpha_counts(
  prot_tea,
  alpha_col = "alpha_tea",
  pval_col  = "pval_tea",
  method_name = "Tea intake",
  filter_dist1 = TRUE
)

all_counts <- dplyr::bind_rows(counts_noedu, counts_tea)

all_counts <- all_counts %>%
  dplyr::mutate(
    method_clean = factor(
      method,
      levels = c("N.children no edu", "Tea intake")
    ),
    cat_clean = dplyr::recode(
      as.character(cat),
      "coefficient < 0 (p<0.05)" = "Stabilizing",
      "coefficient > 0 (p<0.05)" = "Disruptive",
      "Not significant / NA" = "Not significant"
    ),
    cat_clean = factor(
      cat_clean,
      levels = c("Stabilizing",
                 "Disruptive",
                 "Not significant")
    )
  )

manuscript_cols <- c(
  "Stabilizing" = "aquamarine4",
  "Disruptive" = "tomato3",
  "Not significant" = "#BDBDBD"
)

# =========================
# Barplot
# =========================

p_counts <- ggplot(all_counts, aes(x = cat_clean, y = n, fill = cat_clean)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = n),
            size = 3.2,
            color = "black",
            vjust = -0.3) +
  facet_wrap(~ method_clean, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = manuscript_cols) +
  labs(
    x = NULL,
    y = "Number of proteins"
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11,
                               face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.title.y = element_text(size = 11, face = "bold")
  )

# =========================
# Pie charts
# =========================

pie_data <- all_counts %>%
  dplyr::filter(cat_clean %in% c("Stabilizing",
                          "Disruptive")) %>%
  dplyr::group_by(method_clean) %>%
  dplyr::mutate(
    pct = round(100 * n / sum(n), 1),
    label = paste0(pct, "%")
  ) %>%
  dplyr::ungroup()

make_pie <- function(method_name) {
  df <- pie_data %>%
    dplyr::filter(method_clean == method_name)

  ggplot(df, aes(x = "", y = n, fill = cat_clean)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4,
      fontface = "bold",
      color = "white"
    ) +
    scale_fill_manual(
      values = c(
        "Stabilizing" = "aquamarine4",
        "Disruptive" = "tomato3"
      )
    ) +
    theme_void(base_size = 11, base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 35, 0, 35)
    )
}

pie1 <- make_pie("N.children no edu")
pie2 <- make_pie("Tea intake")

pies <- pie1 + pie2 + plot_layout(ncol = 2)

# =========================
# Method labels
# =========================

method_labels <- data.frame(
  method_clean = factor(
    c("N.children no edu", "Tea intake"),
    levels = c("N.children no edu", "Tea intake")
  ),
  x = 1,
  y = 1
)

p_method_labels <- ggplot(method_labels, aes(x = x, y = y, label = method_clean)) +
  geom_text(fontface = "bold", size = 5, family = "Arial") +
  facet_wrap(~ method_clean, nrow = 1) +
  theme_void(base_size = 11, base_family = "Arial") +
  theme(
    strip.text = element_blank(),
    plot.margin = margin(0, 35, 0, 35)
  )

# =========================
# Final plot
# =========================

combined_plot_extra <- (pies / p_method_labels / p_counts) +
  plot_layout(heights = c(1.2, 0.18, 2.5)) +
  plot_annotation(
    title = "Distribution and direction of inferred selection coefficients",
    subtitle = paste0(
      "Proteins classified by quadratic coefficient sign and nominal significance (p < 0.05).\n"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5,
                                   margin = margin(b = 10))
    )
  )

combined_plot_extra
